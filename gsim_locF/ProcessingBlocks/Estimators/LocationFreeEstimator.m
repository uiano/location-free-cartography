classdef LocationFreeEstimator < Estimator
    %LocationFreeEstimator Location-Free cartography class.
    %   This class allows to train and query a kernel
    %   machine desgined to estimate the received power from
    %   the location features received at sensors.
    
    properties
        kernel % function handle defining the kernel
        regularizationParameter % lambda (ridge regression)
        Xenb % locations of the sources
        enableMissingData
        completionMethodSelect % 1 for modified EGM, and 2 for Manopt
        receiverSensitivity % missing features, can be changes   ... use it for threshold
        gamma;% parameter used in the kernel of missing features: low rank and missing data
        desiredRank
        regParEval %regularization parameter for completing missing features when evaluating the map
        evalOption % choose the option of determining missing features in the map evaluation
    end
    
    methods
        function [coefficientsOut, intercept,ues_many_misses,...
                avgNumberOfMissingFeat,combin_sources,orthBasis,meanFeature,featCovarianceMat,completedmeasurements] = train(obj, measurements,...
                powerAtLocationAllsources)
            % given a matrix including measurements (distances and
            % received power at each location), optimize the coefficients
            % of the kernel machine that estimates the power.
            n_sources = size( measurements,1)-1; % n_sources  symbolizes number of features from all sources
            n_ues=size(measurements, 2);
            %             featRank=4;% is the rank of the feature space
            
            
            
            
            if  obj.enableMissingData==true  % incorporates the missing feature case.
                real_n_sources=size(obj.Xenb,2);
                combin_sources=combnk(1:real_n_sources,2);
                if real_n_sources < 6
                    combin_sources=flipud(combin_sources);
                end
                identify_misses=(powerAtLocationAllsources<obj.receiverSensitivity);
                X=measurements(1:end-1,:);
                %                 if obj.receiverSensitivity==0
                %                     Q=orth(X);
                %                     completedmeasurements =Q(:,1:featRank)'*X;
                %                 else
                
                o=zeros(n_sources,n_ues);
                ues_many_misses=[];
                allNumberOfMissingFeat=zeros(1,n_ues);
                for ind_n=1:n_ues
                    miss_ind=find(identify_misses(:,ind_n)==1);
                    for ind_missing=1:length(miss_ind)
                        miss_indices_on_features=union(find(combin_sources(:,1)==miss_ind(ind_missing)),...
                            find(combin_sources(:,2)==miss_ind(ind_missing)));
                        X(miss_indices_on_features,ind_n)=NaN; %NaN symbolizes a missing features
                    end
                    o(:,ind_n)=double(~isnan(X(:,ind_n)));
                    numberOfMissingFeat=sum(isnan(X(:,ind_n)));
                    allNumberOfMissingFeat(ind_n)=numberOfMissingFeat;
                    if numberOfMissingFeat> size(combin_sources,1)-obj.desiredRank % Removes the measurements having a large number of missing features
                        % there must be at least 4 features available.
                        ues_many_misses=[ues_many_misses, ind_n];
                    end
                end
                avgNumberOfMissingFeat=mean(allNumberOfMissingFeat);
                X(:,ues_many_misses)=[];
                o(:,ues_many_misses)=[];
                X_nan=X;
                X(isnan(X))=0;
                
                if obj.completionMethodSelect==1
                    completedmeasurementsToProject =singularValProj(X,o,obj.desiredRank);
                elseif obj.completionMethodSelect==2
                    completedmeasurementsToProject=low_rank_matrix_completion(X,o,obj.desiredRank);
                else
                    opts = [];
                    [X_lma,Y_lma,~] = lmafit_mc_adp(X_nan,obj.desiredRank,opts);
                    completedmeasurementsToProject=X_lma*Y_lma;
                end
                
                %Project completed measurements
                orthBasis=orth(completedmeasurementsToProject);
                completedmeasurements =orthBasis'*(completedmeasurementsToProject); %
                featRank=size(orthBasis,2);
                if obj.evalOption==1
                    meanFeature=mean(completedmeasurementsToProject,2);
                    featCovarianceMat=cov(completedmeasurementsToProject');
                else
                    Coeff_B=orthBasis'*completedmeasurementsToProject;
                    meanFeature=mean(Coeff_B,2);
                    featCovarianceMat=cov(Coeff_B');
                end
                
                %                 end
                n_ues=size(completedmeasurements,2);
                measurements(:,ues_many_misses)=[];
            else
                ues_many_misses=[];avgNumberOfMissingFeat=0;combin_sources=0;orthBasis=0;meanFeature=0;featCovarianceMat=0;completedmeasurements=0;
            end
            
            Ke1=zeros(n_ues,n_ues); % Kernel matrix
            for i = 1:n_ues
                for k=1:n_ues
                    
                    if  (obj.enableMissingData==true) && (obj.receiverSensitivity~=0)
                        x_input_to_kern=zeros(featRank,1);
                        y_input_to_kern=zeros(featRank,1);
                        for ind_s_1=1:featRank
                            x_input_to_kern(ind_s_1)=completedmeasurements(ind_s_1,i);
                            y_input_to_kern(ind_s_1)=completedmeasurements(ind_s_1,k);
                        end
                        Ke1(i,k) = feval(obj.kernel, x_input_to_kern, y_input_to_kern);
                        
                    else
                        x_input_to_kern=zeros(n_sources,1);
                        y_input_to_kern=zeros(n_sources,1);
                        for ind_s_1=1:n_sources
                            x_input_to_kern(ind_s_1)=measurements(ind_s_1,i);
                            y_input_to_kern(ind_s_1)=measurements(ind_s_1,k);
                        end
                        Ke1(i,k) = feval(obj.kernel, x_input_to_kern, y_input_to_kern);
                    end
                end
            end
            zeroMeanMeasurements=measurements;
            intercept = mean(measurements(end,:));
            zeroMeanMeasurements(end,:)=measurements(end,:)-intercept;
            coefficientsOut=(Ke1+((n_ues*obj.regularizationParameter)*eye(n_ues)))\zeroMeanMeasurements(end,:)';
        end
        
       
        
        function [meanErrOnEvalFeat,predictedMeasurements] = estimateGivenEstimatedDistances(obj, coefficients, trainingMeasurements,completedmeasurements,...
                t_distances,allPointAllSourcePower,combin_sources,orthBasis, meanFeat,CovarMat)
            % Estimate received power using Location Free cartography.
      
            ed_dim = size(t_distances);
            featNum=size(combin_sources,1);
            n_ues=size(trainingMeasurements, 2);
            predictedMeasurements=zeros(ed_dim(1:2));
            errAllCompletedEvalFeat=zeros(ed_dim(1:2));
            dist_to_source_input_to_ker = zeros(size(trainingMeasurements, 1)-1,1);
            for kx = 1:ed_dim(1)
                for ky = 1:ed_dim(2)
                    
                    
                    if  obj.enableMissingData==true
                        n_ues=size(completedmeasurements,2);
                        identify_missesEval=(allPointAllSourcePower(:,kx, ky)<obj.receiverSensitivity);
                        evalFeat=permute(t_distances(kx, ky, :),[3 2 1]);
                        miss_ind=find(identify_missesEval==1);
                        for ind_missing=1:length(miss_ind)
                            miss_indices_on_features=union(find(combin_sources(:,1)==miss_ind(ind_missing)),...
                                find(combin_sources(:,2)==miss_ind(ind_missing)));
                            evalFeat(miss_indices_on_features)=NaN; %NaN symbolizes a missing features
                        end
                        numberOfMissingFeatAtEval=sum(isnan(evalFeat));
                        if ( numberOfMissingFeatAtEval >featNum -obj.desiredRank) || (size(orthBasis,2) < obj.desiredRank)
                            predictedMeasurements(kx,ky)=mean(trainingMeasurements(end,:));
                            continue
                        end
                        
                        indOfObserved_Feat=find(~isnan(evalFeat));
                        obervedFeatNum=length(indOfObserved_Feat);
                        rowSelectMat=zeros(obervedFeatNum,featNum);
                        evalFeat(isnan(evalFeat))=0;
                        for ind_obsFeat=1:obervedFeatNum
                            rowSelectMat(ind_obsFeat,indOfObserved_Feat(ind_obsFeat))=1;
                        end
                        if obj.evalOption==1
                            weights=(orthBasis'*rowSelectMat'*rowSelectMat*orthBasis+obj.regParEval*orthBasis'*orthBasis)\...%(CovarMat\
                                (orthBasis'*rowSelectMat'*rowSelectMat*evalFeat+obj.regParEval*orthBasis'*meanFeat);
                        else
                            weights=(orthBasis'*rowSelectMat'*rowSelectMat*orthBasis+obj.regParEval*inv(CovarMat))\...
                                (orthBasis'*rowSelectMat'*rowSelectMat*evalFeat+obj.regParEval*(CovarMat\meanFeat));
                        end
                        dist_to_source_input_to_ker=orthBasis*weights;
                        %                         dist_to_source_input_to_ker(:)=t_distances(kx, ky, :);
                        errAllCompletedEvalFeat(kx, ky)=norm(weights-orthBasis'*permute(t_distances(kx, ky, :),[3 2 1]));
                    else
                        errAllCompletedEvalFeat(kx, ky)=0;
                        
                        dist_to_source_input_to_ker(:)=t_distances(kx, ky, :);
                        
                    end
                    
                    row_kernel=zeros(1,n_ues);
                    for ind_nues=1:n_ues
                        if  obj.enableMissingData==true
                            row_kernel(ind_nues)=feval(obj.kernel, weights, completedmeasurements(:,ind_nues)); % here the weights=orthBasis'*dist_to_source_input_to_ker
                            %  =orthBasis'*orthBasis*weights=weights;
                        else
                            row_kernel(ind_nues)=feval(obj.kernel, dist_to_source_input_to_ker, trainingMeasurements(1:end-1,ind_nues));
                        end
                    end
                    
                    
                    predictedMeasurements(kx,ky)=row_kernel*coefficients+mean(trainingMeasurements(end,:));
                    
                end
            end
            meanErrOnEvalFeat=mean(mean(errAllCompletedEvalFeat));
        end
        
    end %methods
end

