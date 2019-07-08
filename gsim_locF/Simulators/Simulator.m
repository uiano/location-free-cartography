classdef Simulator
    % Wireless Cartography Simulator.
    % Contains the generator, sampler, and estimators corresponding
    % to a single experiment with the purpose of calculating the
    % estimation error of the location-based and location-free
    % cartography estimators.
    
    properties
        generator MapGenerator = PrecomputedMapGenerator();
        sampler SpectrumMapSampler
        locFreeEstimator LocationFreeEstimator
        locBasedEstimator LocationBasedEstimator
        locEstimator
        featureExtractor;
        
        simLocBased = true % Only if set to true will evaluate the error in the location based estimation.
        features  %(1,:) vector containing the indices of the features used to estimate the C2M.
        tdoaEnable=true;
        PCAenable
        PCAdimension
    end
    
    methods
        
        function [SqErr_locFree, SqErr_locBased, meanErrOnEvalFeat, locFreeMapEstimate,UEFeatures,locBasedEstimate, UELocations,avgNumberOfMissingFeat] = ... % SqErr_sampleAverage
                simulate(obj, N_ue, allPointEstimatedLocation, power_out,trueMap,...
                evaluationGrid_x, evaluationGrid_y, source_loc, allLocFreeFeatures, estimatedPilotSignals)
            
            % Generate the map
            %             [trueMap, evaluationGrid_x, evaluationGrid_y, source_loc, estimatedDistances,estimatedCentersOfMass] = obj.generator.generateMap();
            
            % Sample the map
            obj.sampler.map = trueMap;
            obj.sampler.powerAllSources=power_out;
            
            UELocations=zeros(size(source_loc,1),N_ue);
            savedUELocationIndices=zeros(1,N_ue);
            for ind_ue=1:N_ue
                n_ueOutNearF=0;
                while (n_ueOutNearF < size(source_loc,2))
                    UELocationIndex = randi(numel(trueMap));
                    while ismember(UELocationIndex,savedUELocationIndices)
                        UELocationIndex = randi(numel(trueMap));
                    end
                    
                    UELocations_check =[evaluationGrid_x(UELocationIndex); evaluationGrid_y(UELocationIndex)];
                    ue_to_source_dist=pdist2(UELocations_check', source_loc');
                    n_ueOutNearF=sum(ue_to_source_dist > obj.generator.refDistance);
                end
                savedUELocationIndices(ind_ue)=UELocationIndex;
                UELocations(:,ind_ue)=UELocations_check;
            end
            
            obj.sampler.UELocations = UELocations;
            obj.sampler.Xenb=source_loc;
            obj.sampler.evalGrid_x = evaluationGrid_x;
            obj.sampler.evalGrid_y = evaluationGrid_y;
            
            
            %         if obj.FeatureInd==1
            %             estimatedCentersOfMassToConsider=estimatedPilotSignals(:,:, obj.features);
            
            % Location Free: (simulate measurements taken by UEs & train and test estimator
            % remove this line since source locations are not needed!
            [ ~,measurement_pilots,powerAtLocationAllsources, measurements_power]=obj.sampler.sampleGivenPilotSignals(estimatedPilotSignals);
            extractedLocfreeFeatures=obj.featureExtractor.locFreeExtract(measurement_pilots);
            extractedLocfreeFeaturesToConsider=extractedLocfreeFeatures(obj.features,:);
            
            
            obj.locFreeEstimator.Xenb = source_loc;
            
            % Location Free evaluation: load (previously saved) estimated fe and
            % test
%             loadOut=load (allLocFreeFeaturesFileName);
%             allLocFreeFeatures=loadOut.(allLocFreeFeaturesVariableName);
            evalGridSize=sqrt(size(allLocFreeFeatures,2));
            
            estimatedCentersOfMassSelected=allLocFreeFeatures(obj.features,:);
            eval_featNum=size(estimatedCentersOfMassSelected,1);
            estimatedCentersOfMassToConsider=reshape(estimatedCentersOfMassSelected',[evalGridSize,evalGridSize,eval_featNum]);
            
            if obj.PCAenable==true
                zeroMeanExtractedLocfreeFeaturesToConsider=extractedLocfreeFeaturesToConsider;
                PCAEestimatedCentersOfMass=obj.pcaComputation(zeroMeanExtractedLocfreeFeaturesToConsider);
                measurementsLF= [PCAEestimatedCentersOfMass;measurements_power];
                UEFeatures=measurementsLF(1:end-1,:);
                [locFreeCoefficients,~,~,...
                    avgNumberOfMissingFeat,combin_sources,orthBasis,meanFeature,...
                    featCovMat,completedmeasurements] = obj.locFreeEstimator.train(measurementsLF,powerAtLocationAllsources);
                
                % Evaluation with PCA
                %                 estimatedCentersOfMassSelected=estimatedCentersOfMassSelected;
                eval_featPCA=obj.pcaComputation(estimatedCentersOfMassSelected);
                eval_newFeatNum=size(eval_featPCA,1);
                evalGrid_PCA=reshape(eval_featPCA',[evalGridSize,evalGridSize,eval_newFeatNum]);
                
                
                [meanErrOnEvalFeat,locFreeMapEstimate] = obj.locFreeEstimator.estimateGivenEstimatedDistances(...
                    locFreeCoefficients,[extractedLocfreeFeaturesToConsider;measurements_power],completedmeasurements,...
                    estimatedCentersOfMassToConsider,power_out,combin_sources,orthBasis,mean(mean(trueMap)), meanFeature,featCovMat);
                
            else
                measurementsLF=[extractedLocfreeFeaturesToConsider;measurements_power];
                UEFeatures=measurementsLF(1:end-1,:);
                [locFreeCoefficients,~,ues_many_misses,avgNumberOfMissingFeat,combin_sources,orthBasis,meanFeature,...
                    featCovMat,completedmeasurements]= obj.locFreeEstimator.train(measurementsLF,powerAtLocationAllsources);
                
                [meanErrOnEvalFeat,locFreeMapEstimate] = obj.locFreeEstimator.estimateGivenEstimatedDistances(... % sugg: change method name to: evaluateEstimatedMap[given features]
                    locFreeCoefficients, measurementsLF,completedmeasurements,...
                    estimatedCentersOfMassToConsider,power_out,combin_sources,orthBasis,mean(mean(trueMap)), meanFeature,featCovMat);
            end
            
            
            SqErr_locFree = obj.compute_SqErr(trueMap, locFreeMapEstimate );
            SqErr_sampleAverage = obj.compute_SqErr(trueMap, ones(size(trueMap))*mean(measurementsLF(end,:)));
            
            
            if obj.simLocBased==true
                extractedLocBasedFeatures=obj.featureExtractor.locBasedExtract(measurement_pilots);
                measurementsLB=[extractedLocBasedFeatures;measurements_power];
                % Location Based: train estimator
                obj.locBasedEstimator.Xenb = source_loc;
                obj.locEstimator.Xenb = source_loc;
                if obj.tdoaEnable==false
                    estimatedLocation=obj.locBasedEstimator.estimateLocationFromDistances (measurementsLB);
                else
                    estimatedLocation=obj.locEstimator.estimateLocationIRWSRDLS(measurementsLB);
                end
                [locBasedCoefficients, locBasedIntercept] = obj.locBasedEstimator.train(estimatedLocation,measurementsLB);
                
                
                % Location Based: load (previously saved) estimated locations and
                % test
%                 loadOut=load (estimatedLocationsFileName);
%                 allPointEstimatedLocation=loadOut.(estimatedLocationsVariableName);
                estGrid_x=reshape(allPointEstimatedLocation(1,:),size(evaluationGrid_x));
                estGrid_y=reshape(allPointEstimatedLocation(2,:),size(evaluationGrid_y));
                
                locBasedEstimate = obj.locBasedEstimator.estimate(locBasedCoefficients, estimatedLocation,estGrid_x,estGrid_y, locBasedIntercept); %evaluateEstimateMap
                
                % Compare the estimates with the true map
                SqErr_locBased= obj.compute_SqErr(trueMap, locBasedEstimate);
            else
                locBasedEstimate = 0;
                SqErr_locBased = 0;
            end
            
        end
        
        function [allPointEstimatedLocation,allLocFreeFeatures,allLocBasedFeatures] = precomputeEstimatedLocations(obj)
            [power_out,trueMap, evaluationGrid_x, evaluationGrid_y, source_loc, estimatedDistances, estimatedPilotSignals] = obj.generator.generateMap();
            
            allPointSampler=obj.sampler;
            allPointSampler.map = trueMap;
            allPointSampler.powerAllSources=power_out;
            allPointSampler.Xenb=source_loc;
            allPointSampler.evalGrid_x = evaluationGrid_x;
            allPointSampler.evalGrid_y = evaluationGrid_y;
            allPointSampler.UELocations=[flipud(evaluationGrid_x(:)'); flipud(evaluationGrid_y(:)')]; %every point location
            [allnoiseFreePilots,~,~,allPower_measurements]=allPointSampler.sampleGivenPilotSignals(estimatedPilotSignals);
            %           estimatedDistancesTest=sign(estimatedDistances).* sqrt(estimatedDistances.^2-(obj.generator.zt-obj.generator.zr)^.2);
            % we use the difference in height between the
            % eNBs and UEs to calculate the (2D) distance over the ground
            % (norm of the location difference projected over the z=0 plane)
            % instead of the (3D) distance
            if obj.simLocBased==true
                allLocBasedFeatures= obj.featureExtractor.locBasedExtract(allnoiseFreePilots);
                if obj.generator.estimateLocFreeSpace==1
                    estimatedDistances=permute(estimatedDistances,[3 1 2]);
                    allPointMeasurements= [reshape(estimatedDistances,[size(estimatedDistances,1)
                        size(estimatedDistances,2)*size(estimatedDistances,3)]');allPower_measurements];
                else
                    allPointMeasurements= [obj.featureExtractor.locBasedExtract(allnoiseFreePilots);allPower_measurements];
                end
                allLocFreeFeatures=obj.featureExtractor.locFreeExtract(allnoiseFreePilots);
                
                obj.locBasedEstimator.Xenb = source_loc;
                obj.locEstimator.Xenb = source_loc;
                if obj.tdoaEnable==false
                    allPointEstimatedLocation=obj.locBasedEstimator.estimateLocationFromDistances(allPointMeasurements);
                else
                    allPointEstimatedLocation=obj.locEstimator.estimateLocationIRWSRDLS(allPointMeasurements);
                    %              allPointEstimatedLocation=allPointSampler.UELocations;
                end
            else
                
                allLocFreeFeatures=obj.featureExtractor.locFreeExtract(allnoiseFreePilots);
                allLocBasedFeatures=0;
                allPointEstimatedLocation=0;
            end
        end
        
        
        function [U1,newFeat]=pcaComputationFirst(obj,centerOfMass)
            Option=1;
            if Option==1
                meanofFeat=mean(centerOfMass,2);
                centeredFeat=centerOfMass-meanofFeat;
                [U,~,~] = svd(centeredFeat);
                U1=U(:,1:obj.PCAdimension);
                newFeat=U1'*centeredFeat;
            else
                nFeatures=size(centerOfMass,1); %featSize=size(centerOfMass,2);
                %             meanCenterOfMass=mean(centerOfMass,3);
                CovFeat=zeros(nFeatures);
                for ind_feat1=1:nFeatures
                    for ind_feat2=1:nFeatures
                        tempCov=cov(centerOfMass(ind_feat1,:),....
                            centerOfMass(ind_feat2,:));
                        CovFeat(ind_feat1,ind_feat2)=tempCov(1,2);
                    end
                end
                [Eigvector,~]=eig(CovFeat);
                flipped_Eigvector=fliplr(Eigvector);
                pcaSpace=flipped_Eigvector(:,1:obj.PCAdimension); % space corresponding to the... largest eigen value
                newFeat=pcaSpace'*(centerOfMass-mean(centerOfMass,2));
            end
        end
        
        function newFeat=pcaComputation(obj,centerOfMass)
            Option=1;
            if Option==1
                meanofFeat=mean(centerOfMass,2);
                centeredFeat=centerOfMass-meanofFeat;
                [U,~,~] = svd(centeredFeat);
                U1=U(:,1:obj.PCAdimension);
                newFeat=U1'*centeredFeat;
            else
                nFeatures=size(centerOfMass,1); %featSize=size(centerOfMass,2);
                %             meanCenterOfMass=mean(centerOfMass,3);
                CovFeat=zeros(nFeatures);
                for ind_feat1=1:nFeatures
                    for ind_feat2=1:nFeatures
                        tempCov=cov(centerOfMass(ind_feat1,:),....
                            centerOfMass(ind_feat2,:));
                        CovFeat(ind_feat1,ind_feat2)=tempCov(1,2);
                    end
                end
                [Eigvector,~]=eig(CovFeat);
                flipped_Eigvector=fliplr(Eigvector);
                pcaSpace=flipped_Eigvector(:,1:obj.PCAdimension); % space corresponding to the... largest eigen value
                newFeat=pcaSpace'*(centerOfMass-mean(centerOfMass,2));
            end
        end
        
    end %methods
    
    methods (Static)
        function SqErr = compute_SqErr(X,Y)
            SqErr=sum(sum((X-Y).^2, 'omitnan'));
        end
        
    end
    
end %classdef
