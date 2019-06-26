%
%  This file contains experiments for spectrum cartography. The
%  experiments 1001 up to 100... correspond to the paper  Location-free
%  spectrum cartography
%

classdef LocFCartogrExperiments < ExperimentFunctionSet
    
    properties
        
        b_loadData = 0 % set this property to 1 if you want the experiment 
        % methods to load data from the savedResults folder.
        
        % You may define here constants that you will use in the
        % experiments of this file
        
    end
    
    methods
        
        % This first experiment plots fig 1a and 1b
        function F=experiment_101(obj,niter)
            % Figure 1(a) in the paper
            
            b_toyExample = 0; %activate this flag if you want a shorter experiment
            st = dbstack;
            namestr = st.name;
            simulationNumber =str2double(namestr(35:end));
            results_filename = sprintf('./savedResults/dataTest%d.mat', simulationNumber);
            if obj.b_loadData && exist(results_filename, 'file')
                fprintf('%s found. This experiment has been run before but we are not checking the integrity of the saved files. Loading and plotting results...\n', results_filename);
                load(sprintf('./savedResults/dataTest%d.mat', simulationNumber))
            else
                fprintf('%s not found. Running your experiment...\n', results_filename);
                c = CartographySimulations;
                c.f = 800e6;
                c.pilotNoiseSTD =1e-5;
                c.powerNoiseSTD  =0.5;
                c.maxSamplesPerPilot=10;
                c.receiverBandwidth = 20e6;
                c.kernelSigmaLF =14.8;
                c.kernelSigmaLB =0.501;
                c.lambdaLF=1.91e-4;
                c.lambdaLB=3e-3;
                if b_toyExample
                    c.gridSize = [30 30];
                else
                    c.gridSize = [150 150];
                end
                c.numberOfSources=5;
                
                c.computeOnlyFeatures = true;
                
                selectedWalls=0;
                
                x_wall_file_name='./modelFiles/x_coord_4walls.mat';
                y_wall_file_name='./modelFiles/y_coord_4walls.mat';
                
                c.n_runs = niter;% 10
                c.N_ue =[50,150,250,350,450,550,750,1000,1500,2000,2400];
                if b_toyExample
                    c.N_ue = [20 30];
                end
                
                c.simLocBased=true;  % false disables LocB cartography
                c.testDiffSetOfFeatures=false;
                c.PCAenabler=false;
                c.considerMissingData=false;
                
                if  c.testDiffSetOfFeatures==true
                    c.featureCell=c.featureCombinations(10,4); % generate all combinations of features
                    c.N_feat_used=1:length(c.featureCell);
                    c.N_feat =c.N_feat_used;
                else
                    c.N_feat_used=10;
                    c.N_feat=1;
                end
                if  c.PCAenabler==true
                    c.N_feat =4; %  PCA dimension
                elseif c.considerMissingData==true
                    c.N_feat =-85:3:-41; % receiver sensitivy range in dBm
                end
                c.completionMethodSelect=1; % 1 for modified EGM or PGD, and else for Manopt
                c.desiredRank=4;
                c.regularParEval=5.42e-0;
                c.evalcoeffOption=2; % 1 for mean and covariance between completed features
                c.inParallel =true;
                c.forcePrecomputeDistances = false;
                
                [allPointEstimatedLocation,~, evaluationGrid_x,evaluationGrid_y, source_loc, ~, ~,...
                    ~,~,~,~,~,~,~,~,~,~]=c.run(selectedWalls, x_wall_file_name,...
                    y_wall_file_name,simulationNumber);
                
                save (sprintf('./savedResults/dataTest%d', simulationNumber))
                beep
                beep
            end
            
            %Plotting results
            smoothnessToMapLoc=zeros(c.gridSize(1),c.gridSize(2), 2);
            mycolorbar_caxis = [-10, 50];
            m_X = evaluationGrid_x(:,1);
            m_Y =evaluationGrid_y(1,:);
            source_X=source_loc(1,:);
            source_Y=source_loc(2,:);
            for coordNumber=1:2
                smoothnessToMapLoc(:,:,coordNumber)=(reshape((allPointEstimatedLocation(coordNumber,:)),[c.gridSize(1),c.gridSize(2)]));
                m_Z=smoothnessToMapLoc(:,:,coordNumber)';
                if coordNumber==1
                    ch_tit='$\hat{x}$';
                else
                    ch_tit='$\hat{y}$';
                end
                
                m_multiplot(1,coordNumber,1)= GFigure('m_X',m_X,'m_Y',m_Y,'m_Z',m_Z,'ch_plotType3D','imagesc',...
                    'ch_xlabel','x [m]','ch_ylabel','y [m]', 'ch_title',ch_tit, 'v_caxis', mycolorbar_caxis, 'ch_interpreter', 'latex');
                %               this  parts loads the building structure
                m_multiplot(1,coordNumber,2) = GFigure('m_X',source_X,'m_Y',source_Y,'c_styles',{'k*'});
                
            end
            F = GFigure('m_multiplot',m_multiplot);
        end
        
        % This experiment plots fig 1c and 1d
        function F=experiment_102(obj,niter)
            b_toyExample = 0; %!
            st = dbstack;
            namestr = st.name;
            simulationNumber =str2double(namestr(35:end));
            results_filename = sprintf('./savedResults/dataTest%d.mat', simulationNumber);
            if obj.b_loadData && exist(results_filename, 'file')
                fprintf('%s found. This experiment has been run before but we are not checking the integrity of the saved files. Loading and plotting results...\n', results_filename);
                load(sprintf('./savedResults/dataTest%d.mat', simulationNumber))
            else
                fprintf('%s not found. Running your experiment...\n', results_filename);
                c = CartographySimulations;
                c.f = 800e6;
                c.pilotNoiseSTD =1e-5;
                c.powerNoiseSTD  =0.5;
                c.maxSamplesPerPilot=10;
                c.receiverBandwidth = 20e6;
                c.kernelSigmaLF =14.8;
                c.kernelSigmaLB =0.501;
                c.lambdaLF=1.91e-4;
                c.lambdaLB=3e-3;
                if b_toyExample
                    c.gridSize = [30 30];
                else
                    c.gridSize = [150 150];
                end
                c.numberOfSources=5;
                
                c.computeOnlyFeatures = true;
                
                selectedWalls=1:5;
                
                x_wall_file_name='./modelFiles/x_coord_4walls.mat';
                y_wall_file_name='./modelFiles/y_coord_4walls.mat';
                
                c.n_runs = niter;% 10
                
                c.N_ue =[50,150,250]; %,350,450,550,750,1000,1500,2000,2400];
                if b_toyExample
                    c.N_ue = [20 40];
                end
                
                c.simLocBased=true;  % false disables LocB cartography
                c.testDiffSetOfFeatures=false;
                c.PCAenabler=false;
                c.considerMissingData=false;
                
                if  c.testDiffSetOfFeatures==true
                    c.featureCell=c.featureCombinations(10,4); % generate all combinations of features
                    c.N_feat_used=1:length(c.featureCell);
                    c.N_feat =c.N_feat_used;
                else
                    c.N_feat_used=10;
                    c.N_feat=1;
                end
                if  c.PCAenabler==true
                    c.N_feat =2:3; %  PCA dimension
                elseif c.considerMissingData==true
                    c.N_feat =-85:3:-41; % receiver sensitivy range in dBm
                end
                c.completionMethodSelect=2; % 1 for modified EGM or PGD, and else for Manopt
                c.desiredRank=4;
                c.regularParEval=5.42e-0;
                c.evalcoeffOption=2; % 1 for mean and covariance between completed features
                c.inParallel =true;
                c.forcePrecomputeDistances = false;
                
                [allPointEstimatedLocation,~, evaluationGrid_x,evaluationGrid_y, source_loc, ~, ~,...
                    ~,~,~,~,~,~,~,~,~,~]=c.run(selectedWalls, x_wall_file_name,...
                    y_wall_file_name,simulationNumber);
                
                save (sprintf('./savedResults/dataTest%d', simulationNumber))
            end
            
            %Plotting results
            smoothnessToMapLoc=zeros(c.gridSize(1),c.gridSize(2), 2);
            mycolorbar_caxis = [-10, 50];
            m_X = evaluationGrid_x(:,1);
            m_Y =evaluationGrid_y(1,:);
            source_X=source_loc(1,:);
            source_Y=source_loc(2,:);
            for coordNumber=1:2
                smoothnessToMapLoc(:,:,coordNumber)=(reshape((allPointEstimatedLocation(coordNumber,:)),[c.gridSize(1),c.gridSize(2)]));
                m_Z=smoothnessToMapLoc(:,:,coordNumber)';
                if coordNumber==1
                    ch_tit='$\hat{x}$';
                else
                    ch_tit='$\hat{y}$';
                end
                
                m_multiplot(1,coordNumber,1)= GFigure('m_X',m_X,'m_Y',m_Y,'m_Z',m_Z,'ch_plotType3D','imagesc',...
                    'ch_xlabel','x [m]','ch_ylabel','y [m]', 'ch_title',ch_tit, 'v_caxis', mycolorbar_caxis,  'ch_interpreter', 'latex');
                %               this  parts loads the building structure
                load(x_wall_file_name);
                load(y_wall_file_name);
                
                for ind_wall=1:length(X)/2  % Plots the builing structure;
                    walls_X=X(2*ind_wall-1:2*ind_wall);
                    walls_Y=Y(2*ind_wall-1:2*ind_wall);
                    m_multiplot(1,coordNumber,ind_wall+1) = GFigure('m_X',walls_X','m_Y',walls_Y','c_styles',{'k-'});
                end
                m_multiplot(1,coordNumber,length(X)/2+2) = GFigure('m_X',source_X,'m_Y',source_Y,'c_styles',{'k*'});
                
            end
            F = GFigure('m_multiplot',m_multiplot);
        end
        
        % This experiment plots figs. 4 (for the number of sources L=5), 5 and 6
        % Set niter argument of gism to 50(number of Monte Carlo runs in
        % the experiment)
        function F=experiment_401(obj,niter)
            st = dbstack;
            namestr = st.name;
            simulationNumber =str2double(namestr(35:end));
            results_filename = sprintf('./savedResults/dataTest%d.mat', simulationNumber);
            if obj.b_loadData && exist(results_filename, 'file')
                fprintf('%s found. This experiment has been run before but we are not checking the integrity of the saved files. Loading and plotting results...\n', results_filename);
                load(sprintf('./savedResults/dataTest%d.mat', simulationNumber))
            else
                fprintf('%s not found. Running your experiment...\n', results_filename);
                c = CartographySimulations;
                c.f = 800e6;
                c.pilotNoiseSTD =1e-5;
                c.powerNoiseSTD  =0.5;
                c.maxSamplesPerPilot=10;
                c.receiverBandwidth = 20e6;
                c.kernelSigmaLF =37;
                c.kernelSigmaLB =0.501;
                c.lambdaLF=1.91e-4;
                c.lambdaLB=3e-3;
                c.gridSize = [150 150];
                c.numberOfSources=5;
                
                c.computeOnlyFeatures = false;
                
                selectedWalls=1:5;
                
                x_wall_file_name='./modelFiles/x_coord_4walls.mat';
                y_wall_file_name='./modelFiles/y_coord_4walls.mat';
                
                c.n_runs = niter;% 10
                
                c.N_ue =[100, 300]; %,350,450,550,750,1000,1500,2000,2400];
                
                
                c.simLocBased=true;  % false disables LocB cartography
                c.testDiffSetOfFeatures=false;
                c.PCAenabler=false;
                c.considerMissingData=false;
                
                if  c.testDiffSetOfFeatures==true
                    c.featureCell=c.featureCombinations(10,4); % generate all combinations of features
                    c.N_feat_used=1:length(c.featureCell);
                    c.N_feat =c.N_feat_used;
                else
                    c.N_feat_used=10;
                    c.N_feat=1;
                end
                if  c.PCAenabler==true
                    c.N_feat =2:3; %  PCA dimension
                elseif c.considerMissingData==true
                    c.N_feat =-85:3:-41; % receiver sensitivy range in dBm
                end
                c.completionMethodSelect=2; % 1 for modified EGM or PGD, and else for Manopt
                c.desiredRank=4;
                c.regularParEval=5.42e-0;
                c.evalcoeffOption=2; % 1 for mean and covariance between completed features
                c.inParallel =true;
                c.forcePrecomputeDistances = false;
                
                [~,allLocFreeFeatures, evaluationGrid_x,evaluationGrid_y, source_loc, trueMap, ~,SqErr_AllRuns_locFree,~, ~, locFreeEstimate,...
                    locBasedEstimate,UEFeatures,UELocations,~, ~,~]=c.run(selectedWalls,...
                    x_wall_file_name,y_wall_file_name,simulationNumber);
                
                
                
                save (sprintf('./savedResults/dataTest%d', simulationNumber))
            end
            
            % Plotting results
            mycolorbar_caxis = [-50, -25];
            N_feat=1:length(c.N_feat);
            [~, best_run] = min(SqErr_AllRuns_locFree (end,:,max(N_feat)));
            Nue_index = size(locFreeEstimate,1); % using the maximum number of ue
            Nfeat_index=1;
            
            % Plotting fig. 4
            m_Y_ue_sing_val=svd(UEFeatures{Nue_index,best_run, Nfeat_index});
            F(1)=GFigure('m_Y', m_Y_ue_sing_val','c_styles',{'k-'},  'ch_xlabel','m','ch_ylabel','\sigma_{m}',...
                'ch_interpreter','tex');
            
            % Plotting fig. 5
            m_X = evaluationGrid_x(:,1);
            m_Y =evaluationGrid_y(1,:);
            m_Z_trueMap=trueMap';
            m_X_ue=UELocations{Nue_index,best_run, Nfeat_index}(1,:);
            m_Y_ue=UELocations{Nue_index,best_run, Nfeat_index}(2,:);
            
            ch_tit_map=sprintf('True map at %dMhz', c.f/1e6);
            
            m_multiplot(1,1,1)= GFigure('m_X',m_X,'m_Y',m_Y,'m_Z',m_Z_trueMap,'ch_plotType3D','imagesc',...
                'ch_xlabel','x [m]','ch_ylabel','y [m]', 'ch_title',ch_tit_map, 'v_caxis', mycolorbar_caxis,  'ch_interpreter', 'latex');
            %               this  parts loads the building structure
            load(x_wall_file_name);
            load(y_wall_file_name);
            
            for ind_wall=1:length(X)/2  % Plots the builing structure;
                walls_X=X(2*ind_wall-1:2*ind_wall);
                walls_Y=Y(2*ind_wall-1:2*ind_wall);
                m_multiplot(1,1,ind_wall+1) = GFigure('m_X',walls_X','m_Y',walls_Y', 'c_styles',{'-'}, 'm_colorset', 'w');
            end
            m_multiplot(1,1,length(X)/2+2) = GFigure('m_X',m_X_ue,'m_Y',m_Y_ue,'c_styles',{'k+'});
            
            m_Z_locB=locBasedEstimate{Nue_index, best_run, Nfeat_index}';
            ch_tit_locB='Location-based estimated map';
            
            m_multiplot(1,2,1)= GFigure('m_X',m_X,'m_Y',m_Y,'m_Z',m_Z_locB,'ch_plotType3D','imagesc',...
                'ch_xlabel','x [m]','ch_ylabel','y [m]', 'ch_title',ch_tit_locB, 'v_caxis', mycolorbar_caxis,  'ch_interpreter', 'latex');
            %               this  parts loads the building structure
            
            for ind_wall=1:length(X)/2  % Plots the builing structure;
                walls_X=X(2*ind_wall-1:2*ind_wall);
                walls_Y=Y(2*ind_wall-1:2*ind_wall);
                m_multiplot(1,2,ind_wall+1) = GFigure('m_X',walls_X','m_Y',walls_Y','c_styles',{'k+'}, 'c_styles',{'-'}, 'm_colorset', 'w');
            end
            m_multiplot(1,2,length(X)/2+2) = GFigure('m_X',m_X_ue,'m_Y',m_Y_ue,'c_styles',{'k+'});
            
            m_Z_locF=locFreeEstimate{Nue_index, best_run, Nfeat_index}';
            ch_tit_locF='Location-free estimated map';
            
            m_multiplot(1,3,1)= GFigure('m_X',m_X,'m_Y',m_Y,'m_Z',m_Z_locF,'ch_plotType3D','imagesc',...
                'ch_xlabel','x [m]','ch_ylabel','y [m]', 'ch_title',ch_tit_locF, 'v_caxis', mycolorbar_caxis,  'ch_interpreter', 'latex');
            %               this  parts loads the building structure
            for ind_wall=1:length(X)/2  % Plots the builing structure;
                walls_X=X(2*ind_wall-1:2*ind_wall);
                walls_Y=Y(2*ind_wall-1:2*ind_wall);
                m_multiplot(1,3,ind_wall+1) = GFigure('m_X',walls_X','m_Y',walls_Y','c_styles',{'-'}, 'm_colorset', 'w');
            end
            m_multiplot(1,3,length(X)/2+2) = GFigure('m_X',m_X_ue,'m_Y',m_Y_ue,'c_styles',{'k+'});
            
            F(2) = GFigure('m_multiplot',m_multiplot);
            
            % Plotting fig. 6
            mycolorbar_caxis_feat = [-60, 60];
            source_X=source_loc(1,:);
            source_Y=source_loc(2,:);
            num_feat=size(allLocFreeFeatures,1);
            for ind_feat=1:num_feat
                m_Z_feat=reshape(allLocFreeFeatures(ind_feat,:),[c.gridSize(1),c.gridSize(2)])';
                ch_tit=sprintf('$\\mathbf{\\varphi}_{%1d}$',ind_feat);
                if ind_feat>num_feat/2
                    row_ind=2;
                    ind_feat_plot=ind_feat-num_feat/2;
                else
                    row_ind=1;
                    ind_feat_plot=ind_feat;
                end
                m_multiplot2(row_ind,ind_feat_plot,1)= GFigure('m_X',m_X,'m_Y',m_Y,'m_Z',m_Z_feat,'ch_plotType3D','imagesc',...
                    'ch_xlabel','x [m]','ch_ylabel','y [m]', 'ch_title',ch_tit, 'v_caxis', mycolorbar_caxis_feat,...
                    'ch_interpreter', 'latex','b_colorbar', '1');
                for ind_wall=1:length(X)/2  % Plots the builing structure;
                    walls_X=X(2*ind_wall-1:2*ind_wall);
                    walls_Y=Y(2*ind_wall-1:2*ind_wall);
                    m_multiplot2(row_ind,ind_feat_plot,ind_wall+1) = GFigure('m_X',walls_X','m_Y',walls_Y','c_styles',{'k+'}, 'c_styles',{'k-'});
                end
                m_multiplot2(row_ind,ind_feat_plot,length(X)/2+2) = GFigure('m_X',source_X,'m_Y',source_Y,'c_styles',{'k*'});
            end
            
            F(3) = GFigure('m_multiplot',m_multiplot2);
        end
        
        % This experiment plots the dotdashed curves corresponding to the number
        % of transmitters L=4 in Fig. 7. Set niter argument of gism to
        % 100(number of Monte Carlo runs in the experiment)
        function F=experiment_701(obj,niter)
            st = dbstack;
            namestr = st.name;
            simulationNumber =str2double(namestr(35:end));
            results_filename = sprintf('./savedResults/dataTest%d.mat', simulationNumber);
            if obj.b_loadData && exist(results_filename, 'file')
                fprintf('%s found. This experiment has been run before but we are not checking the integrity of the saved files. Loading and plotting results...\n', results_filename);
                load(sprintf('./savedResults/dataTest%d.mat', simulationNumber))
            else
                fprintf('%s not found. Running your experiment...\n', results_filename);
                c = CartographySimulations;
                c.f = 800e6;
                c.pilotNoiseSTD =1e-5;
                c.powerNoiseSTD  =0.5;
                c.maxSamplesPerPilot=10;
                c.receiverBandwidth = 20e6;
                c.kernelSigmaLF =37;
                c.kernelSigmaLB =0.501;
                c.lambdaLF=1.91e-4;
                c.lambdaLB=3e-3;
                c.gridSize = [50 50];
                c.numberOfSources=4;
                
                c.computeOnlyFeatures = false;
                
                selectedWalls=1:5;
                
                x_wall_file_name='./modelFiles/x_coord_4walls.mat';
                y_wall_file_name='./modelFiles/y_coord_4walls.mat';
                
                c.n_runs = niter;% 10
                
                c.N_ue=[50,150,250,350,450,550,750,1000,1500,2000,2400];
                
                
                c.simLocBased=true;  % false disables LocB cartography
                c.testDiffSetOfFeatures=false;
                c.PCAenabler=false;
                c.considerMissingData=false;
                
                if  c.testDiffSetOfFeatures==true
                    c.featureCell=c.featureCombinations(10,4); % generate all combinations of features
                    c.N_feat_used=1:length(c.featureCell);
                    c.N_feat =c.N_feat_used;
                else
                    c.N_feat_used=size(combnk(1:c.numberOfSources,2),1);
                    c.N_feat=1;
                end
                if  c.PCAenabler==true
                    c.N_feat =2:3; %  PCA dimension
                elseif c.considerMissingData==true
                    c.N_feat =-85:3:-41; % receiver sensitivy range in dBm
                end
                c.completionMethodSelect=2; % 1 for modified EGM or PGD, and else for Manopt
                c.desiredRank=4;
                c.regularParEval=5.42e-0;
                c.evalcoeffOption=2; % 1 for mean and covariance between completed features
                c.inParallel =true;
                c.forcePrecomputeDistances = false;
                
                [~,~, ~,~, ~, trueMap,~,...
                    SqErr_AllRuns_locFree,SqErr_AllRuns_locBased,~,~,~,~,~,...
                    ~, NMSE_locFree,NMSE_locBased]=c.run(selectedWalls,...
                    x_wall_file_name,y_wall_file_name,simulationNumber);
                
                
                
                save (sprintf('./savedResults/dataTest%d', simulationNumber))
            end
            
            % Plotting results
            N_feat=1:length(c.N_feat);
            trueMapAverage = mean(mean(trueMap));
            myStd_locB=std(SqErr_AllRuns_locBased(:,:,N_feat),  0, 2);
            error_locB=3*permute(myStd_locB,[1 3 2])./(sqrt(c.N_ue).*sum(sum((trueMap-trueMapAverage).^2)))';
            myStd_locF=std(SqErr_AllRuns_locFree(:,:,N_feat),  0, 2);
            error_locF= 3*permute(myStd_locF,[1 3 2])./(sqrt(c.N_ue).*sum(sum((trueMap-trueMapAverage).^2)))';
            err_bars= [error_locB'; error_locF'];
            
            c_legend{1} = sprintf('LocB, L= %d',c.numberOfSources);
            c_styles{1} = '>r-.';
            c_legend{2} = sprintf('LocF, L= %d',c.numberOfSources);
            c_styles{2} = 'sb-.';
            
            m_X = c.N_ue;
            m_Y =[permute(NMSE_locBased(:,:,N_feat),[1 3 2])'; permute(NMSE_locFree(:,:,N_feat),[1 3 2])'];
            F = GFigure('m_X',m_X,'m_Y',m_Y,'ch_xlabel','Number of sensor locations, N','ch_ylabel','NMSE',...
                'c_legend',c_legend, 'c_styles', c_styles,'m_errorBars', err_bars);
            
        end
        
        % This experiment plots the solid curves corresponds to the number
        % of transmitters L=7 in Fig. 7. Set niter argument of gism to
        % 100(number of Monte Carlo runs in the experiment)
        function F=experiment_702(obj,niter)
            st = dbstack;
            namestr = st.name;
            simulationNumber =str2double(namestr(35:end));
            results_filename = sprintf('./savedResults/dataTest%d.mat', simulationNumber);
            if obj.b_loadData && exist(results_filename, 'file')
                fprintf('%s found. This experiment has been run before but we are not checking the integrity of the saved files. Loading and plotting results...\n', results_filename);
                load(sprintf('./savedResults/dataTest%d.mat', simulationNumber))
            else
                fprintf('%s not found. Running your experiment...\n', results_filename);
                c = CartographySimulations;
                c.f = 800e6;
                c.pilotNoiseSTD =1e-5;
                c.powerNoiseSTD  =0.5;
                c.maxSamplesPerPilot=10;
                c.receiverBandwidth = 20e6;
                c.kernelSigmaLF =37;
                c.kernelSigmaLB =0.501;
                c.lambdaLF=1.91e-4;
                c.lambdaLB=3e-3;
                c.gridSize = [50 50];
                c.numberOfSources=7;
                
                c.computeOnlyFeatures = false;
                
                selectedWalls=1:5;
                
                x_wall_file_name='./modelFiles/x_coord_4walls.mat';
                y_wall_file_name='./modelFiles/y_coord_4walls.mat';
                
                c.n_runs = niter;% 10
                
                c.N_ue=[50,150,250,350,450,550,750,1000,1500,2000,2400];
                
                
                c.simLocBased=true;  % false disables LocB cartography
                c.testDiffSetOfFeatures=false;
                c.PCAenabler=false;
                c.considerMissingData=false;
                
                if  c.testDiffSetOfFeatures==true
                    c.featureCell=c.featureCombinations(10,4); % generate all combinations of features
                    c.N_feat_used=1:length(c.featureCell);
                    c.N_feat =c.N_feat_used;
                else
                    c.N_feat_used=size(combnk(1:c.numberOfSources,2),1);
                    c.N_feat=1;
                end
                if  c.PCAenabler==true
                    c.N_feat =2:3; %  PCA dimension
                elseif c.considerMissingData==true
                    c.N_feat =-85:3:-41; % receiver sensitivy range in dBm
                end
                c.completionMethodSelect=2; % 1 for modified EGM or PGD, and else for Manopt
                c.desiredRank=4;
                c.regularParEval=5.42e-0;
                c.evalcoeffOption=2; % 1 for mean and covariance between completed features
                c.inParallel =true;
                c.forcePrecomputeDistances = false;
                
                [~,~, ~,~, ~, trueMap,~,...
                    SqErr_AllRuns_locFree,SqErr_AllRuns_locBased,~,~,~,~,~,...
                    ~, NMSE_locFree,NMSE_locBased]=c.run(selectedWalls,...
                    x_wall_file_name,y_wall_file_name,simulationNumber);
                
                
                
                save (sprintf('./savedResults/dataTest%d', simulationNumber))
            end
            
            % Plotting results
            N_feat=1:length(c.N_feat);
            trueMapAverage = mean(mean(trueMap));
            myStd_locB=std(SqErr_AllRuns_locBased(:,:,N_feat),  0, 2);
            error_locB=3*permute(myStd_locB,[1 3 2])./(sqrt(c.N_ue).*sum(sum((trueMap-trueMapAverage).^2)))';
            myStd_locF=std(SqErr_AllRuns_locFree(:,:,N_feat),  0, 2);
            error_locF= 3*permute(myStd_locF,[1 3 2])./(sqrt(c.N_ue).*sum(sum((trueMap-trueMapAverage).^2)))';
            err_bars= [error_locB'; error_locF'];
            
            c_legend{1} = sprintf('LocB, L= %d',c.numberOfSources);
            c_styles{1} = '>r-';
            c_legend{2} = sprintf('LocF, L= %d',c.numberOfSources);
            c_styles{2} = 'sb-';
            
            m_X = c.N_ue;
            m_Y =[permute(NMSE_locBased(:,:,N_feat),[1 3 2])'; permute(NMSE_locFree(:,:,N_feat),[1 3 2])'];
            F = GFigure('m_X',m_X,'m_Y',m_Y, 'c_legend',c_legend,...
                'ch_xlabel','Number of sensor locations, N','ch_ylabel','NMSE',...
                'c_styles', c_styles,'m_errorBars', err_bars);
            
        end
        
        % This experiment plots Fig. 8. Set niter argument of gism to
        % 200(number of Monte Carlo runs in the experiment)
        function F=experiment_801(obj,niter)
            st = dbstack;
            namestr = st.name;
            simulationNumber =str2double(namestr(35:end));
            results_filename = sprintf('./savedResults/dataTest%d.mat', simulationNumber);
            if obj.b_loadData && exist(results_filename, 'file')
                fprintf('%s found. This experiment has been run before but we are not checking the integrity of the saved files. Loading and plotting results...\n', results_filename);
                load(sprintf('./savedResults/dataTest%d.mat', simulationNumber))
            else
                fprintf('%s not found. Running your experiment...\n', results_filename);
                c = CartographySimulations;
                c.f = 800e6;
                c.pilotNoiseSTD =1e-5;
                c.powerNoiseSTD  =0.5;
                c.gridSize = [50 50];
                c.numberOfSources=5;
                
                c.computeOnlyFeatures = false;
                
                x_wall_file_name='./modelFiles/x_coord_5walls.mat';
                y_wall_file_name='./modelFiles/y_coord_5walls.mat';
                
                c.n_runs = niter;
                
                c.N_ue=300;
                
                
                c.simLocBased=true;  % false disables LocB cartography
                c.testDiffSetOfFeatures=false;
                c.PCAenabler=false;
                c.considerMissingData=false;
                
                if  c.testDiffSetOfFeatures==true
                    c.featureCell=c.featureCombinations(10,4); % generate all combinations of features
                    c.N_feat_used=1:length(c.featureCell);
                    c.N_feat =c.N_feat_used;
                else
                    c.N_feat_used=size(combnk(1:c.numberOfSources,2),1);
                    c.N_feat=1;
                end
                if  c.PCAenabler==true
                    c.N_feat =2:3; %  PCA dimension
                elseif c.considerMissingData==true
                    c.N_feat =-85:3:-41; % receiver sensitivy range in dBm
                end
                c.completionMethodSelect=2; % 1 for modified EGM or PGD, and else for Manopt
                c.desiredRank=4;
                c.regularParEval=5.42e-0;
                c.evalcoeffOption=2; % 1 for mean and covariance between completed features
                c.inParallel =true;
                c.forcePrecomputeDistances = false;
                
                n_walls=5;
                wallCell={combnk(1:n_walls,0),combnk(1:n_walls,1),combnk(1:n_walls,2),combnk(1:n_walls,3),...
                    combnk(1:n_walls,4),combnk(1:n_walls,5)};
                n_subetWalls=length(wallCell);
                
                myBW=[50e6 100e6 200e6 700e6]; % bandwidth in Hz
                kernelSigmaLF=[27,41,53,28];
                lambdaLF=[3.81e-4,6.1e-5,1.1e-5, 5e-4];
                kernelSigmaLB=[10.1,8.9,9,7];
                lambdaLB=[1.8e-3,9.1e-4,7.1e-4,2.1e-4];
                
                NMSE_diffSetWallsBW=zeros(2,n_subetWalls,length(myBW));
                for ind_bw=1:length(myBW)
                    c.receiverBandwidth = myBW(ind_bw);
                    c.maxSamplesPerPilot=c.receiverBandwidth/(2*1e6);
                    c.kernelSigmaLF =kernelSigmaLF(ind_bw);
                    c.kernelSigmaLB =kernelSigmaLB(ind_bw);
                    c.lambdaLF= lambdaLF(ind_bw);
                    c.lambdaLB=lambdaLB(ind_bw);
                    
                    NMSE_diffSetWalls=zeros(2,n_subetWalls);
                    
                    for ind_subset=1:n_subetWalls
                        currentSubset=wallCell{ind_subset};
                        if isempty(currentSubset)
                            currentSubset=0;
                        end
                        NMSE_oneSubSet=zeros(2,size(currentSubset,1));
                        for ind_currentSub=1:size(currentSubset,1)
                            selectedWalls=currentSubset(ind_currentSub,:);
                            [~,~, ~,~, ~, ~,~,~,~,~,~,~,~,~,...
                                ~, NMSE_locFree,NMSE_locBased]=c.run(selectedWalls,...
                                x_wall_file_name,y_wall_file_name,simulationNumber);
                            
                            NMSE_oneSubSet(:,ind_currentSub)=[NMSE_locFree;NMSE_locBased];
                        end
                        NMSE_diffSetWalls(:,ind_subset)=mean(NMSE_oneSubSet,2);
                    end
                    
                    NMSE_diffSetWallsBW(:,:,ind_bw)=NMSE_diffSetWalls;
                end
                
                save (sprintf('./savedResults/dataTest%d', simulationNumber))
            end
            
            % Plotting results
            c_legend{1} = 'LocF';
            c_styles{1} = 'bs-';
            c_legend{2} = 'LocB';
            c_styles{2} = 'rv-';
            
            plotNum=length(myBW);
            for ind_bw=1:plotNum
                m_X =0:n_subetWalls-1;
                m_Y=NMSE_diffSetWallsBW(:,:,ind_bw);
                ch_tit=sprintf('$B=%d $ $MHz$', myBW(ind_bw)/1e6);
                m_multiplot(1,ind_bw)= GFigure('m_X',m_X,'m_Y',m_Y,...
                    'ch_xlabel','Number of Walls','ch_ylabel','NMSE', 'ch_title',ch_tit,...
                    'v_ylim',[0,0.9],'c_legend',c_legend, 'c_styles', c_styles, 'ch_interpreter', 'latex');
            end
            F = GFigure('m_multiplot',m_multiplot);
        end
        
        % This experiment plots Fig. 9. Set niter argument of gism to
        % 100(number of Monte Carlo runs in the experiment)
        function F=experiment_901(obj,niter)
            st = dbstack;
            namestr = st.name;
            simulationNumber =str2double(namestr(35:end));
            results_filename = sprintf('./savedResults/dataTest%d.mat', simulationNumber);
            if obj.b_loadData && exist(results_filename, 'file')
                fprintf('%s found. This experiment has been run before but we are not checking the integrity of the saved files. Loading and plotting results...\n', results_filename);
                load(sprintf('./savedResults/dataTest%d.mat', simulationNumber))
            else
                fprintf('%s not found. Running your experiment...\n', results_filename);
                c = CartographySimulations;
                c.f = 800e6;
                c.pilotNoiseSTD =1e-5;
                c.powerNoiseSTD  =0.5;
                c.maxSamplesPerPilot=10;
                c.receiverBandwidth = 20e6;
                c.kernelSigmaLF =37;
                c.kernelSigmaLB =0.501;
                c.lambdaLF=1.91e-4;
                c.lambdaLB=3e-3;
                c.gridSize = [50 50];
                c.numberOfSources=5;
                
                c.computeOnlyFeatures = false;
                
                selectedWalls=1:5;
                
                x_wall_file_name='./modelFiles/x_coord_4walls.mat';
                y_wall_file_name='./modelFiles/y_coord_4walls.mat';
                
                c.n_runs = niter;% 10
                
                c.N_ue=[200,300];
                
                
                c.simLocBased=false;  % false disables LocB cartography
                c.testDiffSetOfFeatures=true;
                c.PCAenabler=false;
                c.considerMissingData=false;
                
                if  c.testDiffSetOfFeatures==true
                    c.featureCell=c.featureCombinations(10,4); % generate all combinations of features
                    c.N_feat_used=1:length(c.featureCell);
                    c.N_feat =c.N_feat_used;
                else
                    c.N_feat_used=size(combnk(1:c.numberOfSources,2),1);
                    c.N_feat=1;
                end
                if  c.PCAenabler==true
                    c.N_feat =2:3; %  PCA dimension
                elseif c.considerMissingData==true
                    c.N_feat =-85:3:-41; % receiver sensitivy range in dBm
                end
                c.completionMethodSelect=2; % 1 for modified EGM or PGD, and else for Manopt
                c.desiredRank=4;
                c.regularParEval=5.42e-0;
                c.evalcoeffOption=2; % 1 for mean and covariance between completed features
                c.inParallel =true;
                c.forcePrecomputeDistances = false;
                
                [~,~, ~,~, ~, ~,~,~,~,~,~,~,~,~,...
                    ~, NMSE_locFree,~]=c.run(selectedWalls,...
                    x_wall_file_name,y_wall_file_name,simulationNumber);
                
                save (sprintf('./savedResults/dataTest%d', simulationNumber))
            end
            
            % Plotting results
            NMSE_locFreeResh=permute(NMSE_locFree,[1 3 2]);
            sartingFeatNum=4; maxFeatNum=10;
            featRange=sartingFeatNum:maxFeatNum;
            NMSE_diffM=zeros(length(c.N_ue),length(featRange));
            startin_point=1;
            
            c_legend{1} = sprintf('N=%d', c.N_ue(1));
            c_styles{1} = 'bv-';
            c_legend{2} = sprintf('N=%d', c.N_ue(2));
            c_styles{2} = 'ro-';
            for indfeat=1:length(featRange)
                averagingRange=size(combnk(1:maxFeatNum,featRange(indfeat)),1);
                end_point=startin_point+averagingRange-1;
                NMSE_diffM(:,indfeat)=mean(NMSE_locFreeResh(:,startin_point:end_point),2);
                startin_point=end_point+1;
            end
            m_X=featRange;
            m_Y=NMSE_diffM;
            F = GFigure('m_X',m_X,'m_Y',m_Y,...
                'ch_xlabel','M','ch_ylabel','NMSE','c_legend',c_legend,...
                'c_styles', c_styles, 'ch_interpreter', 'latex');
            
        end
        
        % This experiment plots Fig. 10. 
        function F=experiment_1001(obj,niter)
            st = dbstack;
            namestr = st.name;
            simulationNumber =str2double(namestr(35:end));
            results_filename = sprintf('./savedResults/dataTest%d.mat', simulationNumber);
            if obj.b_loadData && exist(results_filename, 'file')
                fprintf('%s found. This experiment has been run before but we are not checking the integrity of the saved files. Loading and plotting results...\n', results_filename);
                load(sprintf('./savedResults/dataTest%d.mat', simulationNumber))
            else
                fprintf('%s not found. Running your experiment...\n', results_filename);
                c = CartographySimulations;
                c.f = 800e6;
                c.pilotNoiseSTD =1e-5;
                c.powerNoiseSTD  =0.5;
                c.maxSamplesPerPilot=10;
                c.receiverBandwidth = 20e6;
                c.kernelSigmaLF =37;
                c.kernelSigmaLB =0.501;
                c.lambdaLF=1.91e-4;
                c.lambdaLB=3e-3;
                c.gridSize = [200 200];
                c.numberOfSources=5;
                
                c.computeOnlyFeatures = true;
                
                selectedWalls=1:5;
                
                x_wall_file_name='./modelFiles/x_coord_more_walls.mat';
                y_wall_file_name='./modelFiles/y_coord_more_walls.mat';
                
                c.n_runs = niter;% 10
                
                c.N_ue= 200;
                
                
                c.simLocBased=false;  % false disables LocB cartography
                c.testDiffSetOfFeatures=false;
                c.PCAenabler=false;
                c.considerMissingData=false;
                
                if  c.testDiffSetOfFeatures==true
                    c.featureCell=c.featureCombinations(10,4); % generate all combinations of features
                    c.N_feat_used=1:length(c.featureCell);
                    c.N_feat =c.N_feat_used;
                else
                    c.N_feat_used=size(combnk(1:c.numberOfSources,2),1);
                    c.N_feat=1;
                end
                if  c.PCAenabler==true
                    c.N_feat =2:3; %  PCA dimension
                elseif c.considerMissingData==true
                    c.N_feat =-85:3:-41; % receiver sensitivy range in dBm
                end
                c.completionMethodSelect=2; % 1 for modified EGM or PGD, and else for Manopt
                c.desiredRank=4;
                c.regularParEval=5.42e-0;
                c.evalcoeffOption=2; % 1 for mean and covariance between completed features
                c.inParallel =true;
                c.forcePrecomputeDistances = false;
                
                [~,allLocFreeFeatures, evaluationGrid_x,evaluationGrid_y,source_loc, ~,~,~,~,~,~,~,~,~,...
                    ~, ~,~]=c.run(selectedWalls,...
                    x_wall_file_name,y_wall_file_name,simulationNumber);
                
                save (sprintf('./savedResults/dataTest%d', simulationNumber))
            end
            
            % Plotting results
            load(x_wall_file_name);
            load(y_wall_file_name);
            m_X =evaluationGrid_x(:,1);
            m_Y =evaluationGrid_y(1,:);
            mycolorbar_caxis_feat = [-60, 60];
            mycolorbar_caxis2 = [-83, 94];
            source_X=source_loc(1,:);
            source_Y=source_loc(2,:);
            num_feat=size(allLocFreeFeatures,1);
            for ind_feat=1:num_feat
                m_Z_feat=reshape(allLocFreeFeatures(ind_feat,:),[c.gridSize(1),c.gridSize(2)])';
                ch_tit=sprintf('$\\mathbf{\\varphi}_{%1d}$',ind_feat);
                if ind_feat>num_feat/2
                    row_ind=2;
                    ind_feat_plot=ind_feat-num_feat/2;
                else
                    row_ind=1;
                    ind_feat_plot=ind_feat;
                end
                m_multiplot(row_ind,ind_feat_plot,1)= GFigure('m_X',m_X,'m_Y',m_Y,'m_Z',m_Z_feat,'ch_plotType3D','imagesc',...
                    'ch_xlabel','x [m]','ch_ylabel','y [m]', 'ch_title',ch_tit, 'v_caxis', mycolorbar_caxis_feat,...
                    'ch_interpreter', 'latex','b_colorbar', '1');
                for ind_wall=1:length(X)/2  % Plots the builing structure;
                    walls_X=X(2*ind_wall-1:2*ind_wall);
                    walls_Y=Y(2*ind_wall-1:2*ind_wall);
                    m_multiplot(row_ind,ind_feat_plot,ind_wall+1) = GFigure('m_X',walls_X','m_Y',walls_Y','c_styles',{'k+'}, 'c_styles',{'k-'});
                end
                m_multiplot(row_ind,ind_feat_plot,length(X)/2+2) = GFigure('m_X',source_X,'m_Y',source_Y,'c_styles',{'k*'});
            end
            
            F(1) = GFigure('m_multiplot',m_multiplot);
            
            meanofFeat=mean(allLocFreeFeatures,2);
            centeredFeat=allLocFreeFeatures-meanofFeat;
            [U,~,~] = svd(centeredFeat);
            U1=U(:,1:c.numberOfSources-1);
            newNFeat=size(U1,2);
            newFeat=U1'*centeredFeat;
            for ind_newFeat=1:newNFeat
                ch_tit=sprintf('$\\bar\\mathbf{\\varphi}_{%1d}$',ind_newFeat);
                m_Z_feat=reshape(newFeat(ind_newFeat,:),[c.gridSize(1),c.gridSize(2)])';
                m_multiplot2(1,ind_newFeat,1)= GFigure('m_X',m_X,'m_Y',m_Y,'m_Z',m_Z_feat,'ch_plotType3D','imagesc',...
                    'ch_xlabel','x [m]','ch_ylabel','y [m]', 'ch_title',ch_tit, 'v_caxis', mycolorbar_caxis2,...
                    'ch_interpreter', 'latex','b_colorbar', '1');
                for ind_wall=1:length(X)/2  % Plots the builing structure;
                    walls_X=X(2*ind_wall-1:2*ind_wall);
                    walls_Y=Y(2*ind_wall-1:2*ind_wall);
                    m_multiplot2(1,ind_newFeat,ind_wall+1) = GFigure('m_X',walls_X','m_Y',walls_Y','c_styles',{'k+'}, 'c_styles',{'k-'});
                end
                m_multiplot2(1,ind_newFeat,length(X)/2+2) = GFigure('m_X',source_X,'m_Y',source_Y,'c_styles',{'k*'});
            end
            F(2) = GFigure('m_multiplot',m_multiplot2);
        end
        
        % This experiment plots Fig. 11. Set niter argument of gism to
        % 200(number of Monte Carlo runs in the experiment)
        function F=experiment_1101(obj,niter)
            st = dbstack;
            namestr = st.name;
            simulationNumber =str2double(namestr(35:end));
            results_filename = sprintf('./savedResults/dataTest%d.mat', simulationNumber);
            if obj.b_loadData && exist(results_filename, 'file')
                fprintf('%s found. This experiment has been run before but we are not checking the integrity of the saved files. Loading and plotting results...\n', results_filename);
                load(sprintf('./savedResults/dataTest%d.mat', simulationNumber))
            else
                fprintf('%s not found. Running your experiment...\n', results_filename);
                c = CartographySimulations;
                c.f = 800e6;
                c.pilotNoiseSTD =1e-5;
                c.powerNoiseSTD  =0.5;
                c.maxSamplesPerPilot=10;
                c.receiverBandwidth = 20e6;
                c.kernelSigmaLF =26;
                c.kernelSigmaLB =0.501;
                c.lambdaLF=1.6e-3;
                c.lambdaLB=3e-3;
                c.gridSize = [50 50];
                c.numberOfSources=5;
                
                c.computeOnlyFeatures = false;
                
                selectedWalls=1:5;
                
                x_wall_file_name='./modelFiles/x_coord_4walls.mat';
                y_wall_file_name='./modelFiles/y_coord_4walls.mat';
                
                c.n_runs = niter;
                
                c.N_ue=[50,150,250,350,450,550,750,1000,1500,2000,2400];
                
                
                c.simLocBased=false;
                c.testDiffSetOfFeatures=false;
                c.PCAenabler=false;
                c.considerMissingData=false;
                
                if  c.testDiffSetOfFeatures==true
                    c.featureCell=c.featureCombinations(10,4); % generate all combinations of features
                    c.N_feat_used=1:length(c.featureCell);
                    c.N_feat =c.N_feat_used;
                else
                    c.N_feat_used=size(combnk(1:c.numberOfSources,2),1);
                    c.N_feat=1;
                end
                if  c.PCAenabler==true
                    c.N_feat =2:3; %  PCA dimension
                elseif c.considerMissingData==true
                    c.N_feat =-85:3:-41; % receiver sensitivy range in dBm
                end
                c.completionMethodSelect=2; % 1 for modified EGM or PGD, and else for Manopt
                c.desiredRank=4;
                c.regularParEval=5.42e-0;
                c.evalcoeffOption=2; % 1 for mean and covariance between completed features
                c.inParallel =true;
                c.forcePrecomputeDistances = false;
                
                [~,~, ~,~, ~, ~,~,~,~,~,~,~,~,~,...
                    ~, NMSE_locFree_no_pca,~]=c.run(selectedWalls,...
                    x_wall_file_name,y_wall_file_name,simulationNumber);
                
                %overwrite the necessary properties to run the
                %experiments with PCA with dimension 2 and 3
                c.kernelSigmaLF =23;
                c.lambdaLF=1.64e-2;
                c.PCAenabler=true;
                c.N_feat =2:3;
                
                [~,~, ~,~, ~, ~,~,~,~,~,~,~,~,~,...
                    ~, NMSE_locFree_pca2_3,~]=c.run(selectedWalls,...
                    x_wall_file_name,y_wall_file_name,simulationNumber);
                
                %overwrite the necessary properties to run the
                %experiments with PCA with dimension 4
                c.kernelSigmaLF =25;
                c.lambdaLF=9.1e-4;
                c.PCAenabler=true;
                c.N_feat =4;
                
                [~,~, ~,~, ~, ~,~,~,~,~,~,~,~,~,...
                    ~, NMSE_locFree_pca4,~]=c.run(selectedWalls,...
                    x_wall_file_name,y_wall_file_name,simulationNumber);
                
                
                save (sprintf('./savedResults/dataTest%d', simulationNumber))
            end
            
            % Plotting results
            m_X = c.N_ue;
            m_Y =[NMSE_locFree_no_pca';permute(NMSE_locFree_pca2_3,[3 1 2]);...
                NMSE_locFree_pca4'];
            c_legend{1} = 'M=10';
            c_legend{2} = 'r=2';
            c_legend{3} = 'r=3';
            c_legend{4} = 'r=4';
            c_styles{1} = 'bd-';
            c_styles{2} = 'kx-';
            c_styles{3} = 'ms-';
            c_styles{4} = 'r>-';
            
            F = GFigure('m_X',m_X,'m_Y',m_Y,...
                'ch_xlabel','Number of sensor locations, N','ch_ylabel','NMSE','c_legend',c_legend,...
                'c_styles', c_styles, 'ch_interpreter', 'latex');
            
        end
        
        % This experiment plots Fig. 12. Set niter argument of gism to
        % 200(number of Monte Carlo runs in the experiment)
        function F=experiment_1201(obj,niter)
            st = dbstack;
            namestr = st.name;
            simulationNumber =str2double(namestr(35:end));
            results_filename = sprintf('./savedResults/dataTest%d.mat', simulationNumber);
            if obj.b_loadData && exist(results_filename, 'file')
                fprintf('%s found. This experiment has been run before but we are not checking the integrity of the saved files. Loading and plotting results...\n', results_filename);
                load(sprintf('./savedResults/dataTest%d.mat', simulationNumber))
            else
                fprintf('%s not found. Running your experiment...\n', results_filename);
                c = CartographySimulations;
                c.f = 800e6;
                c.pilotNoiseSTD =1e-5;
                c.powerNoiseSTD  =0.5;
                c.maxSamplesPerPilot=10;
                c.receiverBandwidth = 20e6;
                c.kernelSigmaLF =37;
                c.kernelSigmaLB =0.501;
                c.lambdaLF=1.91e-4;
                c.lambdaLB=3e-3;
                c.gridSize = [50 50];
                c.numberOfSources=5;
                
                c.computeOnlyFeatures = false;
                
                selectedWalls=1:5;
                
                x_wall_file_name='./modelFiles/x_coord_4walls.mat';
                y_wall_file_name='./modelFiles/y_coord_4walls.mat';
                
                c.n_runs = niter;
                
                c.N_ue=[100, 400, 700, 1500, 2400];
                
                
                c.simLocBased=false;
                c.testDiffSetOfFeatures=false;
                c.PCAenabler=false;
                c.considerMissingData=true;
                
                if  c.testDiffSetOfFeatures==true
                    c.featureCell=c.featureCombinations(10,4); % generate all combinations of features
                    c.N_feat_used=1:length(c.featureCell);
                    c.N_feat =c.N_feat_used;
                else
                    c.N_feat_used=size(combnk(1:c.numberOfSources,2),1);
                    c.N_feat=1;
                end
                if  c.PCAenabler==true
                    c.N_feat =2:3; %  PCA dimension
                elseif c.considerMissingData==true
                    c.N_feat =-85:3:-41; % receiver sensitivy range in dBm
                end
                c.completionMethodSelect=1; % 1 for modified EGM or PGD, and else for Manopt
                c.desiredRank=4;
                c.regularParEval=5.42e-0;
                c.evalcoeffOption=2; % 1 for mean and covariance between completed features
                c.inParallel =true;
                c.forcePrecomputeDistances = false;
                
                [~,~, ~,~, ~, ~,~,~,~,~,~,~,~,~,...
                    averageMisses, NMSE_locFree_svp,~]=c.run(selectedWalls,...
                    x_wall_file_name,y_wall_file_name,simulationNumber);
                
                %overwrite the necessary properties to run the
                %experiments with Manopt as matrix completion technique
                
                c.completionMethodSelect=2;
                [~,~, ~,~, ~, ~,~,~,~,~,~,~,~,~,...
                    ~, NMSE_locFree_lrgeo,~]=c.run(selectedWalls,...
                    x_wall_file_name,y_wall_file_name,simulationNumber);
                
                save (sprintf('./savedResults/dataTest%d', simulationNumber))
            end
            
            % Plotting results
            matAverageMisses=cell2mat(averageMisses);
            perAverageMisses=permute(matAverageMisses,[1 3 2]);
            averageMisses2=mean(perAverageMisses,3);
            m_X = c.N_feat;
            m_Y =averageMisses2(randi(length(c.N_ue)),:);
            m_multiplot(1,1) = GFigure('m_X',m_X,'m_Y',m_Y,...
                'ch_xlabel','Sensor sensitivity, \Gamma [dBm]','ch_ylabel','Average number of missing features',...
                'c_styles', {'k--'}, 'ch_interpreter', 'latex');
            
            m_Y_nmse =[permute(NMSE_locFree_svp,[1 3 2]); permute(NMSE_locFree_lrgeo,[1 3 2])];
            
            for ind_ue=1:length(c.N_ue)
                c_legend_svp{ind_ue} = sprintf('SVP, N= %d',c.N_ue(ind_ue));
                c_legend_lrgeo{ind_ue} = sprintf('LRGeomCG, N= %d',c.N_ue(ind_ue));
            end
            c_legend=[c_legend_svp,c_legend_lrgeo];
            
            c_styles{1} = 'bs-'; c_styles{2} = 'r<-';  c_styles{3} = 'm*-';  c_styles{4} = 'k+-';
            c_styles{5} = 'gv-'; c_styles{6} = 'bs--'; c_styles{7} = 'r<--'; c_styles{8} = 'm*--';
            c_styles{9} = 'k+--';c_styles{10} = 'gv--';
            
            m_multiplot(2,1)=GFigure('m_X',m_X,'m_Y',m_Y_nmse,...
                'ch_xlabel','Sensor sensitivity, \Gamma [dBm]','ch_ylabel','NMSE','c_legend',c_legend,...
                'c_styles', c_styles, 'ch_interpreter', 'latex');
            F = GFigure('m_multiplot',m_multiplot);
        end
    end
    
    
end
