classdef CartographySimulations
    %CartographySimulations Each static method is a complete simulation.
    %   This class includes several static methods, each of them is a
    %   simulation. Apart from those, the method CartographySimulations.run is
    %   executed once the simulation parameters (properties) are defined.
    
    properties % default values = baseline testcase
        rng_seed = 6; %seed of the random number generator
        forcePrecomputeDistances = true; % if set to true, it will ignore the precomputed distances file
        inParallel = true; % compute in parallel (use PARFOR)
        
        % properties pertaining the statistical analysis
        
        gridSize 	% Size of the grid
        N_ue        % Vector containing different values for the number of UEs/sensors
        n_runs           % Number of Monte Carlo runs to be evaluated in each experiment
        
        % Number of features. Set N_feat = 3 (default) for comparison of
        % locFree vs locBased;
        % set N_feat = 1:3 for evaluating impact of M
        N_feat_used % total number of features per sensor location
        N_feat % used for PCA dimension or Rx sensistivity in missing feature
        numberOfSources
        
        % Should simulate the location based estimation.
        % set simLocBased = true (default) for comparison of
        % locFree vs locBased;
        % set simLocBased = false for evaluating impact of features M, PCA, and missing features for LocF scheme.
        simLocBased
        % Hyperparameters of the kernel machines
        computeOnlyFeatures
        
        kernelSigmaLF  % Sigma width for the gaussian kernel in LocFree
        kernelSigmaLB% Sigma width for the gaussian kernel in LocBased
        lambdaLF      % Ridge regularization parameter (locFree )
        lambdaLB      % Ridge regularization parameter (locBased)
        
        % Properties of the test scenario
        
        f  %Carrier frequency
        locationNoiseSTD  % Standard deviation of location feature noise
        powerNoiseSTD     % Standard deviation of measured power noise
        pilotNoiseSTD
        maxSamplesPerPilot
        % Bandwidth of the receiver (will determine the resolution
        % of the digital impulse response)
        receiverBandwidth
        % each element of the cell array is a vector indicating which features
        % (i.e. signals from which sources) are considered for the
        featureCell % = {[2], [1 2], [1 2 3],[1 2 3 4],[1 2 4 5 6],[1:6],[1:7],[1 2 3 4 8 9 10],[4:10],[1:10]}; %#ok<NBRAK>
        
        
        testDiffSetOfFeatures;
        PCAenabler % if true, the PCA is taken into account
        considerMissingData % accommodating the misses:  see classfication with low rank and missing data
        completionMethodSelect %selects either ManOpt or other specified method (Poject Gradient Descent, EGM, etc. )
        desiredRank % desired rank of the matrix to be completed
        regularParEval % %regularization parameter for completing missing features when evaluating the map
        %         featMissing=4; % number of missing features, which can be changed
        evalcoeffOption % choose the option of determining missing features in the map evaluation
        % 1 for computing mean  and covariance between completed features;
        % else  or computing mean  and covariance between coefficients correponding
        %to completed features;
        parGamma %parameter used in the kernel of missing features: low rank and missing data
        
    end
    
   
    
    methods (Static)
        
        function featureCell = featureCombinations(featNum,startingFeat)
            featureCell={};
            for indFeat=startingFeat:featNum
                selectedComb=combnk(1:featNum,indFeat);
                for indc=1:size(selectedComb,1)
                    rowSelec=selectedComb(indc,:);
                    featureCell{end+1}=rowSelec;
                end
            end
        end
        
    end %methods (Static)
    
    methods
        function generator = baselineGenerator(obj, selectedWalls, x_wall_file_name, y_wall_file_name)
           
            generator = MultiWallGenerator;
           
            generator.boundary=[
                -10,48
                -10,30
                -0,3
                ];
            x1 = generator.boundary(1,1);
            x2 = generator.boundary(1,2);
            y1 = generator.boundary(2,1);
            y2 = generator.boundary(2,2);
            
            generator.f=obj.f;           % Frequency in Hz
            
            if obj.numberOfSources==4
                generator.xt=  [0, 28, 32, -3];
                generator.yt=  [0, 10, -7, 30];    
            elseif obj.numberOfSources==7
                generator.xt=  [0, 28, 32, -3, 42, 10, 18];
                generator.yt=  [0, 10, -7, 30, 30, 20, 2];   
            else
                generator.xt=  [0, 28, 32 ,-3, 42]; %
                generator.yt=  [0, 10, -7, 30, 30]; %
            end
            generator.ptx=zeros(1, length(generator.xt));   % Transmitter powers in dBW
            
            generator.zt=1.5; %transmitter Height
            generator.x1=x1; % First boundary  along x
            generator.x2=x2; % Second boundary  along x
            generator.path_loss_exp = 2;
            generator.n_gridpoints_x=obj.gridSize(1);
            generator.n_gridpoints_y=obj.gridSize(2);
            generator.zr=1.5; %receiverHeight;
            generator.sampling_period=1/obj.receiverBandwidth; %
            generator.maxSamplesPerPilot=obj.maxSamplesPerPilot;
            generator.street_limits = [y1 y2]; %
            if selectedWalls==0
                generator.estimateLocFreeSpace=1;
            else
                generator.estimateLocFreeSpace=0;
            end
            X_wallcoord=load(x_wall_file_name);
            Y_wallcoord=load(y_wall_file_name);
            arrangedCoord=[2*selectedWalls-1;2*selectedWalls];
            if length(selectedWalls)<=2
                if selectedWalls==0
                    coordToUse=(11:18)';
                else
                    coordToUse=[arrangedCoord(:);(11:18)']; % (11 to 18)-th coordinates accounts for small structure ...
                end                                        % allowing the multiwall method to work with few walls
            else
                coordToUse=arrangedCoord(:);
            end
            generator.X=X_wallcoord.X(coordToUse);
            generator.Y=Y_wallcoord.Y(coordToUse);
        end
        %changed streetlimits from [0 300] to [0 75].
   
        function [allPointEstimatedLocation,allLocFreeFeatures, evaluationGrid_x,evaluationGrid_y, source_loc, trueMap, estimatedPilotSignals,...
                SqErr_AllRuns_locFree,SqErr_AllRuns_locBased, meanErrOnEvalFeat, locFreeEstimate,locBasedEstimate,UEFeatures,UELocations,...
                averageMisses, NMSE_locFree,NMSE_locBased]=run(obj, selectedWalls, x_wall_file_name,y_wall_file_name, simulationNumber) %"do-not-touch"
            % RUN runs simulation using the parameters stored in obj.
            % run (obj, simNum) runs simulation and saves generated data
            % in the savedResults folder, in a file ending by simNum
            
            mySim = Simulator;
            mySim.generator = obj.baselineGenerator(selectedWalls, x_wall_file_name,y_wall_file_name);
            mySim.features=1:obj.N_feat_used;
            
            mySim.sampler = SpectrumMapSampler;
            mySim.sampler.powerNoiseSTD =obj.powerNoiseSTD;
            mySim.sampler.pilotNoiseSTD =obj.pilotNoiseSTD;
            mySim.sampler.maxSamplesPerPilot = obj.maxSamplesPerPilot;
            mySim.PCAenable=obj.PCAenabler;
            mySim.featureExtractor=FeatureExtractor;
            mySim.featureExtractor.sampling_period=mySim.generator.sampling_period;
            gaussianKernelLF = @(x, y) ...
                exp(-norm(x-y).^2/(obj.kernelSigmaLF^2));
            gaussianKernelLB = @(x, y) ...
                exp(-norm(x-y).^2/(obj.kernelSigmaLB^2));
            
            mySim.locFreeEstimator = LocationFreeEstimator;
            mySim.locFreeEstimator.kernel = gaussianKernelLF;
            mySim.locFreeEstimator.regularizationParameter = obj.lambdaLF;
            mySim.locFreeEstimator.enableMissingData = obj.considerMissingData;
            mySim.locFreeEstimator.completionMethodSelect=obj.completionMethodSelect;
            mySim.locFreeEstimator.desiredRank=obj.desiredRank;
            mySim.locFreeEstimator.regParEval=obj.regularParEval;
            mySim.locFreeEstimator.evalOption=obj.evalcoeffOption;
            
            mySim.locBasedEstimator = LocationBasedEstimator;
            mySim.locEstimator = LocationEstimator;
            mySim.locBasedEstimator.kernel = gaussianKernelLB;
            mySim.locBasedEstimator.regularizationParameter = obj.lambdaLB;
            mySim.locBasedEstimator.locationNoiseSTD = obj.locationNoiseSTD;
            
            rng(obj.rng_seed)
            points_filename = sprintf('allPointEstLoc%d.mat', simulationNumber);
            points_filename2 = sprintf('allLocFreeFeatures%d.mat', simulationNumber);
            
            if obj.simLocBased==true
                fprintf('Precomputing %s and %s \n', points_filename, points_filename2);
                [allPointEstimatedLocation,allLocFreeFeatures,~] = mySim.precomputeEstimatedLocations(); %#ok<NASGU>
                %                     save(points_filename, 'allPointEstimatedLocation')
                %                     save(points_filename2,'allLocFreeFeatures')
                fprintf('Precomputation of %s and %s\n is done \n', points_filename,points_filename2)
                %                 end
            else
                disp('This simulation does not involve LocBased estimation. Not precomputing location data')
                mySim.simLocBased = false;
                fprintf('Precomputating LocF features only...%s\n', points_filename2);
                
                %                     else
                %                         error('Fatal (logics)')
                %                     end
                [~,allLocFreeFeatures,~] = mySim.precomputeEstimatedLocations(); %#ok<NASGU>
                allPointEstimatedLocation=0;
                
                %
                %                     save(points_filename2,'allLocFreeFeatures')
                
                fprintf('Precomputation of %s is done \n ', points_filename2)
            end
            
            
            fprintf('Generating the true map based on the specified grid size \n')
            [power_out,trueMap, evaluationGrid_x, evaluationGrid_y, source_loc,~,estimatedPilotSignals] = mySim.generator.generateMap();
            
            if obj.computeOnlyFeatures ==true
                SqErr_AllRuns_locFree=0;
                SqErr_AllRuns_locBased=0;
                meanErrOnEvalFeat=0;
                locFreeEstimate=0;
                locBasedEstimate=0;
                UEFeatures=0;
                UELocations=0;
                averageMisses=0;
                NMSE_locFree=0;
                NMSE_locBased=0;
                return
            end
            % [trueMap, evaluationGrid_x, evaluationGrid_y] = mySim.generator.generateMap(); %#ok<ASGLU>
            fprintf('Will run a total of %d simulations with %d Montecarlo runs each.\n', length(obj.N_feat)*length(obj.N_ue), obj.n_runs)
            for  ind_feat=1:length(obj.N_feat)
                if obj.PCAenabler==true
                    mySim.PCAdimension=obj.N_feat(ind_feat);
                elseif obj.considerMissingData==true
                    mySim.locFreeEstimator.receiverSensitivity=obj.N_feat(ind_feat);
                elseif obj.testDiffSetOfFeatures==true
                    mySim.features=obj.featureCell{obj.N_feat_used(ind_feat)};
                end
                for ind_nues=1:length(obj.N_ue)
                    nues = obj.N_ue(ind_nues);
                    if obj.inParallel
                        parfor ind_runs=1:obj.n_runs
                            [SqErr_AllRuns_locFree(ind_nues, ind_runs, ind_feat), SqErr_AllRuns_locBased(ind_nues, ind_runs, ind_feat), ...
                                meanErrOnEvalFeat{ind_nues, ind_runs, ind_feat}, locFreeEstimate{ind_nues, ind_runs, ind_feat},UEFeatures{ind_nues, ind_runs, ind_feat},...
                                locBasedEstimate{ind_nues, ind_runs, ind_feat},UELocations{ind_nues, ind_runs, ind_feat},averageMisses{ind_nues, ind_runs, ind_feat}] = ... % SqErr_AllRuns_sampleAverage(ind_nues, ind_runs, ind_feat) = ...
                                mySim.simulate(nues, allPointEstimatedLocation,power_out,trueMap,...
                                evaluationGrid_x, evaluationGrid_y, source_loc, allLocFreeFeatures, estimatedPilotSignals);
                        end
                    else
                        disp 'Not running in parallel.'
                        for ind_runs=1:obj.n_runs
                            [SqErr_AllRuns_locFree(ind_nues, ind_runs, ind_feat), SqErr_AllRuns_locBased(ind_nues, ind_runs, ind_feat), ...
                                meanErrOnEvalFeat{ind_nues, ind_runs, ind_feat}, locFreeEstimate{ind_nues, ind_runs, ind_feat},UEFeatures{ind_nues, ind_runs, ind_feat},...
                                locBasedEstimate{ind_nues, ind_runs, ind_feat},UELocations{ind_nues, ind_runs, ind_feat},averageMisses{ind_nues, ind_runs, ind_feat}] = ... % SqErr_AllRuns_sampleAverage(ind_nues, ind_runs, ind_feat) = ...
                                mySim.simulate(nues, allPointEstimatedLocation,power_out,trueMap,...
                                evaluationGrid_x, evaluationGrid_y, source_loc, allLocFreeFeatures, estimatedPilotSignals);
                        end
                    end
                    fprintf('Finished simulation for N=%d, M=%d\n', nues, obj.N_feat(ind_feat));
                end
            end
            trueMapAverage = mean(mean(trueMap));
            NMSE_locFree = mean(SqErr_AllRuns_locFree ,2)./sum(sum((trueMap-trueMapAverage).^2)); %#ok<NOPRT>
            NMSE_locBased= (mean(SqErr_AllRuns_locBased,2)./(sum(sum((trueMap-trueMapAverage).^2)))); %#ok<NOPRT>
            beep
        end
    end %methods
end

