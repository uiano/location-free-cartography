classdef FeatureExtractor
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    properties (Constant)
        c = 3e8;
    end
    properties
        sampling_period
    end
    
    methods
        function  estimatedCentersOfMass  = locFreeExtract(obj,pilotSignals)
            %UNTITLED4 Summary of this function goes here
            %   Detailed explanation goes here
            n_ues=size(pilotSignals,3);
            n_sources=size(pilotSignals,1);
            estimatedCentersOfMass=zeros(nchoosek(n_sources,2),n_ues);
            combin_sources=combnk(1:n_sources,2);
            if n_sources < 6
                combin_sources=flipud(combin_sources);
            end
            for ind_ues=1:n_ues
                allCenterOfMass=zeros(1,nchoosek(n_sources,2));
                for ind_comb=1:size(combin_sources,1)
                    [hD_corr, lags]=xcorr(pilotSignals(combin_sources(ind_comb,2),:,ind_ues),pilotSignals(combin_sources(ind_comb,1),:,ind_ues));
                    lag_CenterOfMass=lags*(abs(hD_corr).^2)'/sum(abs(hD_corr).^2);
                    allCenterOfMass(ind_comb)=lag_CenterOfMass*obj.sampling_period*obj.c;
                end
                estimatedCentersOfMass(:,ind_ues)=allCenterOfMass;
            end
        end
        
        
        function estimatedDistances2 = locBasedExtract(obj,pilotSignals)
            %UNTITLED4 Summary of this function goes here
            %   Detailed explanation goes here
            
            n_ues=size(pilotSignals,3);
            n_sources=size(pilotSignals,1);
            estimatedDistances2=zeros(n_sources-1,n_ues);
            for ind_ue=1:n_ues
                averaDistAllSources=zeros(n_sources-1,1);
                for ind_sourceHd=2:n_sources
                    [hD_corr, lags]=xcorr(pilotSignals(ind_sourceHd,:,ind_ue),pilotSignals(1,:,ind_ue));
                    [~, Ind_max]=max(abs(hD_corr));
                    tdoaOneSource=lags(Ind_max)*obj.sampling_period;
                    averaDistAllSources(ind_sourceHd-1)=tdoaOneSource*obj.c;
                end
                estimatedDistances2(:,ind_ue)=averaDistAllSources;
            end
            
        end
    end
    
end

