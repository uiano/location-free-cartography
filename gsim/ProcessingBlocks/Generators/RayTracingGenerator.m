classdef RayTracingGenerator < MapGenerator
    % A class between MapGenerator and all subclasses that
    % generate a "channel-to-map" (C2M) using raytracing.
    
    properties (Constant)
        c = 3e8;
    end
    properties
        f   %carrier frequency
        xt  % x coordinates of the source locations
        yt  % y coordinates of the source locations
        zt  % z coordinates of the source locations
        n_gridpoints_x % size of evaluation grid in x axis
        n_gridpoints_y % size of evaluation grid in y axis
        maxSamplesPerPilot
        
        % sampling period of the receiver (used in computation of the digital impulse response)
        sampling_period
    end
    
    methods (Abstract)
        [distance6ray, powerDelayProfile,estimatedToaRange,totalReceivedPower, delayCenterOfMass] = calculateGivenReceiverLocation(obj, xr, yr);
        % This method is currently implemented by the following subclasses:
        % StreetCanyonGenerator
        % MultiWallGenerator
    end
    
    methods
        function [power_out,map_out, evaluationGrid_x, evaluationGrid_y, source_loc,estimatedDistances,estimatedPilotSignals] = generateMap(obj)
            % GenerateMap Create a C2M by calling a raytracer for every point in grid
            % [map_out, gridX, gridY, source_loc, estimatedDistances,
            % estimatedCentersOfMass] = generateMap(obj) writes the C2M in
            % map_out, the x and y locations in gridX and gridY
            % respectively...
            source_loc=[obj.xt; obj.yt];
            [evaluationGrid_x, evaluationGrid_y] = ndgrid(linspace(obj.x2, obj.x1, obj.n_gridpoints_x), linspace(obj.street_limits(1),obj.street_limits(2), obj.n_gridpoints_y));
            map_out=zeros(obj.n_gridpoints_x,obj.n_gridpoints_y);
            power_out=zeros(size(source_loc,2),obj.n_gridpoints_x,obj.n_gridpoints_y);
            estimatedPilotSignals=zeros(size(source_loc,2),obj.maxSamplesPerPilot,obj.n_gridpoints_x,obj.n_gridpoints_y);
            estimatedDistances=zeros(obj.n_gridpoints_x,obj.n_gridpoints_y,size(source_loc,2)-1);
            
            ltm =LoopTimeControl(obj.n_gridpoints_x*obj.n_gridpoints_x);
            for ind_xr = 1:obj.n_gridpoints_x
                for ind_yr = 1:obj.n_gridpoints_y
                    [distance6ray,  powerDelayProfile,totalReceivedPowerAllSources, totalReceivedPower,estimatedTdoaRange, H_D] = ...
                        obj.calculateGivenReceiverLocation(evaluationGrid_x(ind_xr, ind_yr), evaluationGrid_y(ind_xr,ind_yr));
                    power_out(:,ind_xr, ind_yr)=totalReceivedPowerAllSources;
                    map_out(ind_xr, ind_yr) = totalReceivedPower;
                    estimatedDistances(ind_xr, ind_yr,:) = estimatedTdoaRange;
                    estimatedPilotSignals(:,:,ind_xr, ind_yr) =H_D; % allCenterOfMass;
                    ltm.go(ind_yr+(ind_xr-1)*obj.n_gridpoints_y);
                end
            end
        end
        
        function [h_D,samplingTime] = digitalImpulseResponse(obj, alpha, delays)
            % calculate digital impulse response given the delays and powers of
            % the analog (time-continuous) channel impulse response
            N_rays=length(alpha);  % Set N_rays to 3 while simulating with the multiwall generator; 6 with the streetcanyon generator
            maxmumSamplingTime=ceil(max(delays)/obj.sampling_period); % 10 added is for ensuring some margin, so that all rays are considered
            samplingTime=0:maxmumSamplingTime;
            impulseAllRays=zeros(length(samplingTime),N_rays);
            %sampling frequency in MHz
            for ind_ray=1:N_rays
                impulseAllRays(:,ind_ray)=alpha(ind_ray)*exp(-1i*2*pi*obj.f*delays(ind_ray))*sinc(samplingTime-delays(ind_ray)/obj.sampling_period);
            end
            h_D=sum(impulseAllRays,2);
        end
        
        function ToA = estimateTimeOfArrival(obj, h_D)
            [~,Ind]=max(abs(h_D));
            ToA=Ind*obj.sampling_period;
        end
    end
    
end

