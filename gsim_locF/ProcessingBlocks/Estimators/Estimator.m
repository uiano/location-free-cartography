classdef Estimator
	%Estimator Superclass of (location-free and location-based) 
	% map estimators.
	
	properties
	end
	
	methods (Abstract)
		train(obj, measurements);
% 		estimate(obj, ...
% 			coefficients, trainingData, evalGrid_x, evalGrid_y);
% 		estimateGivenDistances(obj, ...
% 			coefficients, trainingData, t_distances) ;
	end
	
end

