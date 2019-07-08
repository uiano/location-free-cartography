classdef MapGenerator
	% All generator classes in wsim should inherit from MapGenerator.
	
	properties
	end
	
	methods
		function obj = MapGenerator()
		end
	end
	methods (Abstract)
		[map_out, evaluationGrid_x, evaluationGrid_y, source_loc, ...
			estimatedDistances,estimatedCentersOfMass] = generateMap(obj)
		% This abstract class is currently implemented by the subclasses:
		% RayTracingGenerator
		% PrecomputedMapGenerator
	end
	
end

