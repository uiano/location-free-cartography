classdef MatrixCompletor
	% Implementation of the paper "An Accelerated Gradient Method
	% for Trace Norm Minimization" [Ji & Ye, 2009]
	
	% Developer: Luismi
	
	properties
		L_0
		gamma
		mask
		lambda
		W_0
		
		maxIter = 300;
		maxIter_inner = 300;
		epsilon = 1e-5;
	end
	
	methods
		
		function Z_out = p(obj, mu_in, Y_in)
			W_previous = Y_in;
			m_C = W_previous - 1/mu_in*obj.grad_f(W_previous);
			[U, S, V] = svd(m_C);
			S_lambda = S; 
			S_lambda(:, 1:size(S,1)) = ...
				diag(max(0, diag(S)-obj.lambda/mu_in));
			Z_out = U*S_lambda*V';
		end
		
		function out = F(obj, W_in)
			nuclear_norm = sum(svd(W_in));
			out = obj.f(W_in) + obj.lambda*nuclear_norm;
		end
		
		function out = Q(obj, mu_in, X, Y)
			nuclear_norm = sum(svd(X));
			out = obj.P(mu_in, X, Y) + obj.lambda*nuclear_norm;
		end
		
		function out = P(obj, mu_in, W, W_previous)
			delta_W = W-W_previous;
			out = obj.f(W) + trace(delta_W'*obj.grad_f(W_previous)) + mu_in/2*...
				norm(delta_W, 'fro')^2;
		end
		
		function out = f(obj, W)
			delta_W = W-obj.W_0;
			out = 1/2*norm(obj.mask.*delta_W, 'fro').^2;
		end
		
		function out = grad_f(obj, W)
			delta_W = W - obj.W_0;
			out = obj.mask .*delta_W; 
		end
		
		function [L_out, W_out] = ...
				extendedGradientIteration(obj, L_in, W_in)
			% Single iterate of algorithm 1 in [Ji & Ye, 2009]
			L_bar      = L_in; % set L_bar = L_{k-1}
			W_previous = W_in; % W_{k-1}
			
			for i = 1:obj.maxIter_inner
				p = obj.p( L_bar, W_previous);
				F = obj.F(p);
				Q = obj.Q(L_bar, p, W_previous);
				if F <= Q
					break
				end
				L_bar = obj.gamma*L_bar;
			end			
			
			L_now = L_bar; % set L_k = L_bar and update
			W_now = obj.p(L_now, W_previous);
			
			W_out     = W_now;
			L_out     = L_now;
		end
		
		function [L_out, W_out, alpha_out, Z_out] = ...
				acceleratedGradientIteration(obj, ...
				L_in, W_in, alpha_in, Z_in)
			% Single iterate of algorithm 2 in [Ji & Ye, 2009]
			
			L_bar      = L_in; % set L_bar = L_{k-1}
			Z_previous = Z_in; %Z_{k-1}
			alpha_now  = alpha_in; %alpha_k
			W_previous = W_in; % W_{k-1}
			
			broke = 0;
			for i = 1:obj.maxIter_inner
				p = obj.p( L_bar, Z_previous);
				F = obj.F(p);
				Q = obj.Q(L_bar, p, Z_previous);
				if F <= Q
					broke = 1;
					break
				end
				L_bar = obj.gamma*L_bar;
			end
			if not(broke)
				warning('could not find a proper L')
			end
			
			% PROBABLY A TYPO IN THE PAPER:
			% W IS UPDATED USING Z_NOW, BUT THE CHOICE OF L IS DONE
			% USING Z_PREVIOUS.
			% THE PAPER DOES NOT SPECIFY HOW TO GET Z_NOW 
			% FROM Z_PREVIOUS.
			% IS IT THE SAME VARIABLE IN BOTH UPDATES?
			Z_now = Z_previous;
			
			
			L_now = L_bar; % set L_k = L_bar and update
			W_prime = obj.p(L_now, Z_now);
			
% 			Fs = [obj.F(W_prime), obj.F(Z_now), obj.F(W_previous)];
% 			[~, arg] = min(Fs);
% 			switch arg
% 				case 1
% 					W_now = W_prime;
% 				case 2 
% 					W_now = Z_now;
% 				case 3 
% 					W_now = W_previous;
% 			end
            W_now = W_prime;
			
			alpha_next = (1+sqrt(1+4*alpha_now^2))/2;
			% alpha_{k+1} = (1+sqrt(1+4*alpha_{k}^2))/2;
			Z_next = W_now + (alpha_now - 1) / (alpha_next) ...
				* (W_now - W_previous);

			W_out     = W_now;
			L_out     = L_now;
			alpha_out = alpha_next;
			Z_out     = Z_next;
		end
		
		function [W_out, history_out] = extendedGradientAlgorithm (obj)
			% Algorithm 1 in [Ji & Ye, 2009]

			assert(obj.gamma>1);
			L_in = obj.L_0;
			W_in = obj.W_0;
			fval = obj.F(W_in);
			history(obj.maxIter+1) = struct('L', L_in, 'W', W_in, ...
				'fval', fval);
				%'Z', Z_in, 'alpha', alpha_in, 
			history(1) = history(obj.maxIter+1);
			for i = 1:obj.maxIter
				fval_old = fval;
				[L_out, W_out] = ...
					extendedGradientIteration(obj, ...
					L_in, W_in);
				[L_in, W_in] = ...
					deal(L_out, W_out);
				fval = obj.F(W_out);
				if abs(fval_old -fval)<obj.epsilon
					% the chosen stop criterion comes from Yves' script.
					break
				end
				history(i+1) = struct('L', L_in, 'W', W_in, ...
				'fval', fval);
			end
%  			fprintf('Took %d iterations.\n', i)
			history_out = history(1:i);
		end

		
		function [W_out, history_out] = acceleratedGradientAlgorithm (obj)
			% Algorithm 2 in [Ji & Ye, 2009]

			assert(obj.gamma>1);
			L_in = obj.L_0;
			W_in = obj.W_0;
			Z_in = W_in;
			alpha_in = 1;
			fval = 1e10;
			history(obj.maxIter+1) = struct('L', L_in, 'W', W_in, ...
				'Z', Z_in, 'alpha', alpha_in, 'fval', 0);
			history(1) = history(obj.maxIter+1);
			for i = 1:obj.maxIter
				fval_old = fval;
				[L_out, W_out, alpha_out, Z_out] = ...
					acceleratedGradientIteration(obj, ...
					L_in, W_in, alpha_in, Z_in);
				[L_in, W_in, Z_in, alpha_in] = ...
					deal(L_out, W_out, Z_out, alpha_out);
				fval = obj.F(W_out);
				if abs(fval_old -fval)<obj.epsilon
					% the chosen stop criterion comes from Yves' script.
					break
				end
				history(i+1) = struct('L', L_in, 'W', W_in, ...
				'Z', Z_in, 'alpha', alpha_in, 'fval', fval);
			end
 			fprintf('Took %d iterations.\n', i)
			history_out = history(1:i);
		end
	end
end