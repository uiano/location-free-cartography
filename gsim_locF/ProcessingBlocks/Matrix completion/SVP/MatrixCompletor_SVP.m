classdef MatrixCompletor_SVP
    % Implementation of the singular value projection algorithm for matrix
    % completion
    
    properties
        desiredRank
        gamma_const
        X_0
        mask
        alpha_init  % Initial step size
        delta_const
        maxIter =500;
        epsilon = 1e-3;
    end
    
    methods
        
        function Z_out = truncated_SVD(obj,X_in,eta)
            m_C =X_in + eta;
            [U, S, V] = svd(m_C);
            S(:,obj.desiredRank+1:end)=0;
            Z_out = U*S*V';
        end
        function out = f(obj, X)
            delta_X = X-obj.X_0;
            out = 1/2*norm(obj.mask.*delta_X, 'fro').^2;
        end
        
        function out = grad_f(obj, W)
            delta_W = W - obj.X_0;
            out = obj.mask .*delta_W;
        end
        
        function [alpha_armijo] = armijo(obj,X_in,direction)
            %% ARMIJO
            % DESCRIPTION: function to check whether the provided steplength satisfies Armijo
            % condition: f(X_in+alpha*d)<=f(X_in)+gamma*alpha*gradient(f(X_in))'*direction
            alpha=obj.alpha_init;
            delta=obj.delta_const;
            gamma=obj.gamma_const;
            j = 1;
            while (j>0)
                x_new = X_in+alpha.*direction;
                if obj.f(x_new)<=obj.f(X_in)+gamma*alpha*trace(direction'*obj.grad_f(X_in))
                    j = 0;
                    alpha_armijo = alpha;
                else
                    alpha = alpha*delta;
                end
            end
            
        end
        
        function [X_out,obj_val] = singularValueProjection(obj)
            X_in = obj.X_0;
            eta=-obj.grad_f(X_in);
            fval = obj.f(X_in);
            obj_val=[];
            for i = 1:obj.maxIter
                fval_old = fval;
                alpha_armijo = obj.armijo(X_in,eta);
                eta=-alpha_armijo*obj.grad_f(X_in);
                X_now=obj.truncated_SVD(X_in , eta);
                fval = obj.f(X_now);
                obj_val=[obj_val,fval] ;
                if abs(fval_old -fval)<obj.epsilon % the chosen stop criterion 
                    break
                end
                
                X_in=X_now;
            end
            X_out=X_now;
%             fprintf('Took %d iterations.\n', i)
        end
    
    end
end