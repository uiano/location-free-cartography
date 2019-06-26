classdef LocationEstimator
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Xenb
        locationNoiseSTD
    end
    
    methods
        function [estimatedLocationRemin,estimatedLocationAnnih] = estimateLocationReminAndAnnih(obj,measurements) % Reminimization and annihilation methods are used here
            dim=2; % consider sources and UEs in 2D
            n_ues=size(measurements,2);
            n_sources_tdoa = size( measurements,1)-1; % excludes the reference source
            ue_loc_calc_reminimization=zeros(dim,n_ues);
            ue_loc_calc_annihilationRange=zeros(dim,n_ues);
            S=zeros(n_sources_tdoa,dim);
            T=zeros(n_sources_tdoa,1);
            delta=zeros(n_sources_tdoa, n_ues);
            Z=circshift(eye(n_sources_tdoa),-1,1);
            for i=1:n_ues
                for ind2_s=1:n_sources_tdoa
                    S(ind2_s,:)=(obj.Xenb(:,ind2_s+1)-obj.Xenb(:,1))';
                    delta(ind2_s,i)=norm(obj.Xenb(:,ind2_s+1)-obj.Xenb(:,1)).^2-measurements(ind2_s,i).^2;
                    T(ind2_s)=obj.Xenb(1,1)*(obj.Xenb(1,1)-obj.Xenb(1,ind2_s+1))+obj.Xenb(dim,1)*(obj.Xenb(dim,1)-obj.Xenb(dim,ind2_s+1));
                end
                ue_loc_calc_reminimization(:,i)=1/2*pinv(S)*(eye(n_sources_tdoa)-((measurements(1:n_sources_tdoa,i)*measurements(1:n_sources_tdoa,i)'*...
                    (eye(n_sources_tdoa)-S*pinv(S)))./(2* measurements(1:n_sources_tdoa,i)'*(eye(n_sources_tdoa)-S*pinv(S))*...
                    measurements(1:n_sources_tdoa,i))))*(delta(:,i)-2*T);
                
                P=(eye(n_sources_tdoa)-Z)/(diag(measurements(1:n_sources_tdoa,i)));
                ue_loc_calc_annihilationRange(:,i)=1/2*((S'*P'*P*S)\(P*S)')*(P*delta(:,i)-2*P*T);
            end
            estimatedLocationRemin=ue_loc_calc_reminimization;
            estimatedLocationAnnih=ue_loc_calc_annihilationRange;
        end
        
        function estimatedLocation = estimateLocationTreatRangeAs(obj,measurements) % Consider range as a variable method
            dim=2; % consider sources and UEs in 2D
            n_ues=size(measurements,2);
            n_sources_tdoa = size( measurements,1)-1; % excludes the reference source
            ue_loc_calc_TreatRangeAs=zeros(dim,n_ues);
            S_prime=zeros(n_sources_tdoa,dim+1,n_ues);
            T=zeros(n_sources_tdoa,1);
            delta=zeros(n_sources_tdoa, n_ues);
            G=[eye(2);ones(1,2)];
            for i=1:n_ues
                for ind2_s=1:n_sources_tdoa
                    S_prime(ind2_s,:,i)=[(obj.Xenb(:,ind2_s+1)-obj.Xenb(:,1))' measurements(ind2_s,i)];
                    delta(ind2_s,i)=norm(obj.Xenb(:,ind2_s+1)-obj.Xenb(:,1)).^2-measurements(ind2_s,i).^2;
                    T(ind2_s)=obj.Xenb(1,1)*(obj.Xenb(1,1)-obj.Xenb(1,ind2_s+1))+obj.Xenb(dim,1)*(obj.Xenb(dim,1)-obj.Xenb(dim,ind2_s+1));
                end
                xs_hat=1/2*pinv(S_prime(:,:,i))*(delta(:,i)-2*T);
                h=[(xs_hat(1:2)-obj.Xenb(:,1)).^2 ;xs_hat(3)^2];
                C_eInv=(1/(4*obj.locationNoiseSTD.^2))*diag(1./[(xs_hat(1:2)-obj.Xenb(:,1)).^2 ;xs_hat(3)^2]);
                xs_hat_prime=(G'* C_eInv*G)\(G'* C_eInv*h);
                ue_loc_calc_TreatRangeAs(:,i)= sign(xs_hat(1:2)).*sqrt(xs_hat_prime)+obj.Xenb(:,1);
            end
            estimatedLocation= ue_loc_calc_TreatRangeAs; %ue_loc_calc_reminimization
        end
        
        function estimatedLocationDirecCosines = estimateLocationDirectionCosines(obj,measurements) % using direction cosines
            dim=2; % consider sources and UEs in 2D
            n_ues=size(measurements,2);
            n_sources_tdoa = size( measurements,1)-1; % excludes the reference source
            ue_loc_calc_DirecCos=zeros(dim,n_ues);
            centerOfCoord=[16.5;12.5]; % Center of coordinates with x_limts (-5,38) and y_limits (-3,28) or boundaries of the considered area.
            A=zeros(n_sources_tdoa,4, n_ues);
            B=zeros(n_sources_tdoa, n_ues);
            for i=1:n_ues
                for ind2_s=1:n_sources_tdoa
                    A(ind2_s,:,i)=[-(obj.Xenb(:,ind2_s+1)-centerOfCoord-(obj.Xenb(:,1)-centerOfCoord))',...
                        norm(obj.Xenb(:,ind2_s+1)-centerOfCoord-(obj.Xenb(:,1)-centerOfCoord)).^2, -(measurements(ind2_s,i)/3e8)^2];
                    B(ind2_s,i)= measurements(ind2_s,i)/3e8;
                end
                x=(A(:,:,i)'*A(:,:,i))\(A(:,:,i)'*B(:,i));
                ue_loc_calc_DirecCos(:,i)=x(1:2)./(2*x(3))+(obj.Xenb(:,1)-centerOfCoord);
            end
            estimatedLocationDirecCosines=ue_loc_calc_DirecCos;
        end
        
        function estimatedLocationIRWSRDLS = estimateLocationIRWSRDLS(obj,measurements) % using IRWSRDLS
            dim=2;
            n_ues=size(measurements,2);
            n_sources_tdoa = size(measurements,1)-1; % excludes the reference source
            Iter_max=20;
            eps_wls=1e-9;
            eps_bisec=1e-9;
            ue_loc_calc_IRWSRDLS=zeros(dim,n_ues);
            for i=1:n_ues
                ue_loc_calc_IRWSRDLS(:,i)=LocationEstimator.estimateLocationIRWSRDLS_One(obj.Xenb(:,2:n_sources_tdoa+1)...
                    ,measurements(1:n_sources_tdoa,i),Iter_max,eps_wls,eps_bisec);
            end
            estimatedLocationIRWSRDLS=ue_loc_calc_IRWSRDLS;
        end
    end
    
    methods (Static)
        function x = estimateLocationIRWSRDLS_One(Am,dn,Kmax,epsi1,epsi)
            xk = LocationEstimator.srd_ls(Am,dn,epsi);
            C = [eye(2) zeros(2,1); zeros(1,2) -1];
            m = size(Am,2);
            err = 10; %Starting error
            k = 1;
            wk_0=ones(1,m);
            while err >= epsi1 && k <= Kmax
         
                wk=wk_0(k,:);
                Wk = diag(wk);
                B = Wk*[-2*Am' -2*dn];
                g = zeros(m,1);
                for i = 1:m
                    g(i) = wk(i)*(dn(i)^2 - norm(Am(:,i))^2);
                end
                btg = B'*g;
                B1 = B'*B;
                [P,D,~] = eig(C,B1);
                L =sort(diag(D),'descend');
                alp =-1./L;
                alp0 = alp(3);
                alp1 = alp(2);
                alp2 = alp(1);
                aU = alp0;
                aL = alp1;
                dt = aU - aL;
                while dt > epsi
                    lam = 0.5*(aL + aU);
                    yt = (B1+lam*C)\btg;
                    phit = yt'*C*yt;
                    if phit > 0
                        aL = lam;
                    else
                        aU = lam;
                    end
                    dt = aU - aL;
                end
                lam = 0.5*(aL + aU);
                yt = (B1+lam*C)\btg;
                yn1 = yt(3);
                if yn1 >= 0
                    xk_new = yt(1:2);
                else
                    f = P'*btg;
                    f1 = f(1); f2 = f(2); f3 = f(3);
                    f1s = f1^2; f2s = f2^2; f3s = f3^2;
                    d1 = D(1,1); d1s = d1^2;
                    d2 = D(2,2); d2s = d2^2;
                    d3 = D(3,3); d3s = d3^2;
                    a0 = f1s*d1 + f2s*d2 + f3s*d3;
                    a1 = 2*f1s*d1*(d2+d3) + 2*f2s*d2*(d1+d3) + 2*f3s*d3*(d1+d2);
                    a2 = f1s*d1*(d2s+d3s+4*d2*d3) + f2s*d2*(d1s+d3s+4*d1*d3) ...
                        + f3s*d3*(d1s+d2s+4*d1*d2);
                    a3 = 2*f1s*d1*(d2s*d3+d3s*d2) + 2*f2s*d2*(d1s*d3+d3s*d1) ...
                        + 2*f3s*d3*(d1s*d2+d2s*d1);
                    a4 = f1s*d1*d2s*d3s + f2s*d2*d1s*d3s + f3s*d3*d1s*d2s;
                    b_lu = [a4 a3 a2 a1 a0];
                    rts = roots(b_lu);
                    lamq = [];
                    for i = 1:4
                        if imag(rts(i)) == 0
                            lamq = [lamq real(rts(i))];
                        end
                    end
                    I02 = [];
                    I1 = [];
                    Lq = length(lamq);
                    for i = 1:Lq
                        ti = lamq(i);
                        if ti > alp0
                            I02 = [I02 ti];
                        elseif ti < alp0 && ti > alp1
                            I1 = [I1 0];
                        elseif ti < alp1 && ti > alp2
                            I02 = [I02 ti];
                        end
                    end
                    L02 = length(I02);
                    Yt2 = zeros(3,1);
                    for i = 1:L02
                        yi = (B1+I02(i)*C)\btg;
                        if yi(3) >= 0
                            Yt2 = [Yt2 yi];
                        end
                    end
                    L2s = size(Yt2,2);
                    obj2 = zeros(L2s,1);
                    for i = 1:L2s
                        obj2(i) = norm(B*Yt2(:,i) - g);
                    end
                    [~,ind2] = min(obj2);
                    yt = Yt2(:,ind2);
                    xk_new = yt(1:2);
                end
                err = norm(xk_new - xk);
                xk = xk_new;
                for i = 1:m
                    wk_new(i) =1/abs(dn(i)+ norm(xk)+norm(xk-Am(:,i)));
                end
                wk_0(k+1,:)=wk_new;
                k = k + 1;
            end
            x = xk;
        end
        
        
        function x= srd_ls(Am,dn,epsi)
            m = size(Am,2);
            B = [-2*Am' -2*dn];
            C = [eye(2) zeros(2,1); zeros(1,2) -1];
            B1 = B'*B;
            g = zeros(m,1);
            for i = 1:m
                g(i) = dn(i)^2 - norm(Am(:,i))^2;
            end
            btg = B'*g;
            [P,D,~] = eig(C,B1);
            L =sort(diag(D),'descend');
            alp =-1./L;
            alp0 = alp(3);
            alp1 = alp(2);
            alp2 = alp(1);
            aU = alp0;
            aL = alp1;
            dt = aU - aL;
            while dt > epsi
                lam = 0.5*(aL + aU);
                yt = (B1+lam*C)\btg;
                phit = yt'*C*yt;
                if phit > 0
                    aL = lam;
                else
                    aU = lam;
                end
                dt = aU - aL;
            end
            lam = 0.5*(aL + aU);
            yt = (B1+lam*C)\btg;
            yn1 = yt(3);
            if yn1 >= 0
                x = yt(1:2);
                i1 = 1;
            else
                f = P'*btg;
                f1 = f(1); f2 = f(2); f3 = f(3);
                f1s = f1^2; f2s = f2^2; f3s = f3^2;
                d1 = D(1,1); d1s = d1^2;
                d2 = D(2,2); d2s = d2^2;
                d3 = D(3,3); d3s = d3^2;
                a0 = f1s*d1 + f2s*d2 + f3s*d3;
                a1 = 2*f1s*d1*(d2+d3) + 2*f2s*d2*(d1+d3) + 2*f3s*d3*(d1+d2);
                a2 = f1s*d1*(d2s+d3s+4*d2*d3) + f2s*d2*(d1s+d3s+4*d1*d3) + f3s*d3*(d1s+d2s+4*d1*d2);
                a3 = 2*f1s*d1*(d2s*d3+d3s*d2) + 2*f2s*d2*(d1s*d3+d3s*d1) + 2*f3s*d3*(d1s*d2+d2s*d1);
                a4 = f1s*d1*d2s*d3s + f2s*d2*d1s*d3s + f3s*d3*d1s*d2s;
                rts = roots([a4 a3 a2 a1 a0]);
                lamq = [];
                for i = 1:4
                    if imag(rts(i)) == 0
                        lamq = [lamq real(rts(i))];
                    end
                end
                I02 = [];
                I1 = [];
                L = length(lamq);
                for i = 1:L
                    ti = lamq(i);
                    if ti > alp0
                        I02 = [I02 ti];
                    elseif ti < alp0 && ti > alp1
                        I1 = [I1 0];
                    elseif ti < alp1 && ti > alp2
                        I02 = [I02 ti];
                    end
                end
                L02 = length(I02);
                Yt2 = zeros(3,1);
                for i = 1:L02
                    yi = (B1+I02(i)*C)\btg;
                    if yi(3) >= 0
                        Yt2 = [Yt2 yi];
                    end
                end
                L2s = size(Yt2,2);
                obj2 = zeros(L2s,1);
                for i = 1:L2s
                    obj2(i) = norm(B*Yt2(:,i) - g);
                end
                [~,ind2] = min(obj2);
                yt = Yt2(:,ind2);
                x = yt(1:2);
            end
        end
    end
end
