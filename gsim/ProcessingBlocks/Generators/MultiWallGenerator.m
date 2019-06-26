classdef MultiWallGenerator < RayTracingGenerator
    % Allows to create C2Ms using the multiwall propagation model
    % proposed in the paper by Salaheddin [] et al
    % TODO: [write citation here]
    
    properties
        % Coming from the superclass
        %         c = 3e8;
        %         f
        %         xt   Location of the transmitter (s)
        %         yt
        %         zt
        %         n_gridpoints_x
        %         n_gridpoints_y
        %         sampling_period
        ptx
        zr
        x1
        x2
        %         RG
        %         RW
        estimateLocFreeSpace
        path_loss_exp
        boundary % Defining the boundary of the analysis (something like a boundary condition)
        street_limits
        losFlag = 1;
        optimizationMode = 0;
        reflectionFlag=1;                    % whether or not calculate First reflections
        secondReflectionFlag=1;
        reflectExaggerationFac=1; % Must be 1 unless to emphasize reflection for demonstration purposes                                    % whether or not calculate LoS
        disableIncidentAngle=0;                % 1 Disables the incident angle calculation, if disableIncidentAngle= 1
        solidIncidentAngle=45;                  % if disableIncidentAngle =1, then assign this which overwrites all the incident angles! This is unnecessary feature
        polarizationSwap =1;               % (See notes in HOW THIS WORSK)  % 1, Applies TE to walls and TM to ceiling. 0, applies TM to the walls and TE to the ceiling
        refDistance=1.25;            % Reference distance from Tx, which is approximately 3 wavelength at 800MHz
        FPSLRefLoss=0;
        imageRSSIScale=1;       %  % increase this if number of meshes nodes are small
        
        antennaGainRes=40;
        antennaEffiLoss=0;  %-11.5        % dB antenna efficiency, missmatch and loss all together
        
        demoMode=0;
        
        ceilingEnable=0;  % Allowing to define ceiling and floor
        groundLevel=0;
        ceilingLevel=3;   % Height of the ceiling
        
        X % x coordinates for the walls of the building
        Y % y coordinates for the walls of the building
        
    end
    
    methods
        
        function [distanceNray, powerDelayProfile,totalReceivedPowerAllSources,totalReceivedPower,estimatedTdoaRange, H_D] = calculateGivenReceiverLocation(obj, xr, yr)
            
            n_sources=length(obj.xt);
            averaDistAllSources=zeros(1,n_sources-1);
            %             allCenterOfMass=zeros(1, nchoosek(n_sources,2));
            %             delayCenterOfMass=zeros(1,n_sources-1);
            totalReceivedPowerAllSources=zeros(1,n_sources);
            H_D=zeros(n_sources,obj.maxSamplesPerPilot); % Assuming that we do not have more than... samples of the impulse outputs: padded array
            D_direct=zeros(n_sources,1); % tries the free space propagation in TDOA
            for ind_source=1:n_sources
                %%      This section comes from the multiwall code by  Salaheddin
                
                % Antenna Gain pattern calculation
                TxAntennaGainAE = AntennaTemp (obj.antennaGainRes,obj.demoMode) + obj.antennaEffiLoss;  % TxAntennaGainAE needs to be in dB
                RxAntennaGainAE = TxAntennaGainAE;
                
                % Walls to be defined in a clockwise or counter clockwise manner
                % CLOCK WISE WALL DEFINITION
                
                % Reads the structure from an excel file (see in this code section at the
                % top)
                [wallxyz1, wallxyz2, wallxyz3, wallxyz4,wallX,wallY,wallZ] = CSV23D_V1(obj.X, obj.Y,obj.demoMode,obj.groundLevel,obj.ceilingLevel,[obj.xt(ind_source) obj.yt(ind_source) obj.zt]);
                
                wall.xyz1 = wallxyz1;
                wall.xyz2 = wallxyz2;
                wall.xyz3 = wallxyz3;
                wall.xyz4 = wallxyz4;
                wall.X = wallX;
                wall.Y = wallY;
                wall.Z = wallZ;
                % Define the ceiling of the structure manually if required walls can be
                % defined the same fashion
                if obj.ceilingEnable == 1
                    
                    ceillFloor.xyz1 = [-5,-5,obj.ceilingLevel
                        -5,-5,obj.groundLevel
                        ];
                    
                    ceillFloor.xyz2 = [-5,22.3,obj.ceilingLevel
                        -5,22.3,obj.groundLevel
                        
                        ];
                    
                    ceillFloor.xyz3 = [36.9,22.3,obj.ceilingLevel
                        36.9,22.3,obj.groundLevel
                        ];
                    
                    ceillFloor.xyz4 = [36.9,-5,obj.ceilingLevel
                        36.9,-5,obj.groundLevel
                        ];
                else
                    
                    ceillFloor.xyz1 = [];
                    ceillFloor.xyz2 = [];
                    ceillFloor.xyz3 = [];
                    ceillFloor.xyz4 = [];
                    
                end
                % Adding Ceillilng and Floor to the structure
                for i = 1:size(ceillFloor.xyz1,1)
                    
                    wall.xyz1 = [wall.xyz1;ceillFloor.xyz1(i,:)];
                    wall.xyz2 = [wall.xyz2;ceillFloor.xyz2(i,:)];
                    wall.xyz3 = [wall.xyz3;ceillFloor.xyz3(i,:)];
                    wall.xyz4 = [wall.xyz4;ceillFloor.xyz4(i,:)];
                    
                end
                wallRelativePerm = 20*ones(size(wall.xyz1,1),1); % +size(ceillFloor.xyz1,1) TO DO:  Assign the permittivity manually.
                % Producing Plan as Image
                % finding linear mapping that maps:
                % boundary(1,1) ---->  1
                % boundary(1,2) ---->  mesh.xNodeNum .* imageRSSIScale
                % map = alpha .* (x) + beta, where alpha is the expantion fac & beta is shift
                alpha.x = (1 - obj.n_gridpoints_x.* obj.imageRSSIScale)./(obj.boundary(1,1) - obj.boundary(1,2));
                beta.x  = 1 - (obj.boundary(1,1) .* alpha.x);
                
                alpha.y = (1 - obj.n_gridpoints_y.* obj.imageRSSIScale)./(obj.boundary(2,1) - obj.boundary(2,2));
                beta.y  = 1 - (obj.boundary(2,1) .* alpha.y);
                
                
                imageBoundary(1,:) = round(alpha.x .* obj.boundary(1,:) + beta.x);
                imageBoundary(2,:) = round(alpha.y .* obj.boundary(2,:) + beta.y);
                
                % Shifting all the walls and Tx position
                imageWalls.x = (reshape(round((alpha.x .* [wall.xyz1(:,1); wall.xyz2(:,1);wall.xyz3(:,1);wall.xyz4(:,1)]) + beta.x),size(wall.xyz1,1),4));
                imageWalls.y = (reshape(round((alpha.y .* [wall.xyz1(:,2); wall.xyz2(:,2);wall.xyz3(:,2);wall.xyz4(:,2)]) + beta.y),size(wall.xyz1,1),4));
                imageTx.xy(:,1) = round((alpha.x .* obj.xt(ind_source)) + beta.x);
                imageTx.xy(:,2) = round((alpha.y .* obj.yt(ind_source)) + beta.y);
                
                structImage = false(imageBoundary(2,2),imageBoundary(1,2));
                
                if exist('measurementLocationXyz','var')
                    measLoc.xy(:,1) = round((alpha.x .* measurementLocationXyz(:,1)) + beta.x);
                    measLoc.xy(:,2) = round((alpha.y .* measurementLocationXyz(:,2)) + beta.y);
                    for i = 1:size(measurementLocationXyz,1)
                        structImage(measLoc.xy(i,1),measLoc.xy(i,2)) = 1;
                    end
                end
                
                
                for i = 1:size([obj.xt(ind_source) obj.yt(ind_source) obj.zt],1)
                    structImage(imageTx.xy(i,1),imageTx.xy(i,2)) = 1;
                end
                
                % structImage =  imdilate(structImage,strel('disk',3));
                
                if exist('MATLABStructMode','var')
                    if MATLABStructMode == 1
                        for j = 1:(size(wall.xyz1,1) - size(ceillFloor.xyz1,1))
                            for i = 1:3
                                [wallC,wallR] = bresenham(imageWalls.x(j,i),imageWalls.y(j,i),imageWalls.x(j,i+1),imageWalls.y(j,i+1)); %LOS between Tx &Rx
                                for k = 1:numel(wallC)
                                    structImage(wallC(k),wallR(k)) = 1;
                                end
                            end
                        end
                    end
                else
                    for j = 1:size(wall.xyz1,1)
                        for i = 1:3
                            [wallC,wallR] = bresenham(imageWalls.x(j,i),imageWalls.y(j,i),imageWalls.x(j,i+1),imageWalls.y(j,i+1)); %LOS between Tx &Rx
                            for k = 1:numel(wallC)
                                structImage(wallC(k),wallR(k)) = 1;
                            end
                        end
                    end
                end
                
                
                structImage = structImage(1:imageBoundary(1,2),1:imageBoundary(2,2));
                structImage = imrotate(structImage,90);
                
                
                if obj.demoMode ==1
                    figure
                    imshow(structImage)
                end
                
                %% Defining a Finite Panel (wall)
                
                for i = 1:size(wall.xyz1,1)
                    wall.minMax.x(i,:) = [min([wall.xyz1(i,1),wall.xyz2(i,1),wall.xyz3(i,1),wall.xyz4(i,1)]),...
                        max([wall.xyz1(i,1),wall.xyz2(i,1),wall.xyz3(i,1),wall.xyz4(i,1)])];
                    wall.minMax.y(i,:) = [min([wall.xyz1(i,2),wall.xyz2(i,2),wall.xyz3(i,2),wall.xyz4(i,2)]),...
                        max([wall.xyz1(i,2),wall.xyz2(i,2),wall.xyz3(i,2),wall.xyz4(i,2)])];
                    wall.minMax.z(i,:) = [min([wall.xyz1(i,3),wall.xyz2(i,3),wall.xyz3(i,3),wall.xyz4(i,3)]),...
                        max([wall.xyz1(i,3),wall.xyz2(i,3),wall.xyz3(i,3),wall.xyz4(i,3)])];
                end
                
                %% 3D Formation of the Structure
                if obj.demoMode == 1
                    figure
                    wall.X = [wall.xyz1(:,1)';wall.xyz2(:,1)';wall.xyz3(:,1)';wall.xyz4(:,1)'];
                    wall.Y = [wall.xyz1(:,2)';wall.xyz2(:,2)';wall.xyz3(:,2)';wall.xyz4(:,2)'];
                    wall.Z = [wall.xyz1(:,3)';wall.xyz2(:,3)';wall.xyz3(:,3)';wall.xyz4(:,3)'];
                    wall.C = zeros(size(wall.X));
                    fill3(wall.X, wall.Y, wall.Z,wall.C)
                    hold on
                    for i = 1:size([obj.xt(ind_source) obj.yt(ind_source) obj.zt],1)
                        plot3(Tx.xyz(i,1),Tx.xyz(i,2),Tx.xyz(i,3),'LineStyle','none','Marker','*','Color','Red');
                    end
                    %             pause(eps) % to show the 3D structure
                    xlabel('X')
                    ylabel('Y')
                    zlabel('Z')
                end
                
                if ~isequal(size(wall.xyz1,1), size(wallRelativePerm,1))
                    disp(sprintf(['Error: Geometry & setting dimension mismatch!\nsize(wall.xyz1)=[',num2str(size(wall.xyz1)),...
                        ']\nsize(wallRelativePerm)=[',num2str(size(wallRelativePerm)),']\nFirst dimmensions must be of the same size']))
                    return
                end
                % end
                
                
                %% Claculating Fresnel Coefficients for Walls
                
                for i = 1:size(wall.xyz1,1)
                    
                    [wall.TE.refFac(i,:),wall.TE.transFac(i,:),wall.TM.refFac(i,:),wall.TM.transFac(i,:)] = ...
                        FresnelCoefficients(1,wallRelativePerm(i,1),0:90,0);
                    
                    %     disableIncidentAngle =0; % 1 is disabling incident angle.
                    
                    
                    if obj.disableIncidentAngle == 1
                        %         obj.solidIncidentAngle = 45; % will use this angle for incidence instead
                        %         disp('Angle of incidence is disabled!')
                        wall.TE.refFac(i,:) = repmat(wall.TE.refFac(i,obj.solidIncidentAngle),1,91);
                        wall.TE.transFac(i,:) = repmat(wall.TE.transFac(i,obj.solidIncidentAngle),1,91);
                        wall.TM.refFac(i,:) = repmat(wall.TM.refFac(i,obj.solidIncidentAngle),1,91);
                        wall.TM.transFac(i,:) = repmat(wall.TM.transFac(i,obj.solidIncidentAngle),1,91);
                    end
                    
                    
                end
                if obj.disableIncidentAngle == 1
                    %         obj.solidIncidentAngle = 45; % will use this angle for incidence instead
                    disp('Angle of incidence is disabled!')
                end
                
                
                if obj.optimizationMode == 1
                    Rx.xyz = Rx.xyz(locationIndex,:);
                else
                    %     Rx.xyz = [reshape(X,[],1,1),reshape(Y,[],1,1),reshape(Z,[],1,1)];
                    Rx.xyz =[xr,yr,obj.zr];
                end
                
                
                if obj.demoMode == 1
                    plot3(Rx.xyz(:,1),Rx.xyz(:,2),Rx.xyz(:,3),'LineStyle','none','Marker','.');
                    view(2)
                end
                
                
                % Distance Of TX(s) From Every Mesh Node (RXi), Its vector and unit vector
                for i = 1:size([obj.xt(ind_source) obj.yt(ind_source) obj.zt],1)
                    RxTx.vec.xyz(:,1:3,i) = repmat([obj.xt(ind_source) obj.yt(ind_source) obj.zt],size(Rx.xyz,1),1) - Rx.xyz;
                    RxTx.dist(:,1,i) = sqrt(sum(RxTx.vec.xyz(:,1:3,i).^2,2));
                    RxTx.unitVec.xyz(:,1:3,i) = RxTx.vec.xyz(:,1:3,i) ./ repmat(RxTx.dist(:,1,i),1,3);
                end
                
                
                %% EQUATING THE PANELS (WALLS) IN 3D
                % Find Walls Normals
                wall.normal.xyz = (cross(wall.xyz2 - wall.xyz1,wall.xyz3 - wall.xyz1,2));
                wall.unitNormal.xyz = wall.normal.xyz ./ repmat(sqrt(sum(wall.normal.xyz.^2,2)),1,3);
                
                
                
                for i = 1:size([obj.xt(ind_source) obj.yt(ind_source) obj.zt],1)
                    % Finding Projection of Tx on each panel https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
                    Tx.wallProj.xyz(:,:,i) = repmat((dot((wall.xyz1 - repmat([obj.xt(ind_source) obj.yt(ind_source) obj.zt],size(wall.xyz1,1),1)),wall.unitNormal.xyz,2)...
                        ./dot(wall.unitNormal.xyz,wall.unitNormal.xyz,2)),1,3).* wall.unitNormal.xyz + repmat([obj.xt(ind_source) obj.yt(ind_source) obj.zt],size(wall.unitNormal.xyz,1),1);
                    % Calculating the reflection (mirror) of Tx accross each panel
                    Tx.wallReflec.xyz(:,:,i) = repmat([obj.xt(ind_source) obj.yt(ind_source) obj.zt],size(wall.unitNormal.xyz,1),1) + 2.* (Tx.wallProj.xyz(:,:,i)...
                        - repmat([obj.xt(ind_source) obj.yt(ind_source) obj.zt],size(wall.unitNormal.xyz,1),1));
                end
                
                
                %% Calculating the Second Image of Tx across each wall (every Tx.WallReflec should be images across all walls).
                
                for i = 1:size(wall.xyz1,1) % only works for first Tx
                    
                    % Tx.secondProjWallj.xyz(:,:,x) contains the projections of wall x-th Tx.wallReflec.xyz (x,:,1)
                    Tx.secondProjWallj.xyz(:,:,i) = repmat((dot((wall.xyz1 - repmat(Tx.wallReflec.xyz(i,:,1),size(wall.xyz1,1),1)),wall.unitNormal.xyz,2)...
                        ./dot(wall.unitNormal.xyz,wall.unitNormal.xyz,2)),1,3).* wall.unitNormal.xyz + repmat(Tx.wallReflec.xyz(i,:,1),size(wall.unitNormal.xyz,1),1);
                    % only works for 1 TX so (:,:,1)
                    
                    Tx.secondReflecWallj.xyz(:,:,i) = repmat(Tx.wallReflec.xyz(i,:,1),size(wall.unitNormal.xyz,1),1) + 2.* (Tx.secondProjWallj.xyz(:,:,i)...
                        - repmat(Tx.wallReflec.xyz(i,:,1),size(wall.unitNormal.xyz,1),1));
                    
                end
                %%
                
                % Line Vector Between TxReflection & Rx
                for i = 1:size([obj.xt(ind_source) obj.yt(ind_source) obj.zt],1) %
                    for j = 1:size(Tx.wallReflec.xyz,1)
                        % 4th dimension represents the Tx
                        Rx2TxRefl.vec.xyz(:,1:3,j,i) = repmat(Tx.wallReflec.xyz(j,:,i),size(Rx.xyz,1),1) - Rx.xyz;
                        % 3rd dimension represents Tximage across the wall which the reflection took place
                        Rx2TxRefl.dist(:,1,j,i) = sqrt(sum(Rx2TxRefl.vec.xyz(:,1:3,j,i).^2,2));
                        Rx2TxRefl.unitVec.xyz(:,1:3,j,i) = Rx2TxRefl.vec.xyz(:,1:3,j,i) ./ repmat(Rx2TxRefl.dist(:,1,j,i),1,3);
                        %         Rx2TxRefl.dist(:,1 = sqrt(sum(Rx2TxRefl.xyz.^2,2));
                    end
                end
                
                
                %% Calculating LOS Component
                if obj.losFlag == 1
                    
                    for i = 1:1
                        % beam angle is measured in relation to the origin, unit vectors of
                        % [1,0,0],[0,1,0],[0,0,1]. This depends on the orientation of the TX
                        losBeamAngle.Tx.Ele(:,i) = asind(-RxTx.vec.xyz(:,3,i) ./ sqrt(sum(RxTx.vec.xyz(:,:,i).^2,2)));  % Elevation angle (between beam and Z plane not it's normal) -90<ele<90 degrees
                        losBeamAngle.Tx.Azi(:,i) = atan2(-RxTx.vec.xyz(:,2,i),-RxTx.vec.xyz(:,1,i)) * (180/pi); % Azimuth angle (between x and beam) -180<azi<180
                        
                        losBeamAngle.Tx.Ele(find(isnan(losBeamAngle.Tx.Ele(:,i)) == 1),i) = 0; % if nan turns the beam angle to 0
                        losBeamAngle.Tx.Zen(:,i) = abs(90-losBeamAngle.Tx.Ele(:,i));  % Zenith angle is calculated and used to find the antenna gain
                        losBeamAngle.Tx.Azi(find(isnan(losBeamAngle.Tx.Azi(:,i)) == 1),i) = 0; % Azimuth angle is unlikely to be nan due to use of atan2
                        losBeamAngle.Tx.ZenIndex(:,i) = (losBeamAngle.Tx.Zen(:,i)./180).* (obj.antennaGainRes - 1) + 1; % between 1 to Resolution
                        losBeamAngle.Tx.AziIndex(:,i) = ((losBeamAngle.Tx.Azi(:,i) + 180)./360).* (obj.antennaGainRes - 1) + 1; % between 1 to Resolution
                        
                        losBeamAngle.Rx.Ele(:,i) = -1 .* losBeamAngle.Tx.Ele(:,i);
                        losBeamAngle.Rx.Zen(:,i) = abs(90 - losBeamAngle.Rx.Ele(:,i));
                        losBeamAngle.Rx.Azi(:,i) = atan2(RxTx.vec.xyz(:,2,i),RxTx.vec.xyz(:,1,i)) * (180/pi);
                        losBeamAngle.Rx.Azi(find(isnan(losBeamAngle.Rx.Azi(:,i)) == 1),i) = 0; % Azimuth angle is unlikely to be nan due to use of atan2
                        losBeamAngle.Rx.ZenIndex(:,i) = (losBeamAngle.Rx.Zen(:,i)./180).* (obj.antennaGainRes - 1) + 1; % between 1 to Resolution
                        losBeamAngle.Rx.AziIndex(:,i) = ((losBeamAngle.Rx.Azi(:,i) + 180)./360).* (obj.antennaGainRes - 1) + 1; % between 1 to Resolution
                    end
                    
                    Tx2RxWalljd = zeros(size(wall.xyz1,1),1);
                    Tx2RxWalljxyz = zeros(size(wall.xyz1,1),3);
                    Tx2RxVec = zeros(size([obj.xt(ind_source) obj.yt(ind_source) obj.zt],1),3);
                    Tx2RxIntersectingWalls = zeros(size(Tx2RxWalljd));
                    Rx.LosRssi = zeros(size(Rx.xyz,1),1);
                    
                    for k = 1:size(Rx.xyz,1)
                        for i = 1:size([obj.xt(ind_source) obj.yt(ind_source) obj.zt],1)
                            Tx2RxVec(i,:) = Rx.xyz(k,:) - [obj.xt(ind_source) obj.yt(ind_source) obj.zt]; % i index is not really needed-
                            Tx2RxIntersectingWalls = zeros(size(Tx2RxWalljd));
                            incidentAngle = zeros(size(Tx2RxWalljd));
                            tempFresnelCoeff = ones(size(Tx2RxWalljd));
                            for j = 1:size(wall.xyz1,1)
                                % find intersection with each wall and validate it
                                Tx2RxWalljd(j) = dot(wall.xyz1(j,:) - [obj.xt(ind_source) obj.yt(ind_source) obj.zt], wall.normal.xyz(j,:),2) ./ dot(Tx2RxVec(i,:) , wall.normal.xyz(j,:),2); % Scalar value of the line between TX & Rx
                                if (Tx2RxWalljd(j)<1 && Tx2RxWalljd(j)>0)
                                    Tx2RxWalljxyz(j,:) = Tx2RxWalljd(j) .* Tx2RxVec(i,:) + [obj.xt(ind_source) obj.yt(ind_source) obj.zt]; % Intersection point with wall j
                                    if (prod(wall.minMax.x(j,:) - Tx2RxWalljxyz(j,1),2) < eps) && (prod(wall.minMax.y(j,:)...
                                            - Tx2RxWalljxyz(j,2),2) < eps) && (prod(wall.minMax.z(j,:) - Tx2RxWalljxyz(j,3),2) < eps)
                                        % At this point the intersection is definite
                                        Tx2RxIntersectingWalls(j) = 1;
                                        % Angle between the beam and intersecting wall
                                        incidentAngle(j) = acosd(abs(dot(wall.normal.xyz(j,:),Tx2RxVec(i,:),2)./...
                                            (sqrt(sum(wall.normal.xyz(j,:).^2,2)) .* sqrt(sum(Tx2RxVec(i,:).^2,2)))));
                                        if j < (size(wall.xyz1,1) - size(ceillFloor.xyz1,1) + 1) % checks if it's wall or ceiling
                                            if obj.polarizationSwap == 1 %% Polarization check for walls
                                                tempFresnelCoeff (j) = wall.TE.transFac(j,(round(incidentAngle(j))+1));  % Invoking Temporary Transmission Coefficients for walls
                                            elseif obj.polarizationSwap == 0
                                                tempFresnelCoeff (j) = wall.TM.transFac(j,(round(incidentAngle(j))+1));  % Invoking Temporary Transmission Coefficients for walls
                                            end
                                        else
                                            % Polarizaton check for ceiling
                                            if obj.polarizationSwap == 1
                                                tempFresnelCoeff (j) = wall.TM.transFac(j,(round(incidentAngle(j))+1));  % Invoking Temporary Transmission Coefficients for Ceiling
                                            elseif obj.polarizationSwap == 0
                                                tempFresnelCoeff (j) = wall.TE.transFac(j,(round(incidentAngle(j))+1));  % Invoking Temporary Transmission Coefficients for Ceiling
                                            end
                                        end
                                    end
                                end
                            end
                            %Antenna Pattern and ETC can be considered in this one later.
                            
                            % 1 - This only measures the loss not the received power so it only works for
                            %     one TX at the moment.
                            % 2-  log(FresnelCoeff) need to be subtracted as it procudes negative number, log(0:1) < 0
                            Rx.LosRssi(k) = Rx.LosRssi(k) +  10.^((obj.ptx(ind_source)- (obj.FPSLRefLoss + 10*obj.path_loss_exp.* log10(4*pi*((RxTx.dist(k,1,i) >= obj.refDistance)...
                                .* RxTx.dist(k,1,i) + (RxTx.dist(k,1,i)< obj.refDistance)*obj.refDistance) .* obj.f ./ obj.c) - 10.*log10(prod(tempFresnelCoeff))) ...
                                + (TxAntennaGainAE(round(losBeamAngle.Tx.AziIndex(k,i)),round(losBeamAngle.Tx.ZenIndex(k,i)))) ...
                                + (RxAntennaGainAE(round(losBeamAngle.Rx.AziIndex(k,i)),round(losBeamAngle.Rx.ZenIndex(k,i)))))/10)...
                                .* complex(cos(2*pi*obj.f*RxTx.dist(k,1,i)./obj.c + pi) , sin(2*pi*obj.f*RxTx.dist(k,1,i)./obj.c + pi));
                            %               Rx.LosRssi(k) = 10.^(Rx.LosRssi(k)./10) .* complex(cos(2*pi*obj.f*RxTx.dist(k,1,i)./c + pi) , sin(2*pi*obj.f*RxTx.dist(k,1,i)./c + pi)); % converting to linear and multiply by complex carrier
                        end
                    end
                    
                else
                    Rx.LosRssi = zeros(size(Rx.xyz,1),1);
                end
                
                
                %% Calculating Multipath & Reflection Component
                
                if obj.reflectionFlag == 1
                    for k = 1:size(Rx.xyz,1)
                        Rx.distFirstRefl=zeros(1,size([obj.xt(ind_source) obj.yt(ind_source) obj.zt],1));
                        for i = 1:size([obj.xt(ind_source) obj.yt(ind_source) obj.zt],1)
                            Rx.reflecjRssi = zeros(1,size(Tx.wallReflec.xyz,1));
                            firstOrderRef=zeros(1,size(Tx.wallReflec.xyz,1)); % computes reflections from all walls
                            for j = 1:size(Tx.wallReflec.xyz) % this is same as size(wall.xyz1,1)
                                TxRef2Rx.vec.xyz = Rx.xyz(k,:) - Tx.wallReflec.xyz(j,:,i);
                                TxRef2RxRefwallIntd = dot(wall.xyz1(j,:) - Tx.wallReflec.xyz(j,:,i), wall.normal.xyz(j,:),2) ./ dot(TxRef2Rx.vec.xyz, wall.normal.xyz(j,:),2);
                                
                                %                 tempReflecCoeff = 1; % This is in case there is no reflection
                                
                                if (TxRef2RxRefwallIntd < 1 && TxRef2RxRefwallIntd > 0) % d checks the reflection possibility from TxImage j
                                    reflectPointj = TxRef2RxRefwallIntd .* TxRef2Rx.vec.xyz + Tx.wallReflec.xyz(j,:,i);
                                    % xyz check the reflection possibility from TxImage j
                                    if(prod(wall.minMax.x(j,:) - reflectPointj(1,1),2) < eps) && (prod(wall.minMax.y(j,:)...
                                            - reflectPointj(1,2),2) < eps) && (prod(wall.minMax.z(j,:) - reflectPointj(1,3),2) < eps)
                                        % At this point there is a path for reflection, now
                                        % 1- Find the reflection coefficient
                                        % 2- Count the walls between reflection paths
                                        
                                        % 1- Finding Reflection Coefficient
                                        tempReflecAngle = acosd(abs(dot(TxRef2Rx.vec.xyz,wall.normal.xyz(j,:),2) ./ ((sqrt(sum(TxRef2Rx.vec.xyz.^2,2))...
                                            .* sqrt(sum(wall.normal.xyz(j,:).^2))))));
                                        if  j < (size(wall.xyz1,1) - size(ceillFloor.xyz1,1) + 1) % if panel is a wall
                                            if obj.polarizationSwap == 1
                                                tempReflecCoeff = wall.TE.refFac(j,(round(tempReflecAngle)+1));
                                            else
                                                tempReflecCoeff = wall.TM.refFac(j,(round(tempReflecAngle)+1));
                                            end
                                        else % if panel is either ceiling or floor
                                            if obj.polarizationSwap == 1
                                                tempReflecCoeff = wall.TM.refFac(j,(round(tempReflecAngle)+1));
                                            else
                                                tempReflecCoeff = wall.TE.refFac(j,(round(tempReflecAngle)+1));
                                            end
                                        end
                                        
                                        
                                        % 2- now that there is reflection, lets find the walls between the reflection paths
                                        Tx2ReflectPointj = reflectPointj - [obj.xt(ind_source) obj.yt(ind_source) obj.zt];
                                        reflectPointj2Rx = Rx.xyz(k,:) - reflectPointj;
                                        Tx2ReflectPointjDist = sqrt(sum(Tx2ReflectPointj.^2,2));
                                        reflectPointj2RxDist = sqrt(sum(reflectPointj2Rx.^2,2));
                                        
                                        
                                        Tx2ReflectPointIntersecWall  = zeros(size(wall.xyz1,1),1);
                                        reflectPointj2RxIntersecWall = zeros(size(wall.xyz1,1),1);
                                        
                                        % There is reflection so find the antenna gain and
                                        % beam departure angle, departure angle for
                                        % reflections, is the angle between reflection
                                        % point of the reflecting wall and the TX image.
                                        
                                        depBeamAngle.Ele = asind(Tx2ReflectPointj(1,3) ./ sqrt(sum(Tx2ReflectPointj.^2,2)));  % Elevation angle (between beam and Z plane not it's normal) -90<ele<90 degrees
                                        depBeamAngle.Azi = atan2(Tx2ReflectPointj(1,2),Tx2ReflectPointj(1,1)) * (180/pi); % Azimuth angle (between x and beam) -180<azi<180
                                        
                                        if isnan(depBeamAngle.Ele)
                                            depBeamAngle.Ele = 0; % if nan turns the beam angle to 0
                                        end
                                        
                                        depBeamAngle.Zen = abs(90-depBeamAngle.Ele);  % Zenith angle is calculated and used to find the antenna gain
                                        
                                        if isnan(depBeamAngle.Azi)
                                            depBeamAngle.Azi = 0; % Azimuth angle is unlikely to be nan due to use of atan2
                                        end
                                        
                                        depBeamAngle.ZenIndex = (depBeamAngle.Zen /180).* (obj.antennaGainRes - 1) + 1; % between 1 to obj.antennaGainRes
                                        depBeamAngle.AziIndex = ((depBeamAngle.Azi + 180)/360).* (obj.antennaGainRes - 1) + 1; % between 1 to Resolution
                                        
                                        
                                        % Also Calculating the Angle of Arrival
                                        arrBeamAngle.Ele = asin(-reflectPointj2Rx(1,3) ./ sqrt(sum(reflectPointj2Rx.^2,2)));
                                        arrBeamAngle.Azi = atan2(-reflectPointj2Rx(1,2),-reflectPointj2Rx(1,1)) * (180/pi);
                                        
                                        if isnan(arrBeamAngle.Ele)
                                            arrBeamAngle.Ele = 0; % if nan turns the beam angle to 0
                                        end
                                        
                                        arrBeamAngle.Zen = abs(90-arrBeamAngle.Ele);  % Zenith angle is calculated and used to find the antenna gain
                                        
                                        if isnan(arrBeamAngle.Azi)
                                            arrBeamAngle.Azi = 0; % Azimuth angle is unlikely to be nan due to use of atan2
                                        end
                                        
                                        arrBeamAngle.ZenIndex = (arrBeamAngle.Zen /180).* (obj.antennaGainRes - 1) + 1; % between 1 to obj.antennaGainRes
                                        arrBeamAngle.AziIndex = ((arrBeamAngle.Azi + 180)/360).* (obj.antennaGainRes - 1) + 1; % between 1 to Resolution
                                        
                                        
                                        % a- Find Walls Between TX to Reflection point
                                        tempTx2ReflpointTransCoeff = ones(size(wall.xyz1,1),1);
                                        tempReflpoint2RxTransCoeff = ones(size(wall.xyz1,1),1);
                                        for s = 1:size(wall.xyz1,1)
                                            % Finding Scalar value of intersection lines
                                            Tx2ReflectPointjWallsd = dot(wall.xyz1(s,:) - [obj.xt(ind_source) obj.yt(ind_source) obj.zt],wall.normal.xyz(s,:),2)./dot(Tx2ReflectPointj,wall.normal.xyz(s,:),2);
                                            reflectPointj2Rxd = dot(wall.xyz1(s,:) - reflectPointj,wall.normal.xyz(s,:),2) ./ dot(reflectPointj2Rx,wall.normal.xyz(s,:),2);
                                            % Checking for finite plane intersection
                                            if (Tx2ReflectPointjWallsd < 1 && Tx2ReflectPointjWallsd > 0 && abs(Tx2ReflectPointjWallsd - 1) > eps)
                                                Tx2ReflectPointjWallsxyz = Tx2ReflectPointjWallsd .* Tx2ReflectPointj + [obj.xt(ind_source) obj.yt(ind_source) obj.zt];
                                                if(prod(wall.minMax.x(s,:) - Tx2ReflectPointjWallsxyz(1,1),2) < eps) && (prod(wall.minMax.y(s,:)...
                                                        - Tx2ReflectPointjWallsxyz(1,2),2) < eps) && (prod(wall.minMax.z(s,:) - Tx2ReflectPointjWallsxyz(1,3),2) < eps)
                                                    % At this point wall s in between
                                                    Tx2ReflectPointIntersecWall(s) = 1;
                                                    intercepWallsIncAngle.Tx2ReflPoint(s) = acosd(abs(dot(wall.normal.xyz(s,:),Tx2ReflectPointj,2)./...
                                                        (sqrt(sum(wall.normal.xyz(s,:).^2,2)) .* sqrt(sum(Tx2ReflectPointj.^2,2)))));
                                                    if  s < (size(wall.xyz1,1) - size(ceillFloor.xyz1,1) + 1) % Transmission Coeffs for the intercepting walls between TX and refl point
                                                        if obj.polarizationSwap == 1
                                                            tempTx2ReflpointTransCoeff(s) = wall.TE.transFac(round(intercepWallsIncAngle.Tx2ReflPoint(s))+1);
                                                        else
                                                            tempTx2ReflpointTransCoeff(s) = wall.TM.transFac(round(intercepWallsIncAngle.Tx2ReflPoint(s))+1);
                                                        end
                                                    else % if panel is a wall
                                                        if obj.polarizationSwap == 1
                                                            tempTx2ReflpointTransCoeff(s) = wall.TM.transFac(round(intercepWallsIncAngle.Tx2ReflPoint(s))+1);
                                                        else
                                                            tempTx2ReflpointTransCoeff(s) = wall.TE.transFac(round(intercepWallsIncAngle.Tx2ReflPoint(s))+1);
                                                        end
                                                    end
                                                end
                                            end
                                            
                                            % b- Find finite walls Between Reflection point and RX. if below has complicated rule, as the reflectPointj2Rxd tend to
                                            % be smaller than matlab's epsilon and sometimes a little bit smaller than epsilon but bigger than epsm
                                            if (reflectPointj2Rxd < 1 && reflectPointj2Rxd > 0 && abs(reflectPointj2Rxd - 1) > eps && not(reflectPointj2Rxd < 1e-6))
                                                reflectPointj2Rxxyz = reflectPointj2Rxd .* reflectPointj2Rx + reflectPointj;
                                                if(prod(wall.minMax.x(s,:) - reflectPointj2Rxxyz(1,1),2) < eps) && (prod(wall.minMax.y(s,:)...
                                                        - reflectPointj2Rxxyz(1,2),2) < eps) && (prod(wall.minMax.z(s,:) - reflectPointj2Rxxyz(1,3),2) < eps)
                                                    % At this point wall s in between (reflection point to Rx)
                                                    reflectPointj2RxIntersecWall(s,1) = 1;
                                                    intercepWallsIncAngle.ReflPoint2Rx(s) = acosd(abs(dot(wall.normal.xyz(s,:),reflectPointj2Rx,2)./...
                                                        (sqrt(sum(wall.normal.xyz(s,:).^2,2)) .* sqrt(sum(reflectPointj2Rx.^2,2)))));
                                                    if  s < (size(wall.xyz1,1) - size(ceillFloor.xyz1,1) + 1) % Transmission Coeffs for the intercepting walls between TX and refl point
                                                        if obj.polarizationSwap == 1
                                                            tempReflpoint2RxTransCoeff(s) = wall.TE.transFac(round(intercepWallsIncAngle.ReflPoint2Rx(s))+1);
                                                        else
                                                            tempReflpoint2RxTransCoeff(s) = wall.TM.transFac(round(intercepWallsIncAngle.ReflPoint2Rx(s))+1);
                                                        end
                                                    else % if panel is either ceiling or floor
                                                        if obj.polarizationSwap == 1
                                                            tempReflpoint2RxTransCoeff(s) = wall.TM.transFac(round(intercepWallsIncAngle.ReflPoint2Rx(s))+1);
                                                        else
                                                            tempReflpoint2RxTransCoeff(s) = wall.TE.transFac(round(intercepWallsIncAngle.ReflPoint2Rx(s))+1);
                                                        end
                                                    end
                                                    
                                                end
                                            end
                                        end % end of for S
                                        
                                        % found number of walls between reflection paths.
                                        % calculate the received signal at Rx from
                                        % reflection point on wall j
                                        
                                        
                                        %                       % This Considers distance from Tx image to Rx at once
                                        Rx.reflecjRssi(j) = 10^((obj.ptx(ind_source)- (obj.FPSLRefLoss + 10*obj.path_loss_exp*log10(4*pi...
                                            *(Tx2ReflectPointjDist+reflectPointj2RxDist) .* obj.f ./ obj.c)) + (10*log10(prod(tempTx2ReflpointTransCoeff)))...
                                            + (10*log10(tempReflecCoeff)) + (10*log10(prod(tempReflpoint2RxTransCoeff)))...
                                            + (TxAntennaGainAE(round(depBeamAngle.AziIndex),round(depBeamAngle.ZenIndex))) + ...
                                            (RxAntennaGainAE(round(arrBeamAngle.AziIndex),round(arrBeamAngle.ZenIndex))))/10) ...
                                            .* complex(cos(2*pi*obj.f* Rx2TxRefl.dist(k,1,j,i)./obj.c + pi) , sin(2*pi*obj.f*Rx2TxRefl.dist(k,1,j,i)./obj.c + pi));
                                        
                                        Rx.distFirstRefl(i)=Tx2ReflectPointjDist+reflectPointj2RxDist;
                                    end % end of if reflection exist
                                    % Ony where reflections took place
                                end % end of d check for reflection
                                firstOrderRef(j)=Rx.distFirstRefl(i);
                            end% end of for j
                            
                        end % end of for i
                        Rx.ReflecRssi(k,1) = sum(sum(Rx.reflecjRssi,1),2);
                    end
                else
                    Rx.reflecjRssi = zeros(1,size(Tx.wallReflec.xyz,1));
                    Rx.ReflecRssi = zeros(size(Rx.xyz,1),1);
                    firstOrderRef=zeros(1,size(Tx.wallReflec.xyz,1));
                end
                
                
                %% Caclulating Second Reflections (Only works for one Tx).
                if obj.secondReflectionFlag == 1
                    Rx.SecondRefRSSI = zeros(size(Rx.xyz,1),1);
                    for i = 1:size(Rx.xyz,1)
                        Rx.SeconReflWallJRSSI = zeros(1,size(wall.xyz1,1));
                        secondOrderReflec=zeros(1,size(wall.xyz1,1));
                        for j = 1:size(Tx.secondReflecWallj.xyz,3)
                            %Initializing Parameters for Tx to First Reflection ponit
                            Tx2FirstReflPintIntersecWalls = zeros(size(wall.xyz1,1),1); % logs the wall between Tx and first reflection point on wall J
                            Tx2FirstReflPintTransCoeff = ones(size(wall.xyz1,1),1); % logs the Trans coeff of the wall between the Tx and the first reflection point
                            %Initializing Parameters for First to Second Reflection ponit
                            first2SecondReflPintIntersecWalls = zeros(size(wall.xyz1,1),1);
                            first2SecondReflPintTransCoeff = ones(size(wall.xyz1,1),1);
                            %Initializing Parameters for SECOND to Rx path
                            second2RxIntersecWalls = zeros(size(wall.xyz1,1),1);
                            second2RxReflPintTransCoeff = ones(size(wall.xyz1,1),1);
                            Rx.SecondReflWallKRSSI = zeros(size(wall.xyz1,1),1);
                            for k = 1:size(Tx.secondReflecWallj.xyz,1)
                                if (sum(Tx.secondReflecWallj.xyz(k,:,j) ~= [obj.xt(ind_source) obj.yt(ind_source) obj.zt]) ~= 0) % checks if the Tx.secondReflecWallj lies on the Tx
                                    TxSecondRef2Rx.vec.xyz = Rx.xyz(i,:) - Tx.secondReflecWallj.xyz(k,:,j);
                                    TxSecondRef2Rx.dist = sqrt(sum(TxSecondRef2Rx.vec.xyz.^2,2));
                                    
                                    % Find intersection of Tx.secondReflecWallj with wall k
                                    TxSecndRef2wallKIntd = dot(wall.xyz1(k,:) - Tx.secondReflecWallj.xyz(k,:,j), wall.normal.xyz(k,:),2) ...
                                        ./ dot(TxSecondRef2Rx.vec.xyz, wall.normal.xyz(k,:),2);
                                    if (TxSecndRef2wallKIntd < 1 && TxSecndRef2wallKIntd > 0) % check if there is intersection between the TxSecondreflection and the wall K
                                        secondReflectPointK = TxSecndRef2wallKIntd .* TxSecondRef2Rx.vec.xyz + Tx.secondReflecWallj.xyz(k,:,j);
                                        % now check if secondReflectPointK actually lies on the finite plane(wall) K
                                        if (prod(wall.minMax.x(k,:) - secondReflectPointK(1,1),2) < eps) && (prod(wall.minMax.y(k,:)...
                                                - secondReflectPointK(1,2),2) < eps) && (prod(wall.minMax.z(k,:) - secondReflectPointK(1,3),2) < eps)
                                            % At this point there is a path for second reflection, now chekc if there is a valid
                                            % path for first reflection (LOS between first reflection and secondReflectPointK intersects with wall j
                                            TxRefj2secondReflectPointK.vec.xyz = secondReflectPointK - Tx.wallReflec.xyz(j,:);
                                            TxRefj2SecondReflectPointKWalljIntd = dot(wall.xyz1(j,:) - Tx.wallReflec.xyz(j,:), wall.normal.xyz(j,:),2)...
                                                ./ dot(TxRefj2secondReflectPointK.vec.xyz, wall.normal.xyz(j,:),2); % Find the intersection of wall j with TxRef2secondReflectPointK
                                            
                                            if (TxRefj2SecondReflectPointKWalljIntd < 1 && TxRefj2SecondReflectPointKWalljIntd > 0) % check if there is intersection
                                                % now that there is intersection
                                                firstReflecPointj = TxRefj2SecondReflectPointKWalljIntd .* TxRefj2secondReflectPointK.vec.xyz + Tx.wallReflec.xyz(j,:);
                                                
                                                if (prod(wall.minMax.x(j,:) - firstReflecPointj(1,1),2) < eps) && (prod(wall.minMax.y(j,:)...
                                                        - firstReflecPointj(1,2),2) < eps) && (prod(wall.minMax.z(j,:) - firstReflecPointj(1,3),2) < eps) % check if the intersection lies on a finite plane
                                                    % At this point there is a path for second & first reflections,
                                                    % 1- Find the reflection coefficient for both first and second reflections
                                                    % 2- Count the walls between reflection paths
                                                    
                                                    % 1- Finding Reflection Coefficient
                                                    tempSecondReflecAngle = acosd(abs(dot(TxSecondRef2Rx.vec.xyz,wall.normal.xyz(k,:),2)...
                                                        ./ ((sqrt(sum(TxSecondRef2Rx.vec.xyz.^2,2)) .* sqrt(sum(wall.normal.xyz(k,:).^2))))));
                                                    tempFirstReflecAngle  = acosd(abs(dot(TxRefj2secondReflectPointK.vec.xyz,wall.normal.xyz(j,:),2)...
                                                        ./ ((sqrt(sum(TxRefj2secondReflectPointK.vec.xyz.^2,2)) .* sqrt(sum(wall.normal.xyz(j,:).^2))))));
                                                    
                                                    % Second Reflection factors baised on wall K for second reflections
                                                    if  k < (size(wall.xyz1,1) - size(ceillFloor.xyz1,1) + 1) % if panel is a wall
                                                        if obj.polarizationSwap == 1
                                                            tempSecondReflecCoeff = wall.TE.refFac(k,(round(tempSecondReflecAngle)+1));
                                                        else
                                                            tempSecondReflecCoeff = wall.TM.refFac(k,(round(tempSecondReflecAngle)+1));
                                                        end
                                                    else % if panel is either ceiling or floor
                                                        if obj.polarizationSwap == 1
                                                            tempSecondReflecCoeff = wall.TM.refFac(k,(round(tempSecondReflecAngle)+1));
                                                        else
                                                            tempSecondReflecCoeff = wall.TE.refFac(k,(round(tempSecondReflecAngle)+1));
                                                        end
                                                    end
                                                    
                                                    
                                                    % First Reflection factors baised on wall j for second reflections
                                                    if  j < (size(wall.xyz1,1) - size(ceillFloor.xyz1,1) + 1) % if panel is a wall
                                                        if obj.polarizationSwap == 1
                                                            tempFirstReflecCoeff = wall.TE.refFac(j,(round(tempFirstReflecAngle)+1));
                                                        else
                                                            tempFirstReflecCoeff = wall.TM.refFac(j,(round(tempFirstReflecAngle)+1));
                                                        end
                                                    else % if panel is either ceiling or floor
                                                        if obj.polarizationSwap == 1
                                                            tempFirstReflecCoeff = wall.TM.refFac(j,(round(tempFirstReflecAngle)+1));
                                                        else
                                                            tempFirstReflecCoeff = wall.TE.refFac(j,(round(tempFirstReflecAngle)+1));
                                                        end
                                                    end
                                                    
                                                    % The path of second reflection breaks down into 3 parts. Tx to First Reflection point. First to second
                                                    % Reflection point. Second reflection point to RX. For each path, walls in between need ot be checked
                                                    
                                                    Tx2FirstReflPoint.vec.xyz = firstReflecPointj - [obj.xt(ind_source) obj.yt(ind_source) obj.zt];
                                                    distDoubleReflTx2FirstPoint=sqrt(sum(Tx2FirstReflPoint.vec.xyz.^2,2));
                                                    firstReflPoint2SecondReflPoint.vec.xyz  = secondReflectPointK - firstReflecPointj;
                                                    distDoubleReflFirstPoint2SecondPoint=sqrt(sum(firstReflPoint2SecondReflPoint.vec.xyz.^2,2));
                                                    secondReflPoint2Rx.vec.xyz = Rx.xyz(i,:) - secondReflectPointK;
                                                    distDoubleReflSecondPoint2Rx=sqrt(sum(secondReflPoint2Rx.vec.xyz.^2,2));
                                                    %                                                     secondOrderReflec=[secondOrderReflec distDoubleReflTx2FirstPoint+distDoubleReflFirstPoint2SecondPoint+distDoubleReflSecondPoint2Rx]
                                                    % checking number of walls between TX and firstReflectionPoint
                                                    for l = 1:size(wall.xyz1,1) % checking number of walls between TX and firstReflectionPoint
                                                        if (l ~= j) % unecessary, only for safety measures as intersection D for same wall is zero.obj
                                                            wallLIntdTx2FirstReflPoint = dot(wall.xyz1(l,:) - [obj.xt(ind_source) obj.yt(ind_source) obj.zt], wall.normal.xyz(l,:),2)...
                                                                ./ dot(Tx2FirstReflPoint.vec.xyz, wall.normal.xyz(l,:),2);
                                                            if (wallLIntdTx2FirstReflPoint < 1 && wallLIntdTx2FirstReflPoint > 0)
                                                                wallLintTx2Tx2FirstReflPoint = wallLIntdTx2FirstReflPoint .* Tx2FirstReflPoint.vec.xyz ...
                                                                    + [obj.xt(ind_source) obj.yt(ind_source) obj.zt];
                                                                if (prod(wall.minMax.x(l,:) - wallLintTx2Tx2FirstReflPoint(1,1),2)...
                                                                        < eps) && (prod(wall.minMax.y(l,:) - wallLintTx2Tx2FirstReflPoint(1,2),2)...
                                                                        < eps) && (prod(wall.minMax.z(l,:) - wallLintTx2Tx2FirstReflPoint(1,3),2) < eps)
                                                                    % now wall L is in between. Find the angle of incidence and transmission coefficient
                                                                    Tx2FirstReflPintIntersecWalls(l,1) = 1; % logging the wall in between
                                                                    tempWallInterceptingAngle(l) = acosd(abs(dot(wall.normal.xyz(l,:),Tx2FirstReflPoint.vec.xyz,2)...
                                                                        ./(sqrt(sum(wall.normal.xyz(l,:).^2,2)) .* sqrt(sum(Tx2FirstReflPoint.vec.xyz.^2,2))))); % finds the angle between
                                                                    if  l < (size(wall.xyz1,1) - size(ceillFloor.xyz1,1) + 1) % Transmission Coeffs for the intercepting walls between TX and refl point
                                                                        if obj.polarizationSwap == 1
                                                                            Tx2FirstReflPintTransCoeff(l) = wall.TE.transFac(round(tempWallInterceptingAngle(l))+1);
                                                                        else
                                                                            Tx2FirstReflPintTransCoeff(l) = wall.TM.transFac(round(tempWallInterceptingAngle(l))+1);
                                                                        end
                                                                    else % if panel is either ceiling or floor
                                                                        if obj.polarizationSwap == 1
                                                                            Tx2FirstReflPintTransCoeff(l) = wall.TM.transFac(round(tempWallInterceptingAngle(l))+1);
                                                                        else
                                                                            Tx2FirstReflPintTransCoeff(l) = wall.TE.transFac(round(tempWallInterceptingAngle(l))+1);
                                                                        end
                                                                    end
                                                                end % if (prod(wall.minMax.x(j,:) - wallLintTx2Tx2FirstReflPoint(1,1),2) < eps) && (prod(wall.minMax.y(j,:) - wallLintTx2Tx2FirstReflPoint(1,2),2) < eps) && (prod(wall.minMax.z(j,:) - wallLintTx2Tx2FirstReflPoint(1,3),2) < eps)
                                                            end % if (wallLintdTx2FirstReflPoint < 1 && wallLintdTx2FirstReflPoint > 0).
                                                        end % if (l ~= j)
                                                    end % for l = 1:size(wall.xyz1,1)
                                                    
                                                    % checking number of walls between the FIRSTReflectionPoint and SECONDreflection point
                                                    for l = 1:size(wall.xyz1,1) % First Reflection Pint to Second
                                                        wallLIntdFirstRefl2SecondReflPoint = dot(wall.xyz1(l,:) - firstReflecPointj, wall.normal.xyz(l,:),2)...
                                                            ./ dot(firstReflPoint2SecondReflPoint.vec.xyz, wall.normal.xyz(l,:),2);
                                                        if (wallLIntdFirstRefl2SecondReflPoint < 1 && wallLIntdFirstRefl2SecondReflPoint > 0)
                                                            wallLIntFirstRefl2SecondReflPoint = wallLIntdFirstRefl2SecondReflPoint .* firstReflPoint2SecondReflPoint.vec.xyz ...
                                                                + firstReflecPointj;
                                                            if (prod(wall.minMax.x(l,:) - wallLIntFirstRefl2SecondReflPoint(1,1),2) < eps) && (prod(wall.minMax.y(l,:)...
                                                                    - wallLIntFirstRefl2SecondReflPoint(1,2),2) < eps) && (prod(wall.minMax.z(l,:)...
                                                                    - wallLIntFirstRefl2SecondReflPoint(1,3),2) < eps)
                                                                % now wall L is in between. Find the angle of incidence and transmission coefficient
                                                                
                                                                first2SecondReflPintIntersecWalls(l,1) = 1; % logging the wall in between
                                                                tempWallInterceptingAngle(l) = acosd(abs(dot(wall.normal.xyz(l,:),firstReflPoint2SecondReflPoint.vec.xyz,2)...
                                                                    ./(sqrt(sum(wall.normal.xyz(l,:).^2,2)) .* sqrt(sum(firstReflPoint2SecondReflPoint.vec.xyz.^2,2))))); % finds the angle between
                                                                if  l < (size(wall.xyz1,1) - size(ceillFloor.xyz1,1) + 1) % Transmission Coeffs for the intercepting walls between TX and refl point
                                                                    if obj.polarizationSwap == 1
                                                                        first2SecondReflPintTransCoeff(l) = wall.TE.transFac(round(tempWallInterceptingAngle(l))+1);
                                                                    else
                                                                        first2SecondReflPintTransCoeff(l) = wall.TM.transFac(round(tempWallInterceptingAngle(l))+1);
                                                                    end
                                                                else % if panel is either ceiling or floor
                                                                    if obj.polarizationSwap == 1
                                                                        first2SecondReflPintTransCoeff(l) = wall.TM.transFac(round(tempWallInterceptingAngle(l))+1);
                                                                    else
                                                                        first2SecondReflPintTransCoeff(l) = wall.TE.transFac(round(tempWallInterceptingAngle(l))+1);
                                                                    end
                                                                end
                                                                
                                                            end % if (prod(wall.minMax.x(j,:) - wallLintTx2Tx2FirstReflPoint(1,1),2) < eps) && (prod(wall.minMax.y(j,:) - wallLintTx2Tx2FirstReflPoint(1,2),2) < eps) && (prod(wall.minMax.z(j,:) - wallLintTx2Tx2FirstReflPoint(1,3),2) < eps)
                                                        end % if (wallLintdTx2FirstReflPoint < 1 && wallLintdTx2FirstReflPoint > 0).
                                                        
                                                    end % for l = 1:size(wall.xyz1,1)
                                                    
                                                    
                                                    
                                                    % checking number of walls between the SECOND reflection point to the RX.
                                                    for l = 1:size(wall.xyz1,1) % SECOND to RX
                                                        if (l ~= k) % To Avoid the Same wall that the second Reflection is bouncing of! Although this is for safety really
                                                            wallLIntdSecondRef2Rx = dot(wall.xyz1(l,:) - Rx.xyz(i,:), wall.normal.xyz(l,:),2) ...
                                                                ./ dot(secondReflPoint2Rx.vec.xyz, wall.normal.xyz(l,:),2);
                                                            if (wallLIntdSecondRef2Rx < 1 && wallLIntdSecondRef2Rx > 0)
                                                                wallLIntSecondRef2Rx = wallLIntdSecondRef2Rx .* secondReflPoint2Rx.vec.xyz + secondReflectPointK;
                                                                if (prod(wall.minMax.x(l,:) - wallLIntSecondRef2Rx(1,1),2) < eps) && (prod(wall.minMax.y(l,:) ...
                                                                        - wallLIntSecondRef2Rx(1,2),2) < eps) && (prod(wall.minMax.z(l,:) - wallLIntSecondRef2Rx(1,3),2) < eps)
                                                                    % now wall L is in between. Find the angle of incidence and transmission coefficient
                                                                    
                                                                    second2RxIntersecWalls(l,1) = 1; % logging the wall in between
                                                                    tempWallInterceptingAngle(l) = acosd(abs(dot(wall.normal.xyz(l,:),secondReflPoint2Rx.vec.xyz,2)...
                                                                        ./ (sqrt(sum(wall.normal.xyz(l,:).^2,2)) .* sqrt(sum(secondReflPoint2Rx.vec.xyz.^2,2))))); % finds the angle between
                                                                    if  l < (size(wall.xyz1,1) - size(ceillFloor.xyz1,1) + 1) % Transmission Coeffs for the intercepting walls between TX and refl point
                                                                        if obj.polarizationSwap == 1
                                                                            second2RxReflPintTransCoeff(l) = wall.TE.transFac(round(tempWallInterceptingAngle(l))+1);
                                                                        else
                                                                            second2RxReflPintTransCoeff(l) = wall.TM.transFac(round(tempWallInterceptingAngle(l))+1);
                                                                        end
                                                                    else % if panel is either ceiling or floor
                                                                        if obj.polarizationSwap == 1
                                                                            second2RxReflPintTransCoeff(l) = wall.TM.transFac(round(tempWallInterceptingAngle(l))+1);
                                                                        else
                                                                            second2RxReflPintTransCoeff(l) = wall.TE.transFac(round(tempWallInterceptingAngle(l))+1);
                                                                        end
                                                                    end
                                                                    
                                                                end % if (prod(wall.minMax.x(j,:) - wallLintTx2Tx2FirstReflPoint(1,1),2) < eps) && (prod(wall.minMax.y(j,:) - wallLintTx2Tx2FirstReflPoint(1,2),2) < eps) && (prod(wall.minMax.z(j,:) - wallLintTx2Tx2FirstReflPoint(1,3),2) < eps)
                                                            end % if (wallLintdTx2FirstReflPoint < 1 && wallLintdTx2FirstReflPoint > 0).
                                                        end % if (l ~= k)
                                                    end % for l = 1:size(wall.xyz1,1)
                                                    
                                                    % There is reflection so find the antenna gain and
                                                    % beam departure angle, departure angle for
                                                    % reflections, is the angle between reflection
                                                    % point of the reflecting wall and the TX image.
                                                    
                                                    depBeamAngle.Ele = asind(Tx2FirstReflPoint.vec.xyz(1,3) ./ sqrt(sum(Tx2FirstReflPoint.vec.xyz.^2,2)));  % Elevation angle (between beam and Z plane not it's normal) -90<ele<90 degrees
                                                    depBeamAngle.Azi = atan2(Tx2FirstReflPoint.vec.xyz(1,2),Tx2FirstReflPoint.vec.xyz(1,1)) * (180/pi); % Azimuth angle (between x and beam) -180<azi<180
                                                    
                                                    if isnan(depBeamAngle.Ele)
                                                        depBeamAngle.Ele = 0; % if nan turns the beam angle to 0
                                                    end
                                                    
                                                    depBeamAngle.Zen = abs(90-depBeamAngle.Ele);  % Zenith angle is calculated and used to find the antenna gain
                                                    
                                                    if isnan(depBeamAngle.Azi)
                                                        depBeamAngle.Azi = 0; % Azimuth angle is unlikely to be nan due to use of atan2
                                                    end
                                                    
                                                    depBeamAngle.ZenIndex = (depBeamAngle.Zen /180).* (obj.antennaGainRes - 1) + 1; % between 1 to obj.antennaGainRes
                                                    depBeamAngle.AziIndex = ((depBeamAngle.Azi + 180)/360).* (obj.antennaGainRes - 1) + 1; % between 1 to Resolution
                                                    
                                                    
                                                    % Also Calculating the Angle of Arrival
                                                    arrBeamAngle.Ele = asin(-secondReflPoint2Rx.vec.xyz(1,3) ./ sqrt(sum(secondReflPoint2Rx.vec.xyz.^2,2)));
                                                    arrBeamAngle.Azi = atan2(-secondReflPoint2Rx.vec.xyz(1,2),-secondReflPoint2Rx.vec.xyz(1,1)) * (180/pi);
                                                    
                                                    if isnan(arrBeamAngle.Ele)
                                                        arrBeamAngle.Ele = 0; % if nan turns the beam angle to 0
                                                    end
                                                    
                                                    arrBeamAngle.Zen = abs(90-arrBeamAngle.Ele);  % Zenith angle is calculated and used to find the antenna gain
                                                    
                                                    if isnan(arrBeamAngle.Azi)
                                                        arrBeamAngle.Azi = 0; % Azimuth angle is unlikely to be nan due to use of atan2
                                                    end
                                                    
                                                    arrBeamAngle.ZenIndex = (arrBeamAngle.Zen /180).* (obj.antennaGainRes - 1) + 1; % between 1 to obj.antennaGainRes
                                                    arrBeamAngle.AziIndex = ((arrBeamAngle.Azi + 180)/360).* (obj.antennaGainRes - 1) + 1; % between 1 to Resolution
                                                    
                                                    
                                                    % Calculating the Second Reflection of Wall 1:K for First Reflection being of wall J
                                                    
                                                    Rx.SecondReflWallKRSSI(k) = 10^((obj.ptx(ind_source) - (obj.FPSLRefLoss + 10*obj.path_loss_exp*log10(4*pi*(TxSecondRef2Rx.dist) ...
                                                        .* obj.f ./ obj.c)) + (10*log10(prod(Tx2FirstReflPintTransCoeff)))...
                                                        + (10*log10(tempFirstReflecCoeff)) + (10*log10(prod(first2SecondReflPintTransCoeff))) ...
                                                        + (10*log10(tempSecondReflecCoeff)) + (10*log10(prod(second2RxReflPintTransCoeff)))...
                                                        + (TxAntennaGainAE(round(depBeamAngle.AziIndex),round(depBeamAngle.ZenIndex))) + ...
                                                        (RxAntennaGainAE(round(arrBeamAngle.AziIndex),round(arrBeamAngle.ZenIndex))))/10) ...
                                                        .* complex(cos(2*pi*obj.f*TxSecondRef2Rx.dist./obj.c) , sin(2*pi*obj.f*TxSecondRef2Rx.dist./obj.c));
                                                    
                                                end % if (prod(wall.minMax.x(j,:) - secondReflectPointK(1,1),2) < eps) && (prod(wall.minMax.y(j,:) - secondReflectPointK(1,2),2) < eps) && (prod(wall.minMax.z(j,:) - secondReflectPointK(1,3),2) < eps) % check if the intersection lies on a finite plane
                                            end % if (TxRefj2SecondReflectPointKWalljIntd < 1 && TxRefj2SecondReflectPointKWalljIntd > 0)
                                        end %(prod(wall.minMax.x(j,:) - secondReflectPointK(1,1),2) < eps) && (prod(wall.minMax.y(j,:) - secondReflectPointK(1,2),2) < eps) && (prod(wall.minMax.z(j,:) - secondReflectPointK(1,3),2) < eps)
                                    end %(TxSecndRef2wallKIntd < 1 && TxSecndRef2wallKIntd > 0)
                                    % check the validity of the first projection being on the wall
                                end % if Tx.secondReflecWallj ~= Tx.xyz
                            end % for k = size(Tx.secondReflecWallj.xyz,1)
                            secondOrderReflec(j)=TxSecondRef2Rx.dist;
                            Rx.SeconReflWallJRSSI(j) = sum(Rx.SecondReflWallKRSSI);
                        end % for j = 1:size(Tx.secondReflecWallj.xyz,3)
                        Rx.SecondRefRSSI(i,1) = sum(Rx.SeconReflWallJRSSI);
                    end % for i = 1:size(Rx.xyz,1)
                else % if secondReflectionFlag == 1
                    Rx.SeconReflWallJRSSI = zeros(1,size(wall.xyz1,1));
                    Rx.SecondRefRSSI = zeros(size(Rx.xyz,1),1);
                    secondOrderReflec=zeros(1,size(wall.xyz1,1));
                end % if secondReflection == 1
                
                %% Reflection & Line Of Signt Propagation Map
                totalReceivedPowerAllSources(ind_source) = 10*log10(abs(Rx.LosRssi+ (obj.reflectExaggerationFac * Rx.ReflecRssi)+...
                    (obj.reflectExaggerationFac * Rx.SecondRefRSSI)));
                
                
                distanceNray=[RxTx.dist  firstOrderRef  secondOrderReflec]; %  Rx.distFirstRefl  TxSecondRef2Rx.dist
                delays=distanceNray/obj.c;
                
                powerDelayProfile=10*log10([abs(Rx.LosRssi) abs(Rx.reflecjRssi) abs(Rx.SeconReflWallJRSSI)]);   % abs(Rx.ReflecRssi) abs(Rx.SecondRefRSSI)
                
                %
                alpha_coeff = sqrt(db2pow(powerDelayProfile));
                [ h_D, ~] = obj.digitalImpulseResponse(alpha_coeff, delays);
                
                %%
                H_D(ind_source, 1:length(h_D))=h_D; % Get the impulse responses from all sources
                if obj.estimateLocFreeSpace==1
                    D_direct(ind_source)=RxTx.dist; % tries the free space propagation in TDOA
                end
                %
                
            end
            
            
            if obj.estimateLocFreeSpace==1
                for ind_sourceHd=2:n_sources
                    averaDistAllSources(ind_sourceHd-1)=(D_direct(ind_sourceHd)-D_direct(1));% tries the free space propagation in TDO
                end
            end
            
            estimatedTdoaRange=averaDistAllSources;
            totalReceivedPower=10*log10(sum(db2pow(totalReceivedPowerAllSources)));
        end
        
    end
    
end

