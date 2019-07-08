% Measurement Location Assigner 
% Created by: Salaheddin Hosseinzadeh
% Created on: 23 Aug 2016
% Last Revision:
% Notes:


[structFileName,structFileAddress] = uigetfile('.mat');
load([structFileAddress,structFileName]);

%%
inputMode = 'singleIntrupted';
demoMode  = 1;

numWallsOnly = size(wall.xyz1,1) - size(ceillFloor.xyz1,1);

wallX = [wall.xyz1(1:numWallsOnly,1)';wall.xyz2(1:numWallsOnly,1)';wall.xyz3(1:numWallsOnly,1)';wall.xyz4(1:numWallsOnly,1)'];
wallY = [wall.xyz1(1:numWallsOnly,2)';wall.xyz2(1:numWallsOnly,2)';wall.xyz3(1:numWallsOnly,2)';wall.xyz4(1:numWallsOnly,2)'];
wallZ = [wall.xyz1(1:numWallsOnly,3)';wall.xyz2(1:numWallsOnly,3)';wall.xyz3(1:numWallsOnly,3)';wall.xyz4(1:numWallsOnly,3)'];
wallC = ones(size(wallX)); 
% wallC
% wallC(:,1) = [1;.5;.5];

fh = figure
fill3(wallX, wallY, wallZ,wallC)
hold on

plot3(Rx.xyz(:,1),Rx.xyz(:,2),Rx.xyz(:,3),'LineStyle','none','Marker','.');

axis equal
grid on
view(2)


locationCounter = 0;
if ~exist('measurementLocationXyz','var')
    measurementLocationXyz = zeros(1,3);
end

while 1
    choice = input('Choose one of the Add, Delete, Edit, Close actions (A/D/E/C)? ','s');
    if strcmpi(choice,'A')
        locationCounter = locationCounter + 1;
        title(['Selection Location #',num2str(locationCounter)]);
        figure(fh)
        [x,y] = ginput(1);
        locationDist = 
        
        sum(abs(repmat([x,y],size(Rx.xyz,1),1) - Rx.xyz(:,1:2)),2);
        [~,RxIndex] = min(locationDist);
        plot3(Rx.xyz(RxIndex,1),Rx.xyz(RxIndex,2),Rx.xyz(RxIndex,3),'LineStyle','none','Marker','*','Color','Red');
        measurementLocationXyz (locationCounter,:) = Rx.xyz(RxIndex,:);
        disp(['Location #',num2str(locationCounter),'[',num2str(measurementLocationXyz(locationCounter,1))...
            ',',num2str(measurementLocationXyz(locationCounter,2)),']',' added Successfully.']);
        title(['Location #',num2str(locationCounter),'[',num2str(measurementLocationXyz(locationCounter,1))...
            ',',num2str(measurementLocationXyz(locationCounter,2)),']',' added Successfully.']);
    elseif strcmpi(choice,'D')
        % delete last selected point
        figure(fh)
        plot3(measurementLocationXyz(locationCounter,1),measurementLocationXyz(locationCounter,2),measurementLocationXyz(locationCounter,3),'LineStyle','none','Marker','.','Color','blue');
        measurementLocationXyz(locationCounter,:) = [];
        disp(['Measurement location #',num2str(locationCounter),' deleted!']);
        locationCounter = locationCounter - 1;
        measurementLocationXyz        
    elseif strcmpi(choice,'E')
        locEditIndex = input('Which location do you want to edit? ');
        figure(fh)
        [x,y] = ginput(1);
        locationDist = sum(abs(repmat([x,y],size(Rx.xyz,1),1) - Rx.xyz(:,1:2)),2);
        [~,RxIndex] = min(locationDist);
        plot3(Rx.xyz(RxIndex,1),Rx.xyz(RxIndex,2),Rx.xyz(RxIndex,3),'LineStyle','none','Marker','*','Color','Red');
        measurementLocationXyz (locEditIndex,:) = Rx.xyz(RxIndex,:);
        disp(['Location #',num2str(locEditIndex),'[',num2str(measurementLocationXyz(locEditIndex,1))...
            ',',num2str(measurementLocationXyz(locEditIndex,2)),']',' added Successfully.']);
        title(['Location #',num2str(locEditIndex),'[',num2str(measurementLocationXyz(locEditIndex,1))...
            ',',num2str(measurementLocationXyz(locEditIndex,2)),']',' added Successfully.']);
    elseif strcmpi(choice,'C')
        % Close the program
        measurementLocationXyz 
        break

    end
            
    
end
        