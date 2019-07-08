% Assigns Permittivity to each single wall
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

figure
fill3(wallX, wallY, wallZ,wallC)



if demoMode == 1
    for i = 1:numWallsOnly
        figure('Name',['Panel Viewer: Panel #',num2str(i),'/',num2str(size(wall.xyz1,1))])
        wallC = ones(size(wallX));
        wallC(:,i) = 40;
        fill3(wallX,wallY,wallZ,wallC,'faceColor','Flat')
        title(['Panel #',num2str(i),'/',num2str(size(wall.xyz1,1))])
        if strcmpi (inputMode,'singleIntrupted')
            
%             currentPerm  = inputdlg(['Enter permittivity for panel #',num2str(i),':'],...
%                 'Permittivity Assigner');
            currentPerm = input(['Enter permittivity for panel #',num2str(i),': ']);
            wall.relativePerm(i) = currentPerm;
        end
    end
end


wallX = [wall.xyz1(:,1)';wall.xyz2(:,1)';wall.xyz3(:,1)';wall.xyz4(:,1)'];
wallY = [wall.xyz1(:,2)';wall.xyz2(:,2)';wall.xyz3(:,2)';wall.xyz4(:,2)'];
wallZ = [wall.xyz1(:,3)';wall.xyz2(:,3)';wall.xyz3(:,3)';wall.xyz4(:,3)'];
wallC = ones(size(wallX)); 


if demoMode == 1
    for i = numWallsOnly+1:size(wall.xyz1,1)
        figure('Name',['Panel Viewer: Panel #',num2str(i),'/',num2str(size(wall.xyz1,1))])
        wallC = ones(size(wallX));
        wallC(:,i) = 40;
        fill3(wallX,wallY,wallZ,wallC,'faceColor','Flat')
        title(['Panel #',num2str(i),'/',num2str(size(wall.xyz1,1))])
        if strcmpi (inputMode,'singleIntrupted')
            
%             currentPerm  = inputdlg(['Enter permittivity for panel #',num2str(i),':'],...
%                 'Permittivity Assigner');
            currentPerm = input(['Enter permittivity for panel #',num2str(i),': ']);
            wall.relativePerm(i) = currentPerm;
        end
    end
end



        