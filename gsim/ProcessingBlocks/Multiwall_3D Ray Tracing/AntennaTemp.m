%%
% PhiMax = 2*pi; % Max Azimuth Angle
% ThetaMax = -pi/2;  % Max Elevation Angle

function [antennaGainAE] = AntennaTemp (antennaGainRes,demoMode)
Length = 1;
% antennaGainRes = 20;

phi= linspace(0,2*pi,antennaGainRes);
theta= linspace(0,pi,antennaGainRes);
[THETA,PHI]=meshgrid(theta,phi);
R=abs((cos(pi.*Length.*cos(THETA))- cos(pi.*Length))./(sin(THETA+eps)));
x=R.*sin(THETA).*cos(PHI);
y=R.*sin(THETA).*sin(PHI);
z=R.*cos(THETA);

% if demoMode == 1
% 
%      figure
%      mesh(x,y,z);
 antennaGainAE = R;
end

% title(['Antenna Pattern, Resolution = ',num2str(antennaGainRes)]);