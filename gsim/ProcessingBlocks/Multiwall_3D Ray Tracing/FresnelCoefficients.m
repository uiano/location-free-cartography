% Fresnel Transmission and Reflection Coefficients
% Created by: Salaheddin Hosseinzadeh
% Created on: 03 /  May / 2016
% Last Revision:
% Notes:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% e1 & e2 are the epsilon (E) of the first and second material (1st is air)
% theta1 & theta2 are the incident and transmission angles in degree
% sin(theta1)/sin(theta2) = n2/n1, this is how transmission angle is calculated
% n1 & n2 are the reflection index of the materials and if permuability of
% the two medium are equal u1 = u2 then n1 = sqrt(e1)
% [TEReflFac,TETransFac,TMReflFac,TMTransFac] = FresnelCoefficients (e1,e2,theta1)
function [TEReflFac,TETransFac,TMReflFac,TMTransFac] = FresnelCoefficients (e1,e2,theta1,demoMode)

% clear all
e1 = 1; % for air or vacuum
% e2 = 2;


f = 1; % frequncy in GHz

c = 0; % set to zero to cancel conductivity and metal support
d = 1.60;


erC = c.* f.^d;  % imaginary part of relative perm, function of freq

eta = e2 - i.*erC; 

theta1 = 0:90;
theta = theta1;
theta1 = pi*theta1./180; % theta1 is in readian for ease of use in MATLAB



TEReflFac = abs((cos(theta1) - sqrt(eta - sin(theta1).^2))./(cos(theta1) + sqrt(eta - sin(theta1).^2)));
TMReflFac = abs((eta.*cos(theta1) - sqrt(eta - sin(theta1).^2))./(eta.*cos(theta1) + sqrt(eta - sin(theta1).^2)));

TETransFac = abs((2.*cos(theta1)) ./ (cos(theta1) + sqrt(eta - sin(theta1).^2)));
TMTransFac = abs((2.*sqrt(eta).*cos(theta1)) ./ (eta.*cos(theta1) + sqrt(eta - sin(theta1).^2)));


TEReflFac = TEReflFac.^2;
TMReflFac = TMReflFac.^2;

TETransFac = TETransFac.^2 .* sqrt(eta);
TMTransFac = TMTransFac.^2 .* sqrt(eta);


% Converting to power


if demoMode == 1
    figure('Name',['Fresnel Coeffs for Er = ',num2str(e2)])
    plot(theta,TEReflFac);
    hold on
    plot(theta,TETransFac,'r');

    % figure
    plot(theta,TMReflFac,'LineStyle','--')
    hold on
    plot(theta,TMTransFac,'LineStyle','--','Color','red')


    legend('TE Ref','TE Trans','TM Ref','TM Trans')
end




%{
n1 = sqrt(e1);
n2 = sqrt(e2);

% demo = 0;

 
% theta1 = 0:90;
theta = theta1;
theta1 = pi*theta1./180;

theta2 = asin(n1./n2.*sin(theta1));


TEReflFac = abs((n1.*cos(theta1) - n2.*sqrt(1 - (n1./n2 .* sin(theta1)).^2))./...
    (n1.*cos(theta1) + n2.*sqrt(1 - (n1./n2 .* sin(theta1)).^2))); % S Polarization

TETransFac = 1 - TEReflFac;


TMReflFac = abs((n1.*sqrt(1 - (n1/n2 .* sin(theta1)).^2) - n2.*cos(theta1))./...
    (n1.*(sqrt(1 - (n1/n2.*sin(theta1)).^2)) + n2.*cos(theta1))); % P Polarization

TMTransFac = (2.*n1.*cos(theta1))./(n1.*cos(theta2) + n2.*cos(theta1));

if demoMode == 1
    figure('Name',['Fresnel Coeffs for Er = ',num2str(e2)])
    plot(theta,TEReflFac);
    hold on
    plot(theta,TETransFac,'r');

    % figure
    plot(theta,TMReflFac,'LineStyle','--')
    hold on
    plot(theta,TMTransFac,'LineStyle','--','Color','red')


    legend('TE Ref','TE Trans','TM Ref','TM Trans')
end
%}

































