function [R,V] = elementstoRV(alt,e,i,w,Omega,nu)

% This function calculates the position and velocity vectors of a spacecraft
% in inertial coordinates using the six orbital elements
%
%   INPUTS:  a   - semimajor axis (meters)
%            e   - eccentricity (unitless)
%            i   - inclination (radians)
%            w - right ascension of ascending node (radians)
%          Omega - argument of periapsis (radians)
%           nu   - true anomally (radians)
%
%   OUTPUTS: R   - position vector for spacecraft (meters)
%            V   - velocity vector for spacecraft (meters/second)

mu = (1000^3)*398600; % gravitational parameter of earth (m^3s^-2)

p = alt*(1-e^2); % semi-latus rectum
r = p/(1+e*cos(nu)); % length of radius

% R(1,1) = p*(cos(Omega)*cos(w+nu)-sin(Omega)*cos(i)*sin(w+nu));
% R(2,1) = p*(sin(Omega)*cos(w+nu)+cos(Omega)*cos(i)*sin(w+nu));
% R(3,1) = p*sin(i)*sin(w+nu);
% 
% V(1,1) = -sqrt(mu/p)*(cos(Omega)*(sin(w+nu)+e*sin(w))+sin(Omega)*cos(i)*(cos(w+nu)+e*cos(w)));
% V(2,1) = -sqrt(mu/p)*(sin(Omega)*(sin(w+nu)+e*sin(w))+cos(Omega)*cos(i)*(cos(w+nu)+e*cos(w)));
% V(3,1) = -sqrt(mu/p)*(sin(i)*(cos(w+nu)+e*cos(w)));

Rpqw = [r*cos(nu) p*sin(nu) 0]; % radius vector in PWQ coordinate frame
Vpqw = [-(sqrt(mu/p))*sin(nu) sqrt(mu/p)*(e+cos(nu)) 0]; % velocity vector in PWQ coordinate frame

rot = [cos(Omega)*cos(w)-sin(Omega)*sin(w)*cos(i) -cos(Omega)*sin(w)-sin(Omega)*cos(w)*cos(i) sin(Omega)*sin(i);
    sin(Omega)*cos(w)+cos(Omega)*sin(w)*cos(i) -sin(Omega)*sin(w)+cos(Omega)*cos(w)*cos(i) -cos(Omega)*sin(i);
    sin(w)*sin(i) cos(w)*sin(i) cos(i)]; % rotation matrix

R = rot*Rpqw'; % R to inertial coordinate frame
V = rot*Vpqw'; % V to inertial coordinate frame

end

