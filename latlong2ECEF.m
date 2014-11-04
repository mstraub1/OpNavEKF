function ECEF = latlong2ECEF(lat,long,h)
% [x,y,z] = latlong2ECEF(lat,long,h)
% This function calculates the ECEF position of the spacecraft given inputs
% of latitude, longitude, and height
%
%   INPUTS:  lat - latitude coordinates (degrees)
%            long - longitude coordinates (degrees)
%            h - height (meters)
%
%   OUTPUTS: ECEF - x, y, z coordinates of spacecraft in ECEF coordinate
%   frame

lat = (pi/180)*(lat);
long = (pi/180)*(long);

a = 1000*6378; % equatorial radius in m
f = 1/298.257224e3; % flattening parameter (a-b)/a where b is polar radius
% e = 0;
% N = a/sqrt(1-e^2*sind(lat));

C = 1/sqrt(cos(lat)^2+(1-f)^2*sin(lat)^2);
S = (1-f)^2*C;

x = (a*C+h)*cos(lat)*cos(long);
y = (a*C+h)*cos(lat)*sin(long);
z = (a*S+h)*sin(lat);

ECEF = [x y z];

end

