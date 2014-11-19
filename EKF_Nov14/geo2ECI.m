function [ x,y,z ] = geo2ECI(lat,long,height)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% JDtime = 2456964.98264;
Re = 6378000;
% deltaT = 1; % 1 sec
% ang_vel = 7.292e-5; % rad/s
% angle_sid = ang_vel*deltaT; % rad
theta = long;% + angle_sid;
r = (height+Re)*cos(lat);

x = r*cos(theta);
y = r*sin(theta);
z = (height+Re)*sin(lat);

end

