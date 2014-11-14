function [ rcb,rorbit,mu ] = cbinfo(cb)
%[mu]=cbinfor(planet) This function contains planetary data
%
%   Input central body name as a string with first letter capitalized
%
%   OUTPUTS:     rcb   - central body radius (km)
%               rorbit - central body's orbital radius about sun (km)
%                        (NOTE: if cb==Moon, rorbit is the semimajor axis of the
%                        Moon's orbit about Earth
%                 mu   - gravitational parameter (km^3/sec^2)
%
%Created by Andrew Liounis on 1/23/13 9:36am

AU=149.6*10^6;  %creating AU distance unit (km)
if strcmp(cb,'Mercury')==1
    rcb=2439.7; %km
    rorbit=.287*AU;  %km
    mu=2.094*10^4;  %km^3/sec^2
    
elseif strcmp(cb,'Venus')==1
    rcb=6051.8; %km
    rorbit=.723*AU;  %km
    mu=3.249*10^5;  %km^3/sec^2
    
elseif strcmp(cb,'Earth')==1
    rcb=6378.14;  %km
    rorbit=AU;  %km
    mu=3.986004418e5;  %km^3/sec^2
    
elseif strcmp(cb,'Mars')==1
    rcb=3397.2;  %km
    rorbit=1.524*AU;  %km
    mu=4.269*10^4;  %km^3/sec^2
    
elseif strcmp(cb,'Jupiter')==1
    rcb=71492;  %km
    rorbit=5.204*AU;  %km
    mu=1.267*10^8;  %km^3/sec^2
    
elseif strcmp(cb,'Saturn')==1
    rcb=60268;  %km
    rorbit=9.582*AU;  %km
    mu=3.7967*10^7;  %km^3/sec^2
    
elseif strcmp(cb,'Uranus')==1
    rcb=25559;  %km
    rorbit=19.2*AU;  %km
    mu=5.7918*10^6;  %km^3/sec^2
    
elseif strcmp(cb,'Neptune')==1
    rcb=24764;  %km
    rorbit=30.05*AU;  %km
    mu=6.806*10^6;  %km^3/sec^2
    
elseif strcmp(cb,'Pluto')==1
    rcb=1195;  %km
    rorbit=39.24*AU;  %km
    mu=798.04;  %km^3/sec^2
    
elseif strcmp(cb,'Sun')==1
    rcb=6.95508*10^5;  %km
    rorbit=0;
    mu=1.327124*10^11;  %km^3/sec^2 
    
elseif strcmp(cb,'Moon')
    rcb=1737;  %km .4
    rorbit=384388.174;  %km  NOTE  this is the semimajor axis for the Moon's orbit about EARTH
    mu=4902.7949;  %km^3/sec^2
    
else
    error('No central body info available')
end

end

