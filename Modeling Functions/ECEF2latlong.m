function [lat,long,height] = ECEF2latlong(x,y,z)
% [lat,long,height] = ECEF2latlong(x,y,z)
%   Convert ECEF position coordinates to lat/long coordinates

% n = 6905;
% x = Position(1,1);
% y = Position(1,2);
% z = Position(1,3);

a = 1000*6378;
b = 1000*6357;
e = sqrt(1-(b^2/a^2));

long = atan2(y,x);
p = sqrt(x^2+y^2);

for i = 1:4
    if i == 1
        lat_gc(i) = atan2(p,z); % rad
    else
        lat_gc(i) = atan2(p,z); % rad, geocentric
        Rn = a/sqrt(1-e^2*sin(lat_gc(i))^2); % m
        height = p*cos(lat_gc(i))+z*sin(lat_gc(i))-a/sqrt(1-e^2*sin(lat_gc(i))^2); % m
%         height = p/cos(lat_gc(i))-Rn;
        lat_gc(i+1) = atan((z/p)*(1-e^2*(Rn/(Rn+height))^-1)); % rad
    end
end

lat = 180/pi*lat_gc(end); % convert to deg
long = 180/pi*long; % convert to deg
% [lat long height]

% if lat > 90
%     lat_new = 90-lat;
% %     fprintf('Latitude is: %5.2f S \n',lat_new)
% else
%     lat_new = 90-lat;
% %     fprintf('Latitude is: %5.2f N \n',lat_new)
% end
% if long > 0
%     long_new = long;
% %     fprintf('Longitude is: %5.2f E \n',long_new)
% end
% if long < 0
%     long_new = long;
% %     fprintf('Longitude is: %5.2f W \n',long_new)
% end


end
