% Coastline Determination
clc, clear
close all

% Simulate Orbit Parameters/Get Position/lat/long of Coastline Points from MATLAB
load scdatawithrot.mat % converts orbital elements to position and velocity vectors to simulate orbit of spacecraft
load ECEF.mat % converted lat/long data from MATLAB to ECEF coordinates
load coast % provided MATLAB lat/long points to map coastlines (latlong2ECEF.m, checked with ECEF2latlong.m)

% Define constants
rsc = Pos_rot;
Vel = Vel_rot;

%%%%%%%%% Full resolution GSHHG database %%%%%%%%%%%%%%%%%%%
% world = gshhs('gshhs_l.b');
% figure
% worldmap world
% geoshow([world.Lat], [world.Lon])
% figure
% worldmap(map, refvec)
% geoshow(gca,map,refvec,'DisplayType','texturemap');
% demcmap(map)
% S = shaperead('landareas','UseGeoCoords',true);
% geoshow([S.Lat], [S.Lon],'Color','black');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Map image to Ellipsoid
[A,R] = geotiffread('TrueMarble32km.tif');
info = geotiffinfo('TrueMarble32km.tif');
a = 1000*6378.137; % equatorial radius (km)
b = 1000*6356.752; % polar radius (km)
[x,y,z] = ellipsoid(0,0,0,a,a,b);
globe = surf(x,y,-z,'EdgeColor','none');
axis equal, axis off
h1 = globe;
h2 = axesm('globe','Geoid',1000*6378);
% h3 = gridm('GLineStyle','-','Gcolor',[.8 .8 .8]);
h4 = plotm(lat,long);
set(globe,'CData',A,'FaceColor','texturemap');

% Specify Desired FOV and calculate footprint size (GIFOV)
FOV = (pi/180)*30; % deg
GIFOV = 2*height*tan(FOV/2); % km
FOV_new = adjustFOV(GIFOV); % adjust to Matlab's scale [deg]

n = 5*60;
for k = 1:n:length(t) %9000 = Hawaii
    if k+1 > length(t)
        break
    end
    
    if Vel(k,3) < 0
        camup([0 0 -1])
        camva(FOV_new); % deg
        campos([rsc(k,1) rsc(k,2) rsc(k,3)]);
    else
        camup([0 0 1])
        camva(FOV_new); % deg
        campos([rsc(k,1) rsc(k,2) rsc(k,3)]);
    end
    drawnow
%                  waitforbuttonpress
    
    % Find FOV footprint limits in terms of lat/long
    delta_lat = GIFOV/110.54e3; % convert to m [1deg = 110.54 km]
    delta_long = GIFOV/(111.32e3*cos(lat_calc(k)*(pi/180))); % convert to m [1deg = 111.32*cos(lat)]
    
    latbound(k,:) = [lat_calc(k)+delta_lat;
        lat_calc(k)+delta_lat;
        lat_calc(k)-delta_lat;
        lat_calc(k)-delta_lat];
    
    longbound(k,:) = [long_calc(k)+GIFOV/(111.32e3*cos(latbound(k,1)*(pi/180)));
        long_calc(k)-GIFOV/(111.32e3*cos(latbound(k,2)*(pi/180)));
        long_calc(k)-GIFOV/(111.32e3*cos(latbound(k,3)*(pi/180)));
        long_calc(k)+GIFOV/(111.32e3*cos(latbound(k,4)*(pi/180)))];
    
    coordfind = find(long < longbound(k,1) & long < longbound(k,4) & long > longbound(k,3) ...
        & long > longbound(k,2) & lat > latbound(k,3) & lat < latbound(k,2));
    
    if isempty(coordfind) == 1
        fprintf('No Coastline in image at k = %d \n',k)
    else
        coast = [long(coordfind) lat(coordfind)];
        for i = 1:20
            random = randi(size(coast,1));
            coast_plot(i,:) = coast(random,:);
        end
    end
end

%% Plots
figure
hold on
plot(long,lat,'k','LineWidth',1)
plot(long_calc(1,1:end), lat_calc(1,1:end),'.','MarkerSize',0.75,'Color',[0.6, 0.6, 0.6])
for k = 1:n:length(longbound)
    hold on
    plot([longbound(k,:) longbound(k,1)],[latbound(k,:) latbound(k,1)],'b','LineWidth',1.25,'Color',[0.8, 0.8, 0.8]);
end
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

set(findall(gcf,'type','text'),'FontSize',16)

%% Get Position Vector of "known" Landmark
% datacursormode on
% dcmObj = datacursormode(gcf);
% set(dcmObj,'SnapToDataVertex','off','Enable','on')
% 
% disp('Click Point and Press Enter:')
% pause
% location = getCursorInfo(dcmObj);
% ox1 = location.Position(1);
% oy1 = location.Position(2);
% oz1 = location.Position(3);
% o = [ox1; oy1; oz1];
% close, clc

ox1 = Pos_rot(9000,1);
oy1 = Pos_rot(9000,2);
oz1 = Pos_rot(9000,3);
o = [ox1; oy1; oz1];

[latitude1,longitude1,height1] = ECEF2latlong(ox1,oy1,oz1);
% latitude1
% longitude1
ECEF_surf = latlong2ECEF(latitude1,longitude1,0); % position on surface of sphere
ECEF = latlong2ECEF(latitude1,longitude1,height1);
diff = (o'-ECEF_surf) % m

% disp('Click Second Point and Press Enter:')
% pause
% point = getCursorInfo(dcmObj);
% Px2 = point.Position(1);
% Py2 = point.Position(2);
% Pz2 = point.Position(3);
% o2 = [Px2 Py2 Pz2];
% 
% [latitude2,longitude2,height2] = ECEF2latlong(Px2,Py2,Pz2);
% ECEF_surf2 = latlong2ECEF(latitude2,longitude2,0); % position on surface of sphere
% 
% d = latlong2distance(latitude1,longitude1,latitude2,longitude2)
% d1 = sqrt(abs(Px2-ox1)^2+abs(Py2-oy1)^2+abs(Pz2-oz1)^2)

% Map back to ellipsoid
figure
[x,y,z] = ellipsoid(0,0,0,a,a,b);
globe = surf(x,y,-z,'EdgeColor','none');
axis equal, axis off
h1 = globe;
h2 = axesm('globe','Geoid',1000*6378);
% h3 = gridm('GLineStyle','-','Gcolor',[.8 .8 .8]);
h4 = plotm(lat,long);
set(globe,'CData',A,'FaceColor','texturemap');
hold on
plot3(ECEF_surf(1), ECEF_surf(2), ECEF_surf(3),'r*');
% plot3(ECEF_surf2(1), ECEF_surf2(2), ECEF_surf2(3),'r*');
view([ECEF(1), ECEF(2), ECEF(3)])
camup([0 0 1])
% s = ECEF-Pos_rot(n,:);
