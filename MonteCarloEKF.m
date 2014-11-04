% Coastline Determination with Extended Kalman Filter
clc, clear
% close all

% Simulate Orbit Parameters/Get Position/lat/long of Coastline Points from MATLAB
load scdata.mat % from scdata.m
load coord.mat
% load ECEF.mat % converted lat/long data from MATLAB to ECEF coordinates
% load coast % Built-in MATLAB lat/long points to map coastlines (checked with latlong2ECEF.m and ECEF2latlong.m)
% coord = [lat, long]; % put coordinates in matrix

%%%%%%%%%% Full resolution GSHHG database %%%%%%%%%%%%%%%%%%%
% world = gshhs('gshhs_l.b');
% figure
% worldmap world
% geoshow([world.Lat], [world.Lon])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Map image to Ellipsoid
% A2 = imread('map.png'); % import image to map to sphere
% [A,R] = geotiffread('TrueMarble32km.tif'); % import satellite image
% info = geotiffinfo('TrueMarble32km.tif');

% % Model Earth as Ellipsoid
% a = 1000*6378.137; % equatorial radius (m)
% b = 1000*6356.752; % polar radius (m)
% [x,y,z] = ellipsoid(0,0,0,a,a,b);
% globe = surf(x,y,-z,'EdgeColor','none');
% axis equal, axis off
% h1 = globe;
% h2 = axesm('globe','Geoid',1000*6378);
% h3 = gridm('GLineStyle','-','Gcolor',[.8 .8 .8]);
% h4 = plotm(lat,long);
% set(globe,'CData',A,'FaceColor','texturemap');

% Preliminary Data and Such

% Specify Desired FOV and calculate footprint size (GIFOV)
FOV = 30*(pi/180); % rad
GIFOV = 2*height*tan(FOV/2); % km
FOV_new = adjustFOV(GIFOV); % adjust to Matlab's scale [deg]

% Stuff that needs initialized for EKF
rsc = Position; % rename known spacecraft position vector
Vel = Velocity; % convert known velocity vectors from m/s
xtrue = [rsc'; Vel'];

sigma_pos = 100; % 100 m
sigma_vel = .001; % 1 m/s
sigma_theta = 0.01*(pi/180); % 0.01 deg (converted to rad)

T_IC = eye(3,3); % rotation matrix from camera to intertial frame, assuming camera is nadir

close % close ellipsoid image

% Initial Conditions
rsc = xtrue(1:3,1) + 1*sigma_pos*randn(3,1); % use initial position with added noise
rscdot = xtrue(4:6,1) + 1*sigma_vel*randn(3,1); % use initial velocity with added noise
x = [rsc; rscdot]; % create initial state matrix
P = [sigma_pos^2*eye(3), zeros(3); zeros(3), sigma_vel^2*eye(3)]; % create initial P matrix
z = [x; reshape(P,36,1)]; % compile state and covariance into column vector

% Find estimated position of satellite based on camera image
n = 30;
for mc = 1 % for monte carlo analysis
    mc
    for k_old = 1:n:length(t)-1
        
        if k_old+1 > length(t)
            break
        end
        
        if k_old == 1
            k = 1;
        else
            k = ((k_old-1)/n)+1;
        end
        
%         % Find FOV footprint limits in terms of lat/long coordinates
%         delta_lat = GIFOV/110.54e3; % convert to m [1deg = 110.54 km]
%         delta_long = GIFOV/(111.32e3*cos(lat_calc(k_old)*(pi/180))); % convert to m [1deg = 111.32*cos(lat)]
%         top = [lat_calc(k_old)+delta_lat, long_calc(k_old)+delta_long]; % convert to m [1deg = 111.32*cos(lat)]
%         left = [lat_calc(k_old)+delta_lat, long_calc(k_old)-delta_long];
%         bottom = [lat_calc(k_old)-delta_lat, long_calc(k_old)-delta_long];
%         right = [lat_calc(k_old)-delta_lat, long_calc(k_old)+delta_long];
%         box = [top; left; bottom; right];
%         figure, plot(box(1,1),box(1,2),'*',box(2,1),box(2,2),'*',box(3,1),box(3,2),'*',box(4,1),box(4,2),'*')
% dist = acos(sin(lat(k)*(pi/180))*sin(lat(k+1)*(pi/180))+cos(lat(k)*(pi/180))*cos(lat(k+1)*(pi/180))*cos((long(k)-long(k+1))*(pi/180)));

        
%         % Check for coastline points within footprint limits
%         coordfind = find(lat < lat_max & lat > lat_min & long < long_max & long > long_min);
%         if isempty(coordfind) == 0
%             latlong_est = [lat(coordfind(1)) long(coordfind(1))];
%             pos_est(k,:) = latlong2ECEF(latlong_est(1),latlong_est(2),0); % convert lat/long point to ECEF point on surface
%         else
%             pos_est(k,:) = xtrue(1:3,k-1);
%                         fprintf('No Coastline at k = %f. \n',k)
%             continue
%         end
        pos_est(k,:) = latlong2ECEF(llh(k_old,1),llh(k_old,2),llh(k_old,3))';
         
        o = [pos_est(k,1); pos_est(k,3); pos_est(k,3)]; % vector of landmark
        dv = sigma_theta*randn(3,1); % generate random noise for position vector
        dT = v_to_T(dv); % rotation matrix based on random noise in position vector of camera
        h = dT*T_IC*(o-x(1:3))/(sqrt((o-x(1:3))'*(o-x(1:3)))); % map LOS vector to spacecraft-object position vector
        y = h+dv; % find estimate of position measurement
        
        % Integrate State and Covariance from tk-1 to tk
        [t1, z1] = ode45(@integrate,t(k):1:t(k+1),z,options); % integrate z to find xdot and Pdot
        z = z1(end,:);
        
        % Obtain Estimated State Values
        xhatminus(:,k) = [z(1:6)]';
        Pminus = reshape(z(:,7:42),6,6); % reshape into 6x6 matrix
        
        %%% Update Measurement %%%
        s = o-xhatminus(1:3,k); % find LOS vector
        s_norm = norm(s); % normalize LOS vector
        eI = s/s_norm; % find unit vector in direction of inertial frame
        
        % Rotate from body to inertial frame
        ex = -xhatminus(1:3,k)/norm(xhatminus(1:3,k));
        ez = cross(ex,xhatminus(4:6,k)/norm(xhatminus(4:6,k)));
        ey = cross(ez,ex);
        ez = ez/norm(ez);
        ey = ey/norm(ey);
        T_BI = [ex'; ey'; ez'];
        
        H = T_BI*T_IC*(1/s_norm)*[eI*eI'-eye(3) zeros(3)]; % measurement sensitivity matrix
        
        R = sigma_theta^2*(eye(3)-eI*eI'); % measurement covariance matrix
        
        K = (Pminus*H')/((H*Pminus*H'+R+(0.5*trace(R))*(eI*eI'))); % find optimal Kalman gain
        
        xhatplus(:,k) = xhatminus(:,k) + K*(y-h); % state update
        Pplus = (eye(6)-K*H)*Pminus*(eye(6)-K*H)' + K*R*K'; % covariance update
        
        PlotP(k,:) = reshape(Pplus,36,1); % reshape P into vector
        errX(:,k) = [xhatplus(1:3,k)-xtrue(1:3,k_old+1); (xhatplus(4:6,k)-xtrue(4:6,k_old+1))]; % Plot error
        
        x = xhatplus(:,k); % update state measurement
        P = Pplus; % update covariance measurement
        z = [x; reshape(P,36,1)]; % compile state and covariance into column vector
        
    end
    
    % Plot Position Error
    figure(1)
    subplot(3,1,1)
    hold on
    plot(sqrt(PlotP(:,1)),'k')
    plot(-sqrt(PlotP(:,1)),'k')
    plot(errX(1,:),'b')
    title('Error in X Position (m)')
    xlabel('Time (min)')
    ylabel('Magnitude (m)')
    grid on
    
    subplot(3,1,2)
    hold on
    plot(sqrt(PlotP(:,8)),'k')
    plot(-sqrt(PlotP(:,8)),'k')
    plot(errX(2,:),'b')
    title('Error in Y Position (m)')
    xlabel('Time (min)')
    ylabel('Magnitude (m)')
    grid on
    
    subplot(3,1,3)
    hold on
    plot(sqrt(PlotP(:,15)),'k')
    plot(-sqrt(PlotP(:,15)),'k')
    plot(errX(3,:),'b')
    title('Error in Z Position (m)')
    xlabel('Time (min)')
    ylabel('Magnitude (m)')
    grid on
    
    % Plot Velocity Error
    figure(2)
    subplot(3,1,1)
    hold on
    plot(sqrt(PlotP(:,22)),'k')
    plot(-sqrt(PlotP(:,22)),'k')
    plot(errX(4,:),'b')
    title('Error in X Velocity (m/s)')
    xlabel('Time (min)')
    ylabel('Magnitude (mm/s)')
    grid on
    
    subplot(3,1,2)
    hold on
    plot(sqrt(PlotP(:,29)),'k')
    plot(-sqrt(PlotP(:,29)),'k')
    plot(errX(5,:),'b')
    title('Error in Y Velocity (m/s)')
    xlabel('Time (min)')
    ylabel('Magnitude (mm/s)')
    grid on
    
    subplot(3,1,3)
    hold on
    plot(sqrt(PlotP(:,36)),'k')
    plot(-sqrt(PlotP(:,36)),'k')
    plot(errX(6,:),'b')
    title('Error in Z Velocity (m/s)')
    xlabel('Time (min)')
    ylabel('Magnitude (mm/s)')
    grid on
    
end

% set(findall(gcf,'type','text'),'FontSize',18,'fontWeight','bold')
