%%
clc, clear
close all

%%%%%% DETERMINE S/C POSITION AND VELOCITY BASED ON ORBITAL ELEMENTS %%%%%%
% Generate values for spacecraft's state (position, velocity, acceleration)
height = 1000*1000; % height above Earth (m) (assume circular orbit at 1000 km)
Rearth = 1000*6378; % radius of Earth (m)
alt = Rearth+height; % altitude of orbit (m)
e = 0; % eccentricity (rad) (circular orbit = 0)
i = 45*(pi/180); % inclination in rad (polar orbit = 90 deg)
w = 0; % right ascension of ascending node (rad)
Omega = 0; % argument of periapsis (rad)
nu = 0; % true anomaly (rad)
t = [0:24*60*60]'; % 12 hours of data - time vector for ode45 (sec)

% Propagate elements for length of time specified in t vector using ode45
[R,V] = elementstoRV(alt,e,i,w,Omega,nu);
x = [R; V]; % [m, m/s]
options = odeset('abstol',1e-8,'reltol',1e-8);
[t,z_true] = ode45(@earthgravity_m,t,x,options);
Pos_ECI = z_true(:,1:3); % ECI [m]
Vel_ECI = 1000*z_true(:,4:6); % ECI [mm/s]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find FOV footprint limits in terms of lat/long coordinates
% Specify Desired FOV and calculate footprint size (GIFOV)
FOV = 30*(pi/180); % rad
GIFOV = 2*(alt-Rearth)*tan(FOV/2); % km

% Define physical constants
rad2deg = 180/pi;
deg2rad = pi/180;
angvel = 7.292e-5; % angular velocity of Earth (rad/s)

% Load coastline data
load coord.mat;
lat = deg2rad*lat;
long = deg2rad*long;

% Define filter and MC constants
sigma_pos = 500; % 100 m
sigma_vel = 5; % mm/s
sigma_theta = 0.05*deg2rad; % 0.01 deg (converted to rad)

% Set integration options
dt = 30; % frequency of integration data
options = odeset('abstol',1e-8,'reltol',1e-8); % set tolerances for ode45

% Initialize variables
T_IC = eye(3,3); % rotation matrix from camera to inertial frame, assuming camera is nadir
rsc = Pos_ECI(:,1:3); % ECI [m]
Vel = Vel_ECI(:,1:3); % ECI [mm/s]
xtrue = [rsc'; Vel']; % ECI [m, mm/s]

% Create initial P matrix
P0 = [sigma_pos^2*eye(3), zeros(3); zeros(3), (sigma_vel)^2*eye(3)];

% Create true trajectory
t_meas_interval = 10*60;
t_start = 0;
t_end = t_start+t_meas_interval;
x = [xtrue(:,1); reshape(P0,36,1)];

Pos_ECEF = xtrue(1:3,1);
Pos_ECI_idx = xtrue(1:3,1);
idx = 2;

tplot = [];
xplot_true = [];
while t_end <= t(end)
    
    % Integrate over interval
    t_sweep = [t_start:dt:t_end];
    [tout,xout] = ode45(@integrate,t_sweep,x,options);
    
    % Append new interval to plot history
    tplot = [tplot;tout];
    xplot_true = [xplot_true;xout];
    
    % Reset start point
    x = xout(end,:)';
    t_start = t_end;
    t_end = t_start + t_meas_interval;
    
    % Compute ECEF location of new start point
    theta_ECI_to_ECEF = t_start*angvel;
    T_ECI_to_ECEF(:,:,idx) = [cos(theta_ECI_to_ECEF) -sin(theta_ECI_to_ECEF) 0;
        sin(theta_ECI_to_ECEF)  cos(theta_ECI_to_ECEF) 0;
        0 0 1];
    Pos_ECEF(:,idx) = T_ECI_to_ECEF(:,:,idx)*x(1:3);
    Pos_ECI_idx(:,idx) = x(1:3);
    idx = idx + 1;
    
end

NumUpdates = idx - 1;


for mc = 1:15 % for monte carlo analysis
    mc
    % Create perturbed initial state
    rsc    = xtrue(1:3,1) + sigma_pos*randn(3,1); % use initial position with added noise
    rscdot = xtrue(4:6,1) + sigma_vel*randn(3,1); % use initial velocity with added noise
    
    % Build state vector and covariance [m, mm/s]
    x = [rsc; rscdot]; % create initial state matrix
    P = [sigma_pos^2*eye(3), zeros(3); zeros(3), sigma_vel^2*eye(3)]; % create initial P matrix
    z = [x; reshape(P,36,1)]; % compile state and covariance into column vector
    
    t_start = 0;
    t_end = t_start+t_meas_interval;
    
    zplot_est = [];
    
    
    for idx = 2:NumUpdates
        
        % Integrate over interval
        t_sweep = [t_start:dt:t_end];
        [tout,zout] = ode45(@integrate,t_sweep,z,options);
        
        % Append new interval to plot history
        zplot_est = [zplot_est;zout];
        
        xhat_minus = zout(end,1:6)';
        P_minus = reshape( zout(end,7:end), 6, 6 );
        
        % Are measurements present?
        [lat_obs,long_obs,height_obs] =  ECEF2latlong(Pos_ECEF(1,idx),Pos_ECEF(2,idx),Pos_ECEF(3,idx));
        delta_lat = GIFOV/110.54e3; % convert to m [1deg = 110.54 km]
        delta_long = GIFOV/(111.32e3*cos(lat_obs*(pi/180))); % convert to m [1deg = 111.32*cos(lat)]
        top = [lat_obs+delta_lat, long_obs+delta_long]; % convert to m [1deg = 111.32*cos(lat)]
        left = [lat_obs+delta_lat, long_obs-delta_long];
        bottom = [lat_obs-delta_lat, long_obs-delta_long];
        right = [lat_obs-delta_lat, long_obs+delta_long];
        box = [top; left; bottom; right];
        coordfind = find(lat < max(box(:,1)) & lat > min(box(:,1)) & long < max(box(:,2)) & long > min(box(:,2)));
        
        % If coordfind is not empty, then there are measurements
        if isempty(coordfind) == 0
            
            % Find number of points to process (will either be 1 or 2)
            if length(coordfind) > 1
                latlong_est = [lat(coordfind(1:2)) long(coordfind(1:2))];
                number = 2;
            else
                latlong_est = [lat(coordfind(1)) long(coordfind(1))];
                number = 1;
            end
            
            %latlong_est = [lat_obs, long_obs];
            
            
            for num = 1:number
                
                % Compute ECEF location of point
                lat_num = latlong_est(num,1); % + 10*rad2deg*randn(1);
                long_num = latlong_est(num,2); % + 10*rad2deg*randn(1);
                o_ECEF = (height_obs + Rearth)*[ cos(lat_num)*cos(long_num);
                    cos(lat_num)*sin(long_num);
                    sin(lat_num)];
                
                o_ECI = T_ECI_to_ECEF(:,:,idx)'*o_ECEF   ;
                
                
                % Compute measurement
                s_true = o_ECI - Pos_ECI_idx(:,idx);
                eI_true = s_true/norm(s_true);
                y_true = T_IC*eI_true;
                dv = sigma_theta*randn(3,1); % generate random noise for position vector
                dT = v_to_T(dv); % rotation matrix based on random noise in position vector of camera
                y = dT*y_true;
                
                s = o_ECI-xhat_minus(1:3); % find LOS vector
                s_norm = norm(s); % normalize LOS vector
                eI = s/s_norm; % find unit vector in direction of inertial frame
                h = T_IC*eI;
                
                H = T_IC*(1/norm(s))*[eI*eI'-eye(3) zeros(3)]; % measurement sensitivity matrix
                
                R = sigma_theta^2*(eye(3)-eI*eI'); % measurement covariance matrix
                
                K = (P_minus*H')/(H*P_minus*H'+R+(0.5*trace(R))*(eI*eI')); % find optimal Kalman gain
                
                xhat_plus = xhat_minus + K*(y-h); % state update
                P_plus = (eye(6)-K*H)*P_minus*(eye(6)-K*H)' + K*R*K'; % covariance update
            end
            
        else
            
            fprintf('No Coastline at k = %d. \n',idx)
            P_plus = P_minus;
            xhat_plus = xhat_minus;
            
        end
        
        z = [xhat_plus; reshape(P_plus,36,1)];
        t_start = t_end;
        t_end = t_start + t_meas_interval;
        
        
        
        %idx = idx + 1;
        
    end
    
    
    Pos_err = xplot_true(:,1:3) - zplot_est(:,1:3);
    Vel_err = xplot_true(:,4:6) - zplot_est(:,4:6);
    
    for kk = 1:length(tplot)
        r_kk = xplot_true(kk,1:3)';
        v_kk = xplot_true(kk,4:6)';
        xxx = r_kk/norm(r_kk);
        yyy = cross(v_kk,xxx);
        yyy = yyy/norm(yyy);
        zzz = cross(xxx,yyy);
        zzz = zzz/norm(zzz);
        
        T_kk = [xxx yyy zzz]; %LVLH to inertial
        Pos_err_LVLH(kk,:) = ( T_kk'*Pos_err(kk,:)' )';
        Vel_err_LVLH(kk,:) = ( T_kk'*Vel_err(kk,:)' )';
        
        Pkk = reshape( zplot_est(kk,7:end), 6, 6);
        Tblock = [T_kk, zeros(3); zeros(3), T_kk];
        Pkk_LVLH = Tblock'*Pkk*Tblock;
        P_LVLH(kk,:) = reshape( Pkk_LVLH, 1, 36);
        Pmag(kk) = trace( Pkk_LVLH(1:3,1:3) );
        %if kk == 150
        %    Pp = Pkk(1:3,1:3)
        %   [evec,eval] = eig(Pp)
        %    T_kk
        %    pause
        %end
        
    end
    
    
    figure(1)
    subplot(3,1,1), hold on, grid on, plot(tplot./60, Pos_err(:,1), 'Color', 0.6*[1 1 1])
    subplot(3,1,2), hold on, grid on, plot(tplot./60, Pos_err(:,2), 'Color', 0.6*[1 1 1])
    subplot(3,1,3), hold on, grid on, plot(tplot./60, Pos_err(:,3), 'Color', 0.6*[1 1 1])
    
    figure(2)
    subplot(3,1,1), hold on, grid on, plot(tplot./60, Vel_err(:,1), 'Color', 0.6*[1 1 1])
    subplot(3,1,2), hold on, grid on, plot(tplot./60, Vel_err(:,2), 'Color', 0.6*[1 1 1])
    subplot(3,1,3), hold on, grid on, plot(tplot./60, Vel_err(:,3), 'Color', 0.6*[1 1 1])
    
    figure(3)
    subplot(3,1,1), hold on, grid on, plot(tplot./60, Pos_err_LVLH(:,1), 'Color', 0.6*[1 1 1])
    subplot(3,1,2), hold on, grid on, plot(tplot./60, Pos_err_LVLH(:,2), 'Color', 0.6*[1 1 1])
    subplot(3,1,3), hold on, grid on, plot(tplot./60, Pos_err_LVLH(:,3), 'Color', 0.6*[1 1 1])
    
    figure(4)
    subplot(3,1,1), hold on, grid on, plot(tplot./60, Vel_err_LVLH(:,1), 'Color', 0.6*[1 1 1])
    subplot(3,1,2), hold on, grid on, plot(tplot./60, Vel_err_LVLH(:,2), 'Color', 0.6*[1 1 1])
    subplot(3,1,3), hold on, grid on, plot(tplot./60, Vel_err_LVLH(:,3), 'Color', 0.6*[1 1 1])
    
    
    %     n = 60;
    %     xtrue_short = xtrue(:,1:n:end);
    %
    %     xhatminus(:,1) = x;
    %     xhatplus(:,1) = x;
    %     PlotP(1,:) = reshape(P,36,1);
    %
    %     % % Inputing initial camera quaternion
    %     % thetacam = 90*pi/180;
    %     % e_thetacam = [0,1,0]/norm([0,1,0]);
    %     % qcam = [sin(thetacam/2)*e_thetacam cos(thetacam/2)];
    %
    %     z_true = z(1:6)';
    %     err = [];
    
end

% maxstd_pos_x = max(max(sqrt(P_LVLH(:,1))));
% maxstd_pos_y = max(max(sqrt(P_LVLH(:,8))));
% maxstd_pos_z = max(max(sqrt(P_LVLH(:,15))));

maxstd_x = max(max(sqrt(P_LVLH(:,[1 22]))));
maxstd_y = max(max(sqrt(P_LVLH(:,[8 29]))));
maxstd_z = max(max(sqrt(P_LVLH(:,[15 36]))));

figure(1)
subplot(3,1,1), hold on, grid on,
plot(tplot./60, sqrt(zplot_est(:,7)), 'k','Linewidth',1.5)
plot(tplot./60, -sqrt(zplot_est(:,7)), 'k','Linewidth',1.5)

subplot(3,1,2), hold on, grid on,
plot(tplot./60, sqrt(zplot_est(:,14)), 'k','Linewidth',1.5)
plot(tplot./60, -sqrt(zplot_est(:,14)), 'k','Linewidth',1.5)

subplot(3,1,3), hold on, grid on,
plot(tplot./60, sqrt(zplot_est(:,21)), 'k','Linewidth',1.5)
plot(tplot./60, -sqrt(zplot_est(:,21)), 'k','Linewidth',1.5)


figure(2)
subplot(3,1,1), hold on, grid on,
plot(tplot./60, sqrt(zplot_est(:,28)), 'k','Linewidth',1.5)
plot(tplot./60, -sqrt(zplot_est(:,28)),'k','Linewidth',1.5)

subplot(3,1,2), hold on, grid on,
plot(tplot./60, sqrt(zplot_est(:,35)), 'k','Linewidth',1.5)
plot(tplot./60, -sqrt(zplot_est(:,35)), 'k','Linewidth',1.5)

subplot(3,1,3), hold on, grid on,
plot(tplot./60, sqrt(zplot_est(:,42)), 'k','Linewidth',1.5)
plot(tplot./60, -sqrt(zplot_est(:,42)), 'k','Linewidth',1.5)

figure(3)
subplot(3,1,1), hold on, grid on,
plot(tplot./60, 3*sqrt(P_LVLH(:,1)), 'k','Linewidth',1.5)
plot(tplot./60, -3*sqrt(P_LVLH(:,1)), 'k','Linewidth',1.5)
% title('Position Error [m]')
ylabel('x')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
axis([0 tplot(end)/60 -3*maxstd_x 3*maxstd_x])

subplot(3,1,2), hold on, grid on,
plot(tplot./60, 3*sqrt(P_LVLH(:,8)), 'k','Linewidth',1.5)
plot(tplot./60, -3*sqrt(P_LVLH(:,8)), 'k','Linewidth',1.5)
ylabel('y')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
axis([0 tplot(end)/60 -3*maxstd_y 3*maxstd_y])

subplot(3,1,3), hold on, grid on,
plot(tplot./60, 3*sqrt(P_LVLH(:,15)), 'k','Linewidth',1.5)
plot(tplot./60, -3*sqrt(P_LVLH(:,15)), 'k','Linewidth',1.5)
ylabel('z')
xlabel('Time [min]') 
set(findall(gcf,'-property','FontSize'),'FontSize',16)
axis([0 tplot(end)/60 -3*maxstd_z 3*maxstd_z])

figure(4)
subplot(3,1,1), hold on, grid on,
plot(tplot./60, 3*sqrt(P_LVLH(:,22)), 'k','Linewidth',1.5)
plot(tplot./60, -3*sqrt(P_LVLH(:,22)), 'k','Linewidth',1.5)
% title('Velocity Error [mm/s]')
ylabel('x')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
axis([0 tplot(end)/60 -3*maxstd_x 3*maxstd_x])

subplot(3,1,2), hold on, grid on,
plot(tplot./60, 3*sqrt(P_LVLH(:,29)), 'k','Linewidth',1.5)
plot(tplot./60, -3*sqrt(P_LVLH(:,29)), 'k','Linewidth',1.5)
ylabel('y')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
axis([0 tplot(end)/60 -3*maxstd_y 3*maxstd_y])

subplot(3,1,3), hold on, grid on,
plot(tplot./60, 3*sqrt(P_LVLH(:,36)), 'k','Linewidth',1.5)
plot(tplot./60, -3*sqrt(P_LVLH(:,36)), 'k','Linewidth',1.5)
ylabel('z')
xlabel('Time [min]') 
set(findall(gcf,'-property','FontSize'),'FontSize',16)
axis([0 tplot(end)/60 -3*maxstd_z 3*maxstd_z])

% figure(5)
% plot(tplot, sqrt(Pmag)), grid on

