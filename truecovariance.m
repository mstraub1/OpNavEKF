%%
clc, clear
close all

load coord.mat;
load scdata.mat

rsc = Position(:,1:3); % rename known spacecraft position vector
Vel = Velocity(:,1:3); % convert known velocity vectors from m/s
xtrue = [rsc'; Vel']; %

% Define constants
sigma_pos = 500; % 100 m
sigma_vel = .01; % 1 mm/s
sigma_theta = 0.005*(pi/180); % 0.01 deg (converted to rad)
T_IC = eye(3,3); % rotation matrix from camera to inertial frame, assuming camera is nadir
options = odeset('abstol',1e-8,'reltol',1e-8); % set tolerances for ode45

% Create perturbed initial state
rsc = xtrue(1:3,1) + 0*sigma_pos*randn(3,1); % use initial position with added noise
rscdot = xtrue(4:6,1) + 0*sigma_vel*randn(3,1); % use initial velocity with added noise

% Build state vector and covariance [m, mm/s]
x = [rsc; rscdot]; % create initial state matrix
P = [sigma_pos^2*eye(3), zeros(3); zeros(3), (1000*sigma_vel)^2*eye(3)]; % create initial P matrix
z = [x; reshape(P,36,1)]; % compile state and covariance into column vector

n = 60;
xtrue_short = xtrue(:,1:n:end);

xhatminus(:,1) = x;
xhatplus(:,1) = x;
PlotP(1,:) = reshape(P,36,1);

% Inputing initial camera quaternion
thetacam = 90*pi/180;
e_thetacam = [0,1,0]/norm([0,1,0]);
qcam = [sin(thetacam/2)*e_thetacam cos(thetacam/2)];

Pos = z(1:6)';
%%
for mc = 1 % for monte carlo analysis
    mc
    for k_old = 1:n:length(t)
        
        if k_old+1 > length(t)
            break
        end
        
        if k_old == 1
            k = 2;
        else
            k = ((k_old-1)/n)+2;
        end
        
        [t1, z1] = ode45(@integrate,0:60,z,options); % integrate z to find xdot and Pdot
        zminus = z1(end,:);
        z_plot(k,:) = zminus;
        
        % Obtain Estimated State Values
        xhatminus(:,k) = [zminus(1:6)]';
        Pminus = reshape(zminus(:,7:42),6,6); % reshape into 6x6 matrix
        Pminus = 0.5*(Pminus + Pminus');
        
        %%% Update Measurement %%%
        % Find FOV footprint limits in terms of lat/long coordinates
        % Specify Desired FOV and calculate footprint size (GIFOV)
        FOV = 30*(pi/180); % rad
        GIFOV = 2*(alt-Rearth)*tan(FOV/2); % km
        
        % Landmark Observation
        [lat_obs,long_obs,height_obs] =  ECEF2latlong(Position(k,1),Position(k,2),Position(k,3));
        
        delta_lat = GIFOV/110.54e3; % convert to m [1deg = 110.54 km]
        delta_long = GIFOV/(111.32e3*cos(lat_obs*(pi/180))); % convert to m [1deg = 111.32*cos(lat)]
        top = [lat_obs+delta_lat, long_obs+delta_long]; % convert to m [1deg = 111.32*cos(lat)]
        left = [lat_obs+delta_lat, long_obs-delta_long];
        bottom = [lat_obs-delta_lat, long_obs-delta_long];
        right = [lat_obs-delta_lat, long_obs+delta_long];
        box = [top; left; bottom; right];
        %         figure, plot(box(1,1),box(1,2),'*',box(2,1),box(2,2),'*',box(3,1),box(3,2),'*',box(4,1),box(4,2),'*')
        %         dist = acos(sin(lat(k)*(pi/180))*sin(lat(k+1)*(pi/180))+cos(lat(k)*(pi/180))*cos(lat(k+1)*(pi/180))*cos((long(k)-long(k+1))*(pi/180)));
        
        % Check for coastline points within footprint limits
        coordfind = find(lat < max(box(:,1)) & lat > min(box(:,1)) & long < max(box(:,2)) & long > min(box(:,2)));
        if isempty(coordfind) == 0
           if length(coordfind) > 1
                latlong_est = [lat(coordfind(1:2)) long(coordfind(1:2))];
                number = 2;
            else
                latlong_est = [lat(coordfind(1)) long(coordfind(1))];
                number = 1;
            end
            for num = 1:number
                [ox, oy, oz] = geo2ECI(latlong_est(num,1)+10*randn(1),latlong_est(num,2)+10*randn(1),height_obs); % vector of landmark
                o = [ox; oy; oz];
                
                % Compute measurement
                s_true = o-xtrue(1:3,k_old+n);
                eI = s_true/norm(s_true);
                y = T_IC*eI;
                %
                %         dv = sigma_theta*randn(3,1); % generate random noise for position vector
                %         dT = v_to_T(dv); % rotation matrix based on random noise in position vector of camera
                %         y = dT*y_true;
                %
                %         s = o-xhatminus(1:3,k); % find LOS vector
                %         s_norm = norm(s); % normalize LOS vector
                %         eI = s/s_norm; % find unit vector in direction of inertial frame
                h = T_IC*eI;
                
                H = T_IC*(1/norm(s_true))*[eI*eI'-eye(3) zeros(3)]; % measurement sensitivity matrix
                
                R = sigma_theta^2*(eye(3)-eI*eI'); % measurement covariance matrix
                
                K = (Pminus*H')/(H*Pminus*H'+R+(0.5*trace(R))*(eI*eI')); % find optimal Kalman gain
                
                xhatplus(:,k) = xhatminus(:,k) + K*(y-h); % state update
                Pplus = (eye(6)-K*H)*Pminus*(eye(6)-K*H)' + K*R*K'; % covariance update
            end
        else
            fprintf('No Coastline at k = %f. \n',k)
            
            Pplus = Pminus;
            xhatplus = xhatminus;
        end
        
        PlotP(k,:) = reshape(Pplus,36,1); % reshape P into vector
        %         errX(:,k) = [xhatplus(1:3,k)-xtrue_short(1:3,k); (xhatplus(4:6,k)-xtrue_short(4:6,k))]; % Plot error
        
        x = xhatplus(:,k); % update state measurement
        P = Pplus; % update covariance measurement
        z = [x; reshape(P,36,1)]; % compile state and covariance into column vector
        
    end
    %%
    for count = 1:length(PlotP)
        maxstd_pos = max(max(sqrt(z_plot(count,[7,14,21]))));
        maxstd_vel = max(max(sqrt(z_plot(count,[28,35,42]))));
        
        state = z_plot(count,:);
        P2 = reshape(state(:,7:42),6,6);
        P_pos = P2(1:3,1:3);
        P_vel = P2(4:6,4:6);
        
        % Determining rotation matrix from inertial frame to spacecraft frame
        ex = -state(1:3)/norm(state(1:3));
        ez = cross(ex,state(4:6)/norm(state(4:6)));
        ez = ez/norm(ez);
        ey = cross(ez,ex);
        ey = ey/norm(ey);
        T_BI = [ex',ey',ez'];
        
        % Rotating covariances from inertial to spacecraft frame
        P_pos_rot = T_BI'*P_pos*T_BI;
        P_vel_rot = T_BI'*P_vel*T_BI;
        % Rotate state
        %         err_pos_rot(:,count) = T_BI'*errX(1:3,count);
        %         err_vel_rot(:,count) = T_BI'*errX(4:6,count);
        
        P_pos_rot_diag(:,count) = sqrt(diag(P_pos_rot));
        P_vel_rot_diag(:,count) = sqrt(diag(P_vel_rot));
        
    end
    maxstd_pos2 = max(max(P_pos_rot_diag));
    maxstd_vel2 = max(max(P_vel_rot_diag));
    %%
    figure(1)
    subplot(3,1,1)
    hold on
    %     plot(err_pos_rot(1,2:end))
    plot(P_pos_rot_diag(1,2:end),'k')
    plot(-P_pos_rot_diag(1,2:end),'k')
    axis([0 length(PlotP) -maxstd_pos2 maxstd_pos2])
    title('Error in X Position (m)')
    xlabel('Time (min)')
    ylabel('Magnitude (m)')
    
    subplot(3,1,2)
    hold on
    %     plot(err_pos_rot(2,2:end))
    plot(P_pos_rot_diag(2,2:end),'k')
    plot(-P_pos_rot_diag(2,2:end),'k')
    axis([0 length(PlotP) -maxstd_pos2 maxstd_pos2])
    title('Error in Y Position (m)')
    xlabel('Time (min)')
    ylabel('Magnitude (m)')
    
    subplot(3,1,3)
    hold on
    %     plot(err_pos_rot(3,2:end))
    plot(P_pos_rot_diag(3,2:end),'k')
    plot(-P_pos_rot_diag(3,2:end),'k')
    axis([0 length(PlotP) -maxstd_pos2 maxstd_pos2])
    title('Error in Z Position (m)')
    xlabel('Time (min)')
    ylabel('Magnitude (m)')
    
    figure(2)
    subplot(3,1,1)
    hold on
    %     plot(err_vel_rot(1,:))
    plot(P_vel_rot_diag(1,2:end),'k')
    plot(-P_vel_rot_diag(1,2:end),'k')
    axis([0 length(PlotP) -maxstd_vel2 maxstd_vel2])
    title('Error in X Velocity (mm/s)')
    xlabel('Time (min)')
    ylabel('Magnitude (mm/s)')
    
    subplot(3,1,2)
    hold on
    %     plot(err_vel_rot(2,:))
    plot(P_vel_rot_diag(2,2:end),'k')
    plot(-P_vel_rot_diag(2,2:end),'k')
    axis([0 length(PlotP) -maxstd_vel2 maxstd_vel2])
    title('Error in Y Velocity (mm/s)')
    xlabel('Time (min)')
    ylabel('Magnitude (mm/s)')
    
    subplot(3,1,3)
    hold on
    %     plot(err_vel_rot(3,:))
    plot(P_vel_rot_diag(3,2:end),'k')
    plot(-P_vel_rot_diag(3,2:end),'k')
    axis([0 length(PlotP) -maxstd_vel2 maxstd_vel2])
    title('Error in Z Velocity (mm/s)')
    xlabel('Time (min)')
    
    
    
end
