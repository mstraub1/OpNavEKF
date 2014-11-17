%% 
clc, clear
close all

load coord.mat;
load scdata.mat

rsc = Position(:,1:3); % rename known spacecraft position vector
Vel = Velocity(:,1:3); % convert known velocity vectors from m/s
xtrue = [rsc'; Vel']; %

% Define constants
sigma_pos = 100; % 100 m
sigma_vel = .001; % 1 mm/s
sigma_theta = 0.005*(pi/180); % 0.01 deg (converted to rad)
T_IC = eye(3,3); % rotation matrix from camera to inertial frame, assuming camera is nadir
options = odeset('abstol',1e-8,'reltol',1e-8); % set tolerances for ode45

% Create perturbed initial state
rsc = xtrue(1:3,1) + sigma_pos*randn(3,1); % use initial position with added noise
rscdot = xtrue(4:6,1) + sigma_vel*randn(3,1); % use initial velocity with added noise

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
              
        % Integrate State and Covariance from tk-1 to tk
%         [t2,Posout] = ode45(@earthgravity,0:60,Pos(end,:)',options);
%         if k_old == 1
%             T2 = t2;
%             Pos = Posout;
%         else
%             T2 = [T2;t2+T2(end)];
%             Pos = [Pos;Posout];
%         end
        
        [t1, z1] = ode45(@integrate,0:60,z,options); % integrate z to find xdot and Pdot
        
%         if k_old == 1
%             T = t1;
%             ZZ=z1;
%             ERR_Pos = sqrt(diag((Pos(:,1:3)-ZZ(:,1:3))*(Pos(:,1:3)-ZZ(:,1:3))'));
%         else
%             T = [T;t1+T(end)];
%             ZZ = [ZZ;z1];
%             ERR_Pos = [ERR_Pos;sqrt(diag((Posout(:,1:3)-z1(:,1:3))*(Posout(:,1:3)-z1(:,1:3))'))];
%         end
%         
        zminus = z1(end,:);
        z_plot(k,:) = zminus;
        
        % Check error between true and estimate
        %         xest(:,k) = [zminus(1) zminus(2) zminus(3) zminus(4) zminus(5) zminus(6)]';
        %         err_pos(:,k) = xtrue_short(1:6,k)-xest(:,k);
        
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
        %     llh = [lat' long' height'];

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
            latlong_est = [lat(coordfind(1)) long(coordfind(1))];
            [ox, oy, oz] = geo2ECI(latlong_est(1)+10*randn(1),latlong_est(2)+10*randn(1),height_obs); % vector of landmark
            o = [ox; oy; oz];
        else
%             fprintf('No Coastline at k = %f. \n',k)
            [ox, oy, oz] = geo2ECI(lat_obs+10*randn(1),long_obs+10*randn(1),height_obs); % vector of landmark
            o = [ox; oy; oz];
            % continue
        end
        
%         [ox, oy, oz] = geo2ECI(lat_obs+10*randn(1),long_obs+10*randn(1),height_obs); % vector of landmark
%         %     o = [Position(k_old,1),Position(k_old,2),Position(k_old,3)]';
%         o = [ox; oy; oz];

        % Rotate from inertial frame to spacecraft frame
        %ex = -xhatminus(1:3,k)/norm(xhatminus(1:3,k));
        %ez = cross(ex,xhatminus(4:6,k)/norm(xhatminus(4:6,k)));
        %ey = cross(ez,ex);
        %ez = ez/norm(ez);
        %ey = ey/norm(ey);
        %T_BI = [ex'; ey'; ez'];
        %T_IC = q_to_T(qcam)'*T_BI';
        
        % Compute measurement
        s_true = o-xtrue(1:3,k_old+n);
        eI_true = s_true/norm(s_true);
        y_true = T_IC*eI_true;
        
        dv = sigma_theta*randn(3,1); % generate random noise for position vector
        dT = v_to_T(dv); % rotation matrix based on random noise in position vector of camera
        y = dT*y_true;
        
        s = o-xhatminus(1:3,k); % find LOS vector
        s_norm = norm(s); % normalize LOS vector
        eI = s/s_norm; % find unit vector in direction of inertial frame
        h = T_IC*eI;
        
        H = T_IC*(1/s_norm)*[eI*eI'-eye(3) zeros(3)]; % measurement sensitivity matrix
        
        R = sigma_theta^2*(eye(3)-eI*eI'); % measurement covariance matrix
        
        K = (Pminus*H')/(H*Pminus*H'+R+(0.5*trace(R))*(eI*eI')); % find optimal Kalman gain
        
%         -xhatminus(:,k)+xtrue_short(:,k)
%         K*(y-h)
        xhatplus(:,k) = xhatminus(:,k) + K*(y-h); % state update
        Pplus = (eye(6)-K*H)*Pminus*(eye(6)-K*H)' + K*R*K'; % covariance update
        
%         [UU,SS,VV] = svd(Pplus);
%         SS_keep(:,k) = diag(SS);
%         
        PlotP(k,:) = reshape(Pplus,36,1); % reshape P into vector
        errX(:,k) = [xhatplus(1:3,k)-xtrue_short(1:3,k); (xhatplus(4:6,k)-xtrue_short(4:6,k))]; % Plot error
%         innovkeep(:,k) = y-h;
%         ykeep(:,k) = y;
        
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
        err_pos_rot(:,count) = T_BI'*errX(1:3,count);
        err_vel_rot(:,count) = T_BI'*errX(4:6,count);

        P_pos_rot_diag(:,count) = sqrt(diag(P_pos_rot));
        P_vel_rot_diag(:,count) = sqrt(diag(P_vel_rot));

    end
    maxstd_pos2 = max(max(P_pos_rot_diag));
    maxstd_vel2 = max(max(P_vel_rot_diag));
    %%
    
    figure(1)
    subplot(3,1,1)
    hold on
    plot(err_pos_rot(1,2:end))
    plot(P_pos_rot_diag(1,2:end),'k')
    plot(-P_pos_rot_diag(1,2:end),'k')
    axis([0 length(PlotP) -maxstd_pos2 maxstd_pos2])
    title('Error in X Position (m)')
    xlabel('Time (min)')
    ylabel('Magnitude (m)')
    
    subplot(3,1,2)
    hold on
    plot(err_pos_rot(2,2:end))
    plot(P_pos_rot_diag(2,2:end),'k')
    plot(-P_pos_rot_diag(2,2:end),'k')
    axis([0 length(PlotP) -maxstd_pos2 maxstd_pos2])
    title('Error in Y Position (m)')
    xlabel('Time (min)')
    ylabel('Magnitude (m)')
    
    subplot(3,1,3)
    hold on
    plot(err_pos_rot(3,2:end))
    plot(P_pos_rot_diag(3,2:end),'k')
    plot(-P_pos_rot_diag(3,2:end),'k')
    axis([0 length(PlotP) -maxstd_pos2 maxstd_pos2])
    title('Error in Z Position (m)')
    xlabel('Time (min)')
    ylabel('Magnitude (m)')
    
    figure(2)
    subplot(3,1,1)
    hold on
    plot(err_vel_rot(1,:))
    plot(P_vel_rot_diag(1,2:end),'k')
    plot(-P_vel_rot_diag(1,2:end),'k')
    axis([0 length(PlotP) -maxstd_vel2 maxstd_vel2])
    title('Error in X Velocity (mm/s)')
    xlabel('Time (min)')
    ylabel('Magnitude (mm/s)')
    
    subplot(3,1,2)
    hold on
    plot(err_vel_rot(2,:))
    plot(P_vel_rot_diag(2,2:end),'k')
    plot(-P_vel_rot_diag(2,2:end),'k')
    axis([0 length(PlotP) -maxstd_vel2 maxstd_vel2])
    title('Error in Y Velocity (mm/s)')
    xlabel('Time (min)')
    ylabel('Magnitude (mm/s)')
    
    subplot(3,1,3)
    hold on
    plot(err_vel_rot(3,:))
    plot(P_vel_rot_diag(3,2:end),'k')
    plot(-P_vel_rot_diag(3,2:end),'k')
    axis([0 length(PlotP) -maxstd_vel2 maxstd_vel2])
    title('Error in Z Velocity (mm/s)')
    xlabel('Time (min)')
    ylabel('Magnitude (mm/s)')
    %%
%     subplot(3,2,1)
%     hold on,plot(xhatminus(1,:),'bo'),plot(xhatplus(1,:),'r*'),plot(xtrue(1,2:n:end),'kx')
%     title('X Position')
%     
%     subplot(3,2,3)
%     hold on,plot(xhatminus(2,:),'bo'),plot(xhatplus(2,:),'r*'),plot(xtrue(2,2:n:end),'kx')
%     title('Y Position')
%     
%     subplot(3,2,5)
%     hold on,plot(xhatminus(3,:),'bo'),plot(xhatplus(3,:),'r*'),plot(xtrue(3,2:n:end),'kx')
%     title('Z Position')
%     
%     subplot(3,2,2)
%     hold on,plot(xhatminus(4,:),'bo'),plot(xhatplus(4,:),'r*'),plot(xtrue(4,2:n:end),'kx')
%     title('X Velocity')
%     
%     subplot(3,2,4)
%     hold on,plot(xhatminus(5,:),'bo'),plot(xhatplus(5,:),'r*'),plot(xtrue(5,2:n:end),'kx')
%     title('Y Velocity')
%     
%     subplot(3,2,6)
%     hold on,plot(xhatminus(6,:),'bo'),plot(xhatplus(6,:),'r*'),plot(xtrue(6,2:n:end),'kx')
%     title('Z Velocity')
%     legend('xhatminus','xhatplus','xtrue','Location','southoutside','Orientation','horizontal')
%     
    % Rotated Plots

%     % Non-rotated Plots
%     figure(1)
%     subplot(3,1,1)
%     hold on
%     plot(errX(1,:))
%     plot(sqrt(PlotP(:,1)),'k')
%     plot(-sqrt(PlotP(:,1)),'k')
%     title('Error in X Position (m)')
%     xlabel('Time (min)')
%     ylabel('Magnitude (m)')
%     
%     subplot(3,1,2)
%     hold on
%     plot(errX(2,:))
%     plot(sqrt(PlotP(:,8)),'k')
%     plot(-sqrt(PlotP(:,8)),'k')
%     title('Error in Y Position (m)')
%     xlabel('Time (min)')
%     ylabel('Magnitude (m)')
%     
%     subplot(3,1,3)
%     hold on
%     plot(errX(3,:))
%     plot(sqrt(PlotP(:,15)),'k')
%     plot(-sqrt(PlotP(:,15)),'k')
%     title('Error in Z Position (m)')
%     xlabel('Time (min)')
%     ylabel('Magnitude (m)')
%     
%     figure(2)
%     subplot(3,1,1)
%     hold on
%     plot(errX(4,:))
%     plot(sqrt(PlotP(:,22))/1000,'k')
%     plot(-sqrt(PlotP(:,22))/1000,'k')
%     title('Error in X Velocity (m/s)')
%     xlabel('Time (min)')
%     ylabel('Magnitude (m/s)')
%     
%     subplot(3,1,2)
%     hold on
%     plot(errX(5,:))
%     plot(sqrt(PlotP(:,29))/1000,'k')
%     plot(-sqrt(PlotP(:,29))/1000,'k')
%     title('Error in Y Velocity (m/s)')
%     xlabel('Time (min)')
%     ylabel('Magnitude (m/s)')
%     
%     subplot(3,1,3)
%     hold on
%     plot(errX(6,:))
%     plot(sqrt(PlotP(:,36))/1000,'k')
%     plot(-sqrt(PlotP(:,36))/1000,'k')
%     title('Error in Z Velocity (m/s)')
%     xlabel('Time (min)')
%     ylabel('Magnitude (m/s)')
    
%     sig3 = 3*sqrt(ZZ(:,7)+ZZ(:,14)+ZZ(:,21));
%     figure
%     plot(T,sig3,'r-')
%     hold on
%     plot(T2,ERR_Pos,'k-')
%     axis([0 T(end),0,max(sig3)])
%     
%     figure
%     plot(SS_keep(1,:),'k'), hold on
%     plot(SS_keep(2,:))
%     plot(SS_keep(3,:))
%     plot(SS_keep(4,:))
%     plot(SS_keep(5,:))
%     plot(SS_keep(6,:))
%     grid on
       
end