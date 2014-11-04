clc, clear
close all

load coord.mat;
load scdata.mat

rsc = Position(:,1:3); % rename known spacecraft position vector
Vel = Velocity(:,1:3); % convert known velocity vectors from m/s
xtrue = [rsc'; Vel']; %

% Define constants
sigma_pos = 100; % 100 m
sigma_vel = .001; % 1 m/s
sigma_theta = 0.01*(pi/180); % 0.01 deg (converted to rad)
T_IC = eye(3,3); % rotation matrix from camera to inertial frame, assuming camera is nadir
options = odeset('abstol',1e-8,'reltol',1e-8); % set tolerances for ode45

% rsc1 = latlong2ECEF(llh(1,1),llh(1,2),llh(1,3))';
%
theta = (7.292e-5);%*(t(k+1)-t(k));
T_ei = [cos(theta) -sin(theta) 0;
    sin(theta) cos(theta) 0;
    0 0 1];
% rsc = T_ei*rsc1;

rsc = xtrue(1:3,1) + 0*sigma_pos*randn(3,1); % use initial position with added noise
rscdot = xtrue(4:6,1) + 0*sigma_vel*randn(3,1); % use initial velocity with added noise
x = [rsc; rscdot]; % create initial state matrix
P = [sigma_pos^2*eye(3), zeros(3); zeros(3), sigma_vel^2*eye(3)]; % create initial P matrix
z = [x; reshape(P,36,1)]; % compile state and covariance into column vector

n = 60;
for k_old =1:n:length(t)
    
    if k_old+1 > length(t)
        break
    end
    
    if k_old == 1
        k = 1;
    else
        k = ((k_old-1)/n)+1;
    end
    
    % Landmark Observation
    [lat_obs,long_obs,height_obs] =  ECEF2latlong(Position(k,1),Position(k_old,2),Position(k_old,3));
%     llh = [lat' long' height'];
      
    % Integrate State and Covariance from tk-1 to tk
    [t1, z1] = ode45(@integrate,t(k):1:t(k+1),z,options); % integrate z to find xdot and Pdot
    z = z1(end,:);
    
    % Check error between true and estimate
    xest(:,k) = [z(1) z(2) z(3) z(4) z(5) z(6)]';
    err_pos(:,k) = xtrue(1:6,k+1)-xest(:,k);

    % Obtain Estimated State Values
    xhatminus(:,k) = [z(1:6)]';
    Pminus = reshape(z(:,7:42),6,6); % reshape into 6x6 matrix
    
    %%% Update Measurement %%%
    
    theta = (7.292e-5)*(t(k_old+1)-t(k_old));
    T_ei = [cos(theta) -sin(theta) 0;
        sin(theta) cos(theta) 0;
        0 0 1];
    [ox, oy, oz] = geo2ECI(lat_obs,long_obs,0); % vector of landmark
    %     o = [Position(k_old,1),Position(k_old,2),Position(k_old,3)]';
    o = [ox; oy; oz];
    
    s = o-xhatminus(1:3,k); % find LOS vector
    s_norm = norm(s); % normalize LOS vector
    eI = s/s_norm; % find unit vector in direction of inertial frame
    
    dv = sigma_theta*randn(3,1); % generate random noise for position vector
    dT = v_to_T(dv); % rotation matrix based on random noise in position vector of camera
    h = dT*T_IC*(o-x(1:3))/(sqrt((o-x(1:3))'*(o-x(1:3)))); % map LOS vector to spacecraft-object position vector
    y = h+dv; % find estimate of position measurement

    % Rotate from body to inertial frame
    ex = -xhatminus(1:3,k)/norm(xhatminus(1:3,k));
    ez = cross(ex,xhatminus(4:6,k)/norm(xhatminus(4:6,k)));
    ey = cross(ez,ex);
    ez = ez/norm(ez);
    ey = ey/norm(ey);
    T_BI = [ex'; ey'; ez'];
    
    H = T_BI*T_IC*(1/s_norm)*[eI*eI'-eye(3) zeros(3)]; % measurement sensitivity matrix
    
    R = 0.0001*sigma_theta^2*(eye(3)-eI*eI'); % measurement covariance matrix
    
    K = (Pminus*H')/((H*Pminus*H'+R+(0.5*trace(R))*(eI*eI'))); % find optimal Kalman gain
    
    xhatplus(:,k) = xhatminus(:,k) + K*(y-h); % state update
    Pplus = (eye(6)-K*H)*Pminus*(eye(6)-K*H)' + K*R*K'; % covariance update
    
    PlotP(k,:) = reshape(Pplus,36,1); % reshape P into vector
    errX(:,k) = [xhatplus(1:3,k)-xtrue(1:3,k+1); (xhatplus(4:6,k)-xtrue(4:6,k+1))]; % Plot error
    
    x = xhatplus(:,k); % update state measurement
    P = Pplus; % update covariance measurement
    z = [x; reshape(P,36,1)]; % compile state and covariance into column vector
    
end

% errX = err_pos;

% figure
% hold on
% plot(xtrue(1,:),'g')
% plot(xest(1,:),'go')
% plot(xtrue(2,:),'b')
% plot(xest(2,:),'bo')
% plot(xtrue(3,:),'r')
% plot(xest(3,:),'ro')
figure
subplot(3,1,1),plot(errX(1,:))
title('Error in X Position (m)')
xlabel('Time (min)')
ylabel('Magnitude (m)')

subplot(3,1,2),plot(errX(2,:))
title('Error in Y Position (m)')
xlabel('Time (min)')
ylabel('Magnitude (m)')

subplot(3,1,3),plot(errX(3,:))
title('Error in Z Position (m)')
xlabel('Time (min)')
ylabel('Magnitude (m)')

figure
subplot(3,1,1),plot(errX(4,:))
title('Error in X Velocity (m/s)')
xlabel('Time (min)')
ylabel('Magnitude (m/s)')

subplot(3,1,2),plot(errX(5,:))
title('Error in Y Velocity (m/s)')
xlabel('Time (min)')
ylabel('Magnitude (m/s)')

subplot(3,1,3),plot(errX(6,:))
title('Error in Z Velocity (m/s)')
xlabel('Time (min)')
ylabel('Magnitude (m/s)')


