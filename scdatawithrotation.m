% scdatawithrot 
clc, clear
close all

load coord.mat;

% Generate values for spacecraft's state (position, velocity, acceleration)
height = 1000*1000; % height above Earth (m) (assume circular orbit at 500 km)
Rearth = 1000*6378; % radius of Earth (m)
alt = Rearth+height; % altitude of orbit (m)
e = 0; % eccentricity (rad) (circular orbit = 0)
i = 90*(pi/180); % inclination in rad (polar orbit = 90 deg)
w = 0; % right ascension of ascending node (rad)
Omega = 0; % argument of periapsis (rad)
nu = 0; % true anomaly (rad)
t = [0:6*60*60]'; % 12 hours of data - time vector for ode45 (sec)

% Propagate elements for length of time specified in t vector using ode45
[R,V] = elementstoRV(alt,e,i,w,Omega,nu); 
x = [R; V]; % [m, m/s]
options = odeset('abstol',1e-8,'reltol',1e-8);
[t1,Pos] = ode45(@earthgravity,t(1):1:t(end),x,options);
R = Pos(:,1:3); 
V = Pos(:,4:6);

angvel = 7.292e-5; % angular velocity of Earth (rad/s)

% include Earth's rotation in calculations
n = 60;
for m = 1:length(t1)

    if m+1 > length(t1)
        break
    end
    
    if m == 1
        theta(m) = angvel*(t(m+1)-t(m)); % rad
    else
        theta(m) = theta(m-1) + angvel*(t(m+1)-t(m));
    end

rot = [cos(theta(m)) -sin(theta(m)) 0;
        sin(theta(m)) cos(theta(m)) 0;
        0 0 1]
    
R_new(m,:) = rot*R(m,:)';
V_new(m,:) = rot*V(m,:)';

lat_calc(m) = (180/pi)*(asin(R_new(m,3)/norm(R_new(m,:)))); % deg
long_calc(m) = (180/pi)*(atan2(R_new(m,2),R_new(m,1))); % deg

Pos_rot(m,:) = R_new(m,:);
Vel_rot(m,:) = V_new(m,:);

end


% Plot orbit to check
figure, plot3(Pos_rot(:,1),Pos_rot(:,2),Pos_rot(:,3))
xlabel('x')
ylabel('y')
zlabel('z')
grid on
axis equal
% 
figure
plot(long,lat,'k')
hold on
plot(long_calc(1:n:end),lat_calc(1:n:end),'*')
% 
% save data to load later
save('scdatawithrot.mat','t','options','height','Rearth','alt','e','i','w','Omega','nu','Pos_rot','Vel_rot','lat_calc','long_calc')

