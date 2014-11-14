function Rdot = earthgravity(t,R)

% Earth's gravitational parameter
mu = (1000^3)*398600; % gravitational parameter of earth (m3/s2)
R_norm = norm(R(1:3));
adot = (-mu*R(1:3))/(R_norm^3); % m/s^2
Rdot = [R(4:6)/1000; adot*1000]; % [mm/s mm/s^2]

end