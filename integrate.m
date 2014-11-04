function zdot = integrate(t,z)

% Redefine parameters from z matrix
rsc = z(1:3);
rscdot = z(4:6);
P = reshape(z(7:42),6,6);

sigmaQ = 1e-7; % standard deviation of state measurement
Q = [zeros(3,6); zeros(3), sigmaQ^2*eye(3)];

rsc_norm = norm(rsc);

mu = (1000^3)*398600; % m3/s2
rscdoubledot = (-mu/(rsc_norm^3))*rsc; % m/s2
Grr = (-mu/(rsc_norm^3))*eye(3)+(3*mu/(rsc_norm^5))*(rsc')*rsc; % m/s2
% Grr(1,1) = mu*((3*(rsc(1)^2)/(rsc_norm^5))-(1/(rsc_norm)^3));
% Grr(1,2) = mu*(3*rsc(1)*rsc(2)/(rsc_norm^5));
% Grr(1,3) = mu*(3*rsc(1)*rsc(3)/(rsc_norm^5));
% Grr(2,1) = mu*((3*(rsc(1)*rsc(2))/(rsc_norm^5)));
% Grr(2,2) = mu*((3*(rsc(2)^2)/(rsc_norm^5))-(1/(rsc_norm)^3));
% Grr(2,3) = mu*(3*rsc(2)*rsc(3)/(rsc_norm^5));
% Grr(3,1) = mu*(3*rsc(1)*rsc(3)/(rsc_norm^5));
% Grr(3,2) = mu*(3*rsc(2)*rsc(3)/(rsc_norm^5));
% Grr(3,3) = mu*((3*(rsc(3)^2)/(rsc_norm^5))-(1/(rsc_norm)^3));

F = [zeros(3) eye(3); 
     Grr zeros(3)];
Pdot = F*P + P*F' + Q; % calculate covariance for derivatives of rsc

% Compile derivative results into single zdot matrix
zdot(1:3) = rscdot; 
zdot(4:6) = rscdoubledot; % m/s2
zdot(7:42) = reshape(Pdot,36,1);
zdot = zdot';

end


