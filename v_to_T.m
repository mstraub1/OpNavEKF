%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %% 
% Title: v_to_T                                                           %
% Description: Receive a 3x1 rotation vector defined by Euler angles.     %
%              Convert the 3x1 vector to a 3x3 rotation matrix.           %
% Input: 3x1 Euler angle rotation vector in radians                       %
% Output: 3x3 Rotation matrix                                             % 
%                                                                         %
% Developed by: ASEL Lab, WVU                                             % 
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


function [ T ] = v_to_T( v )


if norm(v) < 1e-14
    T = eye(3);
else
    theta = norm(v);
    etheta(:,1) = v/theta;
    
    T = cos(theta)*eye(3) + (1-cos(theta))*(etheta*etheta') - sin(theta)*crossmat(etheta);
end

end

