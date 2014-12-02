%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %% 
% Title: q_to_T                                                           %
% Description: Receive a 4x1 quaternion and convert it to a 3x3 rotation  %
%              matrix.                                                    %
% Input: 4x1 quaternion defined by [3x1 vector, 1x1 scalar]               %
% Output: 3x3 Rotation matrix                                             % 
%                                                                         %
% Developed by: ASEL Lab, WVU                                             % 
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [T] =q_to_T(q)

Qv=[q(1);q(2);q(3)];
Qs=q(4);

T=(Qs^2-Qv'*Qv)*eye(3)+2*(Qv*Qv')-2*Qs*crossmat(Qv);

end