%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %% 
% Title: crossmat                                                         %
% Description: Create a skew-symmetric cross-matrix from a 3x1 vector     %
%                                                                         %
% Input: 3x1 vector                                                       %
% Output: 3x3 skew-symmetric cross-matrix such that transpose(CM)= -CM    % 
%                                                                         %
% Developed by: ASEL Lab, WVU                                             % 
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [ CM ] = crossmat( alpha )

CM=[0,-alpha(3),alpha(2);
    alpha(3),0,-alpha(1);
    -alpha(2),alpha(1),0];
end

