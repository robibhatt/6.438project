function [s_hat] = project_part_II_decoder(...
    x,H,phi_source,phi_code,psi_source,temp)
% [do not modify arguments (input/output) of function]
% INPUT
%   x           - k x 1 compressed data vector 
%   H           - k x n code matrix
%   phi_source  - m x 4 source node potentials from doping
%   phi_code    - n x 2 code node potentials from doping
%   psi_source  - 4 x 4 chain transition matrix (each column sums to 1)
%   temp        - [optional] 1 x L array storing L nonneg integer values
% OUTPUT
%   s_hat       - m x 1 decoded/decompressed result

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
% Implement your decoder (going from x to s_hat) here 






s_hat = zeros(size(H,2)/2,1); % output decoded solution as column vector


end

