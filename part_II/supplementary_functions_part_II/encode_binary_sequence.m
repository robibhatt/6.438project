function [x,H] = encode_binary_sequence(sb,r)
% (do not modify function)

% INPUT:
%   sb - binary source data
%   r - compression rate
% OUTPUT:
%   x - compressed data
%   H - LDPC coding matrix in array form

n = length(sb);
k = floor(r*n);
H = matlab.internal.sparse.repairSparse(ldpc_generate(k,n,3,2,100));
         % By Eric's suggestion
x = mod(H*sb, 2); % compressed data

end

