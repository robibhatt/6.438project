function [phi_source,phi_code] = dope_source(s,doperate)
% [do not modify function]

% INPUT
%   s - source data 
%       [m x 1] 
%   rate_dope - doping rate (0 to 1 scalar)
% OUTPUT
%   phi_source - doped node potential of source graph
%                [m x 4]
%   phi_node - doped node potential of code graph
%              [n x 2]

m = length(s);
npot_template = reshape(1/4 * ones(4,1), [1,1,4]);
npot = repmat(npot_template, [1, m, 1]);
D = randsample(m, ceil(doperate*m));   % random doping
for i=1:4
    temp = npot(:,:,i);
    temp(D) = (s(D) == i-1);
    npot(:,:,i) = temp;
end
p_s = msgs_8to1_gray(npot);
phi_source = squeeze(npot);
phi_code = squeeze(p_s);


end

