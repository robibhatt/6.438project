function [seq,M] = sample_genome_sequence(n,st)
% (do not modify function)
%
% INPUT:
%   n   - length of whole genome
%   st  - strength of chain
% OUTPUT:
%   seq - single-read
%   M   - transition matrix

X = 4; % ATCG
bits = log2(X);
% //
bleed_src = st;
% start with a transition matrix
prow = normcdf(.5:X-1+.5, 0, bleed_src) - normcdf(-.5:X-1-.5, 0, bleed_src);
epot_template_src = toeplitz(prow);
epot_template_src = ...
    epot_template_src ./ repmat(sum(epot_template_src, 1), X, 1);
% this is normalized as a transition matrix...
% find the stationary distribution
[sta,d] = eig(epot_template_src);  % get eigenvector
[ev,ind] = max(diag(d));  
sta = sta(:,ind) / sum(sta(:,ind));
% draw source
seq(1) = find(rand < cumsum(sta), 1) - 1;
for i=2:n
    seq(i) = find(rand < cumsum(epot_template_src(:,seq(i-1)+1)), 1) - 1;
end    
M = epot_template_src;
%
% compute entropy
hcond = zeros(1,X);
for k=1:X
    hcond(k) = entropy(epot_template_src(:,k));
end
h = hcond * sta / bits;
fprintf(['entropy of chain = ' num2str(h) '\n']);

end

function H = entropy(p)
if (length(p) == 1)
    if p==0 || p==1
        H = 0;
        return;
    else
        p = [p, 1-p];
    end
end
if abs(sum(p) - 1) <= 20*eps  % <-- this is a kludge
    temp = p.*log2(p);
    temp(isnan(temp)) = 0;
    H = -sum(temp);
else
    H = zeros(1,length(p));
    for i=1:length(H)
        H(i) = entropy(p(i));
    end
end
end