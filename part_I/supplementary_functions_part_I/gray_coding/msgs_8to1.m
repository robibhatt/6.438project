function M_out = msgs_8to1(M_in)

% converts 8bit messages on an image to 1bit messages on a single vector
% M_in: height x width x M
% M_out: (height*width*bits) x 2

nd = ndims(M_in);
siz = size(M_in);
N = prod(siz(1:nd-1));
M = siz(end);
M_in = reshape(M_in, N, M);
bits = log2(M);
assert(prod(size(M_in)) == N*M);
M_out = zeros(N*bits, 2);

for i=1:bits
    index_mask = bitget([0:M-1], i); % i-th bit of index
    temp1 = sum(M_in(:,index_mask == 0), 2); % sumprob of 0's on i-th bit of index
    temp2 = sum(M_in(:,index_mask == 1), 2); % sumprob of 1's on i-th bit of index
    M_out(i:bits:end,1) = temp1(:);
    M_out(i:bits:end,2) = temp2(:);
end