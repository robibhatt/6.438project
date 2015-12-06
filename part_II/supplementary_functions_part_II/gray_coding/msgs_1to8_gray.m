function M_out = msgs_1to8_gray(M_in, height, width)

% converts 1bit messages on a single vector to 8bit messages on an image
% M_in: (heigh*width*bits) x 2
% M_out: height x width x M

bits = size(M_in,1)/(height*width);
M = 2^bits;

assert(size(M_in,1) == height * width * bits);
temp = reshape(M_in, [bits, height, width, 2]);
M_out = ones(height, width, M);

for i=1:bits
    index_mask = bitget([0:M-1], i); % i-th bit of index
    M_out(:,:,index_mask == 0) = M_out(:,:,index_mask == 0) .* ...
        repmat(reshape(temp(i,:,:,1),height,width),[1,1,M/2]);
    M_out(:,:,index_mask == 1) = M_out(:,:,index_mask == 1) .* ...
        repmat(reshape(temp(i,:,:,2),height,width),[1,1,M/2]);
end

% gray coding
M_out(:,:,gray2bin(0:M-1)+1) = M_out;

if width == 1
    M_out = squeeze(M_out);
end