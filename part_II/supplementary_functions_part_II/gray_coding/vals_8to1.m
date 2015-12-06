function s = vals_8to1(Sb, bits)

% converts 8bit (grayscale) images to 1bit data vectors
if(~exist('bits'))
    bits = 8;
end

[height, width] = size(Sb);
s = zeros(height*width*bits,1);
for i=1:bits
    temp = bitget(Sb,i);
    s(i:bits:end) = temp(:);
end