function s = vals_8to1_gray(Sb, bits)

% converts 8bit (grayscale) images to 1bit GRAYCODED data vectors
if(~exist('bits'))
    bits = 8;
end

[height, width] = size(Sb);
s = zeros(height*width*bits,1);

% do gray coding
Sb = bin2gray(Sb);

for i=1:bits
    temp = bitget(Sb,i);
    s(i:bits:end) = temp(:);
end