function Sb = vals_1to8(s, height, width, bits)

if(~exist('bits'))
    bits = 8;
end

Sb = zeros(height, width);

for i=1:bits
    Sb = Sb + 2^(i-1) * reshape(s(i:bits:end), height, width);    
end