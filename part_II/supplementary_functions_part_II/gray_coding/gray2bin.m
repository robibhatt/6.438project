function out = gray2bin(in)

out = in;
mask = bitshift(in, -1);
while(any(mask ~= 0))
    out = bitxor(out, mask);
    mask = bitshift(mask, -1);
end