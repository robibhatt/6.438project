function out = bin2gray(in)

out = bitxor(bitshift(in, -1), in);