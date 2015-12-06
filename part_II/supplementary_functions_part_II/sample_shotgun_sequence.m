function [seq_shotgun] = ...
    sample_shotgun_sequence(...
    seq,num_shotgun,l_shotgun,...
    deg_max_shotgun,deg_min_shotgun,noise_shotgun...
    )
% (do not modify function)
%
% INPUT:
%   seq - whole genome sequence
%   num_shotgun - number of shotgun reads
%   l_shotgun - length of a shotgun read
%   deg_max_shotgun - max degree of overlap between a read and its adjacent (0 to 1)
%   deg_min_shotgun - min degree of overlap between a read and its adjacent (0 to 1)
%   noise_shotgun - level of noise in shotgun seq (0 to 1)
% OUTPUT:
%   seq_shotgun - array of noisy shotgun reads
%                 [num_shotgun x l_shotgun]

rand('state',sum(100*clock));
m = length(seq);
if((m+num_shotgun)<num_shotgun*l_shotgun)
    error('whole genome length is too short')
end
if(deg_max_shotgun<deg_min_shotgun)
    error('bad overlap parameters');
end
seq_shotgun = zeros(num_shotgun,l_shotgun);
offsets = round(rand(1)*(length(seq)-(num_shotgun*l_shotgun+1)));
inds = offsets+(1:l_shotgun);
for i=randperm(num_shotgun)
    seq_raw = seq(inds);
    seq_noisy_index = rand(1,l_shotgun)<noise_shotgun;
    seq_noisy = seq_raw;
    for l=find(seq_noisy_index) % i.i.d.
        bases = [0,1,2,3];
        val_l = seq_raw(l);
        bases(bases==val_l) = [];
        new_ind = ceil(rand(1)*3);
        seq_noisy(l) = bases(new_ind);
    end
    seq_shotgun(i,:) = seq_noisy;
    inds = round(...
        inds+l_shotgun... % shift to new index
        -l_shotgun*((deg_max_shotgun-deg_min_shotgun)*rand(1)+deg_min_shotgun)...
        );
end


end

