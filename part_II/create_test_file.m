%% 6.438 FALL 2015
%% create_test_file.m
%   create test shotgun data

clc; clear; close all;

% feel free to modify parameters to generate
% various types of source data

% # parameters for simulating data
p = 0.4; % strength of stationarity of source chain 
       % (0: identity transition (same state), inf: uniform transition to diff state)
num_shotgun = 10; % number of shotgun reads
l_shotgun = 100; % length of shotgun reads
deg_max_shotgun = 1; % max degree of overlap of a read with its adjacent one
deg_min_shotgun = 1; % min degree of overlap of a read with its adjacent one
noise_shotgun = 0.01; % noise level of shotgun reads

% # parameters for compression
r = 0.4; % coding rate (= k/n)


% / / / / / / / / / / / / / / / /
% / / / do not modify below / / /
addpath(genpath([pwd '/supplementary_functions_part_II'])); %
res = 10;
[s_whole, psi_source] = sample_genome_sequence(res*l_shotgun*num_shotgun,p);
s_shotgun = sample_shotgun_sequence(...
    s_whole,num_shotgun,l_shotgun,...
    deg_max_shotgun,deg_min_shotgun,noise_shotgun);
s_shotgun_trans = s_shotgun';
s = s_shotgun_trans(:); % source data 
    % (vectorized shotgun reads)
sb = vals_8to1_gray(uint8(s), 2); % source data (in binary form)
fprintf('# non-zero elements in LDPC code = ');
[x,H] = encode_binary_sequence(sb,r); % encode source data

% save source and compressed data that will be used for your algo
% data is called in project_part_II.m
save('test_file', ...
    's','sb','x','H','psi_source', ...
    'num_shotgun','l_shotgun','r');
fprintf('saved compressed dataset in test_file ! \n');
