%% 6.438 FALL 2015
%% project_part_II.m
%   main file for part 2 of project

clc; clear; close all;

% / / / / / / / / / / / / / / / / / / / / / / / / /
% / / do not modify below for final submission / / 
load('test_file'); % call test_file created by create_test_file.m
addpath(genpath([pwd '/supplementary_functions_part_II'])); %
% / / do not modify above for final submission / / 
% / / / / / / / / / / / / / / / / / / / / / / / / /


% Write your code here for choosing
%   - r_dope: doping rate (scalar)
%     There will be penalization in choosing high doping rates
%     in the final calculation of compression ratio.
%   - temp (optional): row vector storing nonnegative integers
%     e.g. temp = [1,10,100] that you may use for your decoder. 
%     There will be penalization in using this extra information 
%     in the final calculation of compression ratio. 
%     Leave it empty if your decoder doesn't use this

r_dope = 0.14; % some doping rate defined here 
                % (can depend on some variables in test_file)
temp = [];


% / / / / / / / / / / / / / / / / / / / / / / / / /
% / / do not modify below for final submission / / 
[phi_source,phi_code] = dope_source(s,r_dope); % dope potentials
[s_hat] = project_part_II_decoder(...
    x,H,phi_source,phi_code,psi_source,temp); % call decoder
if(isempty(temp)); sc = 0; else sc = numel(de2bi(temp)); end
[rH,cH] = size(H);
r_normalized = (rH+cH*r_dope+sc)/cH; % compute normalized rate
fprintf(['* normalized rate = ' num2str(r_normalized) '\n']); % print rate
if(sum(sum(s-s_hat))==0) % lossless?
   fprintf('lossless compression successful\n'); 
   flag_success = 1;
else
   fprintf('compression is not successful... \n'); 
   flag_success = 0;
end
