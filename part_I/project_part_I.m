%% 6.438 FALL 2015
%% MAIN FILE FOR PART I OF PROJECT

% make sure you have unzipped all the given files in the same dir...

clc; close all; clear;

%
compressed_file_name = 'mcoli_rate_high'; % compressed file name (can modify)
%   - mcoli_rate_high
%   - mcoli_rate_moderate
%   - mcoli_rate_low

%
fprintf('Start testing compression method ... \n');
disp(compressed_file_name);
load(compressed_file_name);
load('mcoli'); % ground truth (fixed)
load('mcoli_code_dope'); % doping parameters (fixed)
addpath(genpath([pwd '/supplementary_functions_part_I'])); % add functions in path
% fetch some dimensions
m = length(s);
[k,n] = size(H);

% some initializations for variables you will use (no modification necessary)
M_from_code = ones(n,2); % binary msgs coming out of code graph
M_to_source = ones(m,4); % M_from_code after bit-to-alphabet conversion to enter source graph
M_from_source = ones(m,4); % msgs coming out of source graph
M_to_code = ones(n,2); % M_from_source after alphabet-to-bit converstion to enter code graph
vector_error = []; % a vector storing error for every ite
s_hat = zeros(m,1); % estimate of source data
s_hat_old = zeros(m,1); % previous estimate of source data
% o_source(source_node, direction, value) = 
% message from source_node to source_node+1 (or -1) depending on direction
% evaluated with value
o_source = ones(m,2,4); 
o_code = [];

% start BP
l = 0;
while(1)
    l = l+1;
    errs = 0;
    fprintf(['Ite num = ' num2str(l) '\n']); % print iteration number
    
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    % [1] CODE GRAPH BP ITERATION
    
    % ************************************************************
    % ************ write your code graph BP code here ************
    % ************************************************************
    % (either directly in here or as function call)
    % input variables:     
    %   M_to_code - msg passed from source to code graph
    %               [size n x 2]
    %   o_code - struct of msgs in code graph that is being updated
    %   x - the compressed data
    %       [size k x 1] (given by us)
    %   H - LDPC matrix 
    %       [size k x n] (given by us)
    %   phi_code - doped node potentials
    %       [size n x 2] (given by us)
    % output variable:
    %   M_from_code - product of incoming code graph msgs at every node
    %                 to be passed to source graph (all after msg update)
    %                 [size n x 2] (convert the LLR message to standard
    %                 message)
    %   o_code - struct of updated msgs in code graph
    
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % [2] BIT-TO-ALPHABET CONVERSION FOR SOURCE GRAPH BP
    % (no modification necessary)
    M_to_source = squeeze(msgs_1to8_gray(M_from_code,1,m)); % [size m x 4]

    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % [3] SOURCE GRAPH BP ITERATION
    
    % **************************************************************
    % ************ write your source graph BP code here ************
    % **************************************************************
    % (either directly in here or as function call)
    % input variables:
    %   M_to_source - msg passed from code to source graph
    %                 [size m x 4]
    %   o_source - struct of msgs in source graph that is being updated
    %              [size m x 2 x 4]
    %   phi_source - doped node potentials 
    %                [size m x 4] (given by us)
    %   psi_source - the representative edge potential matrix 
    %                representing probabilistic transitions in markov chain
    %                from i to i+1 where the row sum is 1
    %                [size 4 x 4] (given by us)
    % output variable:
    %   M_from_source - product of incoming source graph msgs at every node
    %                   to be passed to code graph (all after msg update)
    %                   [size m x 4]
    %   o_source - struct of updated msgs in source graph
    new_o_source = ones(m,2,4);
    for from_node = 1:m
        for to_value = 1:4 %1 to 4 map to 0 to 3 because matlab is annoying
            % forward messages
            if from_node == 1
                message = 0;
                for from_value = 1:4
                    message = message + phi_source(from_node, from_value)*psi_source(from_value, to_value)*M_to_source(from_node, from_value);
                end
                new_o_source(from_node, 1, to_value) = message;
            elseif from_node < m
                message = 0;
                for from_value = 1:4
                    message = message + phi_source(from_node, from_value)*psi_source(from_value, to_value)*o_source(from_node - 1, 1, from_value)*M_to_source(from_node, from_value);
                end
                new_o_source(from_node, 1, to_value) = message;
            end
            % backward messages
            if from_node == m
                message = 0;
                for from_value = 1:4
                    message = message + phi_source(from_node, from_value)*psi_source(to_value, from_value)*M_to_source(from_node, from_value);
                end
                new_o_source(from_node, 2, to_value) = message;
            elseif from_node > 1
                message = 0;
                for from_value = 1:4
                    message = message + phi_source(from_node, from_value)*psi_source(to_value, from_value)*o_source(from_node + 1, 2, from_value)*M_to_source(from_node, from_value);
                end
                new_o_source(from_node, 2, to_value) = message;
            end  
        end
        % normalization
        for direction = 1:2
            total_message = 0;
            for to_value = 1:4
                total_message = total_message + new_o_source(from_node, direction, to_value);
            end
            for to_value = 1:4
                new_o_source(from_node, direction, to_value) = new_o_source(from_node, direction, to_value) / total_message;
            end
        end
    end
    % sets M_from_source
    for node = 1:m
        for value = 1:4
            if node == 1
                M_from_source(node, value) = o_source(node + 1, 2, value);
            elseif node < m
                M_from_source(node, value) = o_source(node - 1, 1, value)*o_source(node + 1, 2, value);
            else
                M_from_source(node, value) = o_source(node-1, 1, value);
            end
        end
        % normalization
        total_message = 0;
        for value = 1:4
            total_message = total_message + M_from_source(node, value);
        end
        for value = 1:4
            M_from_source(node, value) = M_from_source(node,value) / total_message;
        end
    end
    % sets o_source
    for node = 1:m
        for direction = 1:2
            for value = 1:4
                o_source(node, direction, value) = new_o_source(node, direction, value);
            end
        end
    end
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % [4] ALPHABET-TO-BIT CONVERSION FOR CODE GRAPH BP
    % (no modification necessary)
    M_to_code = msgs_8to1_gray(M_from_source); % [size n x 2]
    
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        
    % ***************************************************************
    % **** write your code here for decoding from source msgs *******
    % ***************************************************************
    % (either directly in here or as function call)
    % input variables:
    %   M_from_source - msgs coming into nodes within source graph
    %                   [size m x 4]
    %   M_to_source - msgs coming into nodes from code graph
    %                 [size m x 4]
    %   phi_source - doped source graph node potentials
    %                [size m x 4]
    % output variable:
    %   s_hat - decoded solution using marginal mode
    %          [size m x 1]
    % added a comment
    for index = 1:m
        max = 1;
        value = phi_source(index,1)*M_from_source(index,1)*M_to_source(index,1);
        for guess = 2:4
            new_value = phi_source(index,guess)*M_from_source(index,guess)*M_to_source(index,guess);
            if new_value > value
                max = guess;
                value = new_value;
            end
        end
        s_hat(index) = max-1;
    end
    % ************************************************************
    % ****** write your code here for computing error ************
    % ************************************************************
    % (either directly in here or as function call)
    % input variables:
    %   s_hat - decoded solution
    %           [size m x 1]
    %   s - true source data
    %       [size m x 1] (given by us)    
    % output variable:
    %   errs - abs difference error
    %          [scalar variable]
    errs = 0;
    for index = 1:m
        if s_hat(index) ~= s(index)
            errs = errs + 1;
        end
    end
    errs = errs / m;
    vector_error = [vector_error errs];
    fprintf(['... Error = ' num2str(errs) '\n']);
    % terminate if BP gradient doesn't change
    if(l>1 && sum(abs(s_hat - s_hat_old))<0.5)
        break;  % exit BP loop
    end
    s_hat_old = s_hat; % update solution
end


% ******************************************************************************
% ****** write your code here for plotting vector_error vs. l (iteration) ******
% ******************************************************************************


% note: if vector_error converges to 0 exactly, 
%       then you have sucessfully achieved lossless compression 

