Below is a description of given .m files.
Check the header of these files as well for full description of input/output.

= = = = = = = = = = = = = = = = = = = = = = 
[main files]

There are three .m files that you should take a look at in sequential order

1. create_test_file.m
    - Creates noisy shotgun data and saves it in test_file.mat.
    The necessary function calls are already given, and you just have to
    tweak parameters to create new source data for your own testing.
    However, since one of the test datasets used for grading will have the
    default set of parameters in create_test_file.m, 
        p = 0.4; % strength of stationarity of source chain,
        num_shotgun = 10; % number of shotgun reads,
        l_shotgun = 100; % length of shotgun reads,
        deg_max_shotgun = 1; % max degree of overlap of shotgun reads,
        deg_min_shotgun = 0.6; % min degree of overlap of shotgun reads,
        noise_shotgun = 0.01; % noise level of shotgun reads,
    try to optimize your lossless compression algorithm to achieve
    lowest rate on this parameter set first before trying out different
    parameter sets. The values you could try out to further
    optimize your algo can lie in the following range
        p = 0.3 ~ 0.7
        num_shotgun = 3 ~ 20
        l_shotgun = 50 ~ 200
        deg_max_shotgun = 0.5 ~ 1
        deg_min_shotgun = 0 ~ min(0.6,deg_max_shotgun)
        noise_shotgun = 0 ~ 0.05

2. project_part_II.m
    - Your main file that loads source and encoded data created
    by create_test_file.m, and calls project_part_II_decoder.m
    to perform LBP

3. project_part_II_decoder.m
    - The main function that performs the decoding using
    message passing algorithms on code and source graphs.
    You should explore several graphical model structures (especially for the source graph),
    using your knowledge from Part I of the project.
    (Called by project_part_II.m)
      

= = = = = = = = = = = = = = = = = = = = = = 
[supplementary files]

Supplementary functions in older /supplementary_functions_part_II
(do not modify them since in the end we will be grading using the same
functions...)

encode_binary_sequence.m
    - creates LDPC code and encodes source data (in binary form).
    For creating the code, calls .mex files of ldpc_generate in /ldpc.
    In the case your OS doesn't doesn't support the given .mex,
    we will be releasing a .zip file where you can look-up H directly
    from memory. Also, raw ldpc_generate.c code (copyright 1999 IK)
    is given in /ldpc
    (Called by project_part_II.m)

dope_source.m
    - creates randomized doping potentials using source data
    (Called by project_part_II.m)

sample_shotgun_sequence.m
    - generates a set of noisy shotgun reads
    (Called by create_test_file.m)
    
Also, we have the same gray-coding supplementary functions
used for Part I of the project in /gray_coding

Finally, at least for the Part II code submission,
do not modify supplementary functions or script parts
that have claims 'do not modify',
since we will be grading using the same supplementary code.
With your submitted code of project_part_II.m and project_part_II_decoder.m,
we will be generating several shotgun datasets and computing the lowest rate
(=highest compression ratio) achieved by your algorithm
along with its computational complexity.
For each of our test dataset with specified rate, we will run your algorithm
multiple times, s.t. you get a chance to sample a several of doping potentials,
and we will record any success.
