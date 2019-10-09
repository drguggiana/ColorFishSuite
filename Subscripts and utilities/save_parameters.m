% create an analysis parameters variable to load from the different scripts

% define the path to save things to
save_path = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Subscripts and utilities\parameters.mat';

% define the stage 2 parameters
snr_param = struct([]);

% define a field counter
field_c = 1;

snr_param(field_c).strain = 'gc6s';
snr_param(field_c).program = 'p17b';
snr_param(field_c).percentile = 50;
snr_param(field_c).stimThreshold = 0;
field_c = field_c + 1;

snr_param(field_c).strain = 'h2b6s';
snr_param(field_c).program = 'p17b';
snr_param(field_c).percentile = 50;
snr_param(field_c).stimThreshold = 0;
field_c = field_c + 1;

snr_param(field_c).strain = 'syngc6s';
snr_param(field_c).program = 'p17b';
snr_param(field_c).percentile = 50;
snr_param(field_c).stimThreshold = 0;
field_c = field_c + 1;

snr_param(field_c).strain = 'gc6s';
snr_param(field_c).program = 'p6';
snr_param(field_c).percentile = 50;
snr_param(field_c).stimThreshold = 8;
field_c = field_c + 1;

snr_param(field_c).strain = 'SynG6s';
snr_param(field_c).program = 'p8';
snr_param(field_c).percentile = 50;
snr_param(field_c).stimThreshold = 1;
field_c = field_c + 1;

% save the variable
save(save_path,'snr_param')