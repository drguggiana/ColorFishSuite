%% Save tif file from n hdf5 data file

%% Clean up

close all
clearvars
load('paths.mat')
addpath(genpath(paths(1).main_path))
%% Define paths and load dataset

% define the output folder
output_folder = fullfile(paths.registration_path,'Reference_brains');

% get the path of the hdf5 file
% target_path = fullfile(paths.registration_path,'ZBrain_datasets','AnatomyLabelDatabase.hdf5');
target_path = fullfile(paths.registration_path,'ZBrain_datasets','ZBBDatabase.hdf5');

% display the contents
h5disp(target_path)

% define the target dataset
% target_dataset = 'Isl2bGal4-uasDendra_6dpf_MeanImageOf8Fish';
target_dataset = 'ZBB_isl2b-GFP';

% load the dataset
target_data = h5read(target_path,strcat('/',target_dataset));
%% Save the file as tif

% define the output path
output_path = fullfile(output_folder,strcat(target_dataset,'.tif'));

% get the number of frames
z_number = size(target_data,3);

% for all the frames
for z = 1:z_number
    % select depending on first frame
    if z == 1
        imwrite(target_data(:,:,z),output_path,'tif','WriteMode','overwrite')
    else
        imwrite(target_data(:,:,z),output_path,'tif','WriteMode','append')
    end
end