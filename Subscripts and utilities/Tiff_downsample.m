%% Downsample the selected tiff files

%% clean up
clearvars
close all
load('paths.mat')
addpath(genpath(paths(1).main_path))
%% Define miscellaneous constants

% define the dowmsapling factor
downsampling_factor = 3;
%% Load the files and define paths

%get the folder where the image files are
tar_path_all = uipickfiles('FilterSpec','J:\2P');


%get the number of experiments selected
num_exp = length(tar_path_all);

%for all the experiments
for exp_ind = 1:num_exp
    
    %get the path of the current experiment
    tar_path = tar_path_all{exp_ind};

    %define the path for saving the file with the BIC info
    [~,tar_name] = fileparts(tar_path);

    %determine whether it's an RF file
    str_parts = strsplit(tar_name,'_');
    if strcmp(str_parts{3},'p18')||strcmp(str_parts{3},'p18proto')
        rf_var = 1;
        %define the amount of mode taking when aligning
        mode_wind = 4;
    else
        rf_var = 0;
        %define the amount of mode taking when aligning
        mode_wind = 5;
    end
    %% And parse them based on their filenames
    [file_info,stim_num,rep_num,z_num,tar_files] = parser(tar_path);

    % assemble the new folder name
    path_parts = strsplit(tar_path,'_');
    output_folder = strjoin({path_parts{1},path_parts{2},strcat(path_parts{3},'downsample'),path_parts{4}},'_');
    % make it if it doesn't exist yet
    if ~isfolder(output_folder)
        status = mkdir(output_folder);
    end
    
    % for all the files
    for files = 1:length(tar_files)
        % get the current file path
        current_path = fullfile(tar_path,tar_files{files});
        % also the current output path
        current_output = fullfile(output_folder,tar_files{files});
        % load the info only in the first file
        if  files == 1
            % get the stack info
            stack_info = imfinfo(current_path);
        end
        % for all the frames
        for frames = 1:size(stack_info)
            % load the frame
            frame = imread(current_path,frames);
            % downsample
            ds_frame = imresize(frame,1/downsampling_factor);
            % save the frame
            if frames == 1
                imwrite(ds_frame,current_output,'tif','WriteMode','overwrite')
            else
                imwrite(ds_frame,current_output,'tif','WriteMode','append')
            end
        end
    end  
    % copy the TDMS file
    status2 = copyfile(strcat(tar_path,'.TDMS'),strcat(output_folder,'.TDMS'));
    %% Also downsample the anatomy
    
    
end