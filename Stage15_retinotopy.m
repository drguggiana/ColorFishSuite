%% Calculate retinotopy
%% Clean up

clearvars
close all
load('paths.mat')
addpath(genpath(paths(1).main_path))
%% Load transformation for retinotopy stack

registration_path = paths(1).registration_path;

%define the main loading path
tmat_path = fullfile(registration_path,'Registration_info','Anatomy','registration');

% %allocate memory for the affine objects
% aff_cell = cell(num_files,1);

% %go file by file extracting the matrix
% %allocate memory for the matrices
% x_mat = zeros(5,3,num_files);
% %for all the files
% for files = 1:num_files
%     %get the file name
%     [~,f_name,~] = fileparts(reformatted_cell{files});
%     [~,f_name] = fileparts(f_name);
% 
%     %get rid of the extension and add the "list" extension to define the
%     %target path
%     f_name = strcat(f_name,'.list');

    %load the actual matrix
%     temp_cell = importdata(fullfile(tmat_path,f_name,'registration'));
%     temp_cell = importdata(tmat_path);
%%
    fileID = fopen(tmat_path);
    temp_cell = textscan(fileID,'%s','TextType','char');
    fclose(fileID);
    %%
    
    decoded = jsondecode(reshape(char(temp_cell{1}),[],1));
    %and parse it
    %allocate memory for the parsed data
    parse_mat = zeros(5,3);
    %for all the relevant lines
    for plines = 6:10
        %split the string via spaces
        temp_str = strsplit(temp_cell{plines},' ');
        %convert the last 3 to numbers and store
        parse_mat(plines-5,:) = [str2double(temp_str{2});str2double(temp_str{3});str2double(temp_str{4})];
    end
    %store the parse matrix in the main matrix
    x_mat(:,:,files) = parse_mat;

% end   
%% Assign azimuth and elevation values to each voxel of the map

%% Load the actual data

%% Generate visual space maps for each stimulus

