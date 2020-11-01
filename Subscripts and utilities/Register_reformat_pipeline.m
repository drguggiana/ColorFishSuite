%% Run registration and reformat
%% Clean up
clearvars
close all
load('paths.mat')
addpath(genpath(paths(1).main_path))
run(fullfile(paths.registration_path,'Munger_Fiji_scripts','NRRD_properties.m'));
%% Load the target files

% select the target folder containing the pre-registration tifs
% target_folder = uipickfiles('FilterSpec',fullfile(paths.registration_path,'Pre_registration_brains'));
% target_folder = target_folder{1};
% target_folder = fullfile(paths.registration_path,'Pre_registration_brains','p17b_syngc6s');
% get the dataset name
data = load_clusters(paths.stage3_path);

% split the name
name_parts = strsplit(data.name,'_');
% get a list of the files
file_list = dir(fullfile(paths.stage1_path,strcat('*',name_parts{2},'*',name_parts{1},'*anato.tif')));
% % get the info for the first one
% file_info = imfinfo(fullfile(file_list(1).folder,file_list(1).name));

% define the path to save the NRRD files to
nrrd_path = fullfile(paths.registration_path,'Pre_registration_brains',data.name);
%% Correct for fish of origin (p17b_syngc6s in particular)

% get the fish of origin
fish_ori = data.fish_ori;
% correct the data in case there were multiple experiments per fish
% get the number of experiments
num_experiment = sum(diff(fish_ori(:,2))<0)+1;
num_fish = length(unique(fish_ori(:,1)));
% if the numbers are different, correct
if num_fish < num_experiment
    % allocate memory for the new fish_ori
    fish_ori_corrected = zeros(size(fish_ori));
    fish_ori_corrected(:,2) = fish_ori(:,2);
    % find the junctions
    junctions = [0;find(diff(fish_ori(:,2))<0);size(fish_ori,1)];

    % for all the junctions
    for junc = 1:length(junctions)-1
        fish_ori_corrected(junctions(junc)+1:junctions(junc+1),1) = junc;
    end
    % replace the old fish_ori with the new one for registration only
    fish_ori = fish_ori_corrected;
end
%% Loop through the files

% allocate memory for the new coordinates
new_coordinates = cell(length(file_list),1);

for files = 1:length(file_list)
    %% Read the ROI coordinates
    
    % get the vector for this fish
    fish_vector = fish_ori(:,1)==files;
    % get the coordinates
    xy_info = data.xy_seed(fish_vector);
%     xy_info = round(cat(1,xy_info.centroid));
    z_info = data.z_seed(fish_vector);
    %% Generate a stack with the mapped ROI locations and save as tif
    % get the file info
    file_info = imfinfo(fullfile(file_list(files).folder,file_list(files).name));
    
    % allocate memory for the stack
    roi_stack = zeros(file_info(1).Height,file_info(1).Width,size(file_info,1));
    
    % force z_info to start on 1
    if files > 1
        z_info = z_info - z_info(1) + 1;
    end
    
%     % linearize the indexes
%     linear_idx = sub2ind(size(roi_stack),xy_info(:,1),xy_info(:,2),z_info);
    
    % Full ROIS
%     z_coord = z_seed;
    
    % get the number of ROIs
    roi_number = length(z_info);
    % get the xy full coordinates
    [x_full,y_full] = ind2sub(size(data.ave_stack(:,:,1)),cat(1,xy_info.pxlist));
    % allocate memory for the full coordinates
    xyz_coord = zeros(length(x_full),3);
    xyz_coord(:,1) = x_full;
    xyz_coord(:,2) = y_full;
    % allocate memory for the labels
    label_vector = zeros(length(x_full),1);
    % initialize a counter for the idx
    counter = 1;
    % initialize a counter for the labels
    label_counter = 1;
    % for all the rois
    for roi = 1:roi_number
        % get the area of the roi
        area = xy_info(roi).area;
        % put as many copies of the z in those positions
        xyz_coord(counter:counter+area-1,3) = z_info(roi);
        % fill the label vector
        label_vector(counter:counter+area-1) = label_counter;
        % update the counters
        counter = counter + area;
        label_counter = label_counter + 1;
    end
    
    % turn them into linear indexes
    linear_idx = sub2ind(size(roi_stack),xyz_coord(:,1),xyz_coord(:,2),xyz_coord(:,3));
    
    % get the number of rois
    roi_number = size(z_info,1);
    
    % fill up the stack
    roi_stack(linear_idx) = label_vector;
    % convert to int
    roi_stack = uint16(roi_stack);
    
    % get the name of the fish
    tif_name = file_list(files).name;
    % save the stack
    save_tif_path = fullfile(paths.registration_path,'ROI_maps',data.name);
    % for all the frames
    for frames = 1:size(file_info,1)
        if frames == 1
            oa = 'overwrite';
        else
            oa = 'append';
        end
        imwrite(roi_stack(:,:,frames),fullfile(save_tif_path,tif_name),'tif','WriteMode',oa)
    end
    %% Generate the nrrd files
    
    % define the source name
    source_name = file_list(files).name;
    % get the properties
    property = properties(contains({properties.name},file_list(files).name(1:end-4)));
    nrrd_conf = property.string;
    % for both the raw and map
    for stack = 1:2
        
        % select the source path file
        switch stack
            case 1
                source_path = file_list(files).folder;
            case 2
                source_path = save_tif_path;
        end
        

        % define the digit for registration
        digit = strcat('0',num2str(stack));

        
        % assemble the arguments
        arguments = strjoin({source_path,source_name,nrrd_path,nrrd_conf,digit},',');
        % flip all the slashes
        arguments = strrep(arguments,'\','/');
        
        % assemble the path to the macro
        macro_path = fullfile(paths.registration_path,'Munger_Fiji_scripts','tif_2_nrrd.ijm');
        % assemble the command string
        command_string = sprintf('%s -batch "%s" "%s"',paths.imagej_path,macro_path,arguments);
        % run Fiji
        status = system(command_string);
    end
    %% Run the registration
    
    % load the settings for this dataset
    reference = property.reference;
    settings = property.registration;
    target_path = fullfile(nrrd_path,strrep(source_name,'.tif','_01.nrrd'));
    reformat_path = fullfile(nrrd_path,strrep(source_name,'.tif','_02.nrrd'));
    % define the target directory
    reformatted_path = strrep(fullfile('Reformatted_brains',...
        data.name,strcat(source_name(1:end-4),'_reformatted_DIGIT.nrrd')),'\','/');
    
    % assemble the command
    cygwin_command = strcat(paths.cygwin_path,' --login -c',...
        ' "',paths.munger_path,' -3',' ''',reformatted_path,'''',settings,reference,' ''',...
        target_path,''' ''',reformat_path,''' "');
    
%     cygwin_command = strcat('C:\cygwin64\bin\bash.exe --login -c',' "C:/cfolder/usr/local/bin/munger2',' -a -w -r [01,02] -l f -v -T 8 -X 52 -C 4 -G 15 -R 1 -A', ' ''--accuracy 0.4''',' -s ''E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\registration\Reference_brains\ZBB_isl2b-GFP_cut.nrrd'' ''E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\registration\AF_Segmentation_VAST\preregistration_nrrd\7dpf_top_down_nocrop_01.nrrd'' ''E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\registration\AF_Segmentation_VAST\preregistration_nrrd\7dpf_top_down_nocrop_02.nrrd'' "');
    % run the registration
    system(cygwin_command)
    %% Retrieve the reformatted map back to tif
    
    % build the path to the file (there are forbidden characters, hence the
    % tossing and turning)
    reformat_dir = dir(fullfile('C:\cygwin64\home\Drago','Reformatted_brains',data.name,'*02.nrrd'));
    % find the current file
    file_vector = contains({reformat_dir.name},source_name(1:end-4));
    % get the name
    file_name = reformat_dir(file_vector).name;
    % define the target directory
    reformatted_directory = fullfile(paths.registration_path,'Reformatted_brains',data.name);
    
    % assemble the arguments
    arguments = strjoin({reformat_dir(1).folder,file_name,reformatted_directory,},',');
    % flip all the slashes
    arguments = strrep(arguments,'\','/');
    
    % assemble the path to the macro
    macro_path = fullfile(paths.registration_path,'Munger_Fiji_scripts','nrrd_2_tif.ijm');
    % assemble the command string
    command_string = sprintf('%s -batch "%s" "%s"',paths.imagej_path,macro_path,arguments);
    % run Fiji
    status = system(command_string);
    %% Load the tif
    
    % define the path
    tif_path = fullfile(reformatted_directory,strrep(file_name,'.nrrd','.tif'));
    % get the info
    reformat_info = imfinfo(tif_path);
    
    % allocate memory for the stack
    reformatted_stack = zeros(reformat_info(1).Height,reformat_info(1).Width,length(reformat_info));
    
    % for all the frames
    for frames = 1:length(reformat_info)
        % load the stack
        reformatted_stack(:,:,frames) = imread(tif_path,frames);
    end
    
    disp(length(unique(reformatted_stack(:))))
    %%  Save the coordinates in a joint cell
    
    % get a list of the rois that were recovered
    recovered_rois = unique(reformatted_stack(:));
    % get the coordinates of all the recovered ones
    recovered_xyz = find(reformatted_stack);
    % get the values
    recovered_values = reformatted_stack(recovered_xyz);
    % turn the coordinates to 3d
    [rec_x,rec_y,rec_z] = ind2sub(size(reformatted_stack),recovered_xyz);
    recovered_xyz = [rec_x,rec_y,rec_z];
    % allocate memory for the centroid coordinates
    xyz_new = zeros(roi_number,3);
    % for all the original ROIs
    for roi = 1:roi_number
        % if not recovered, save a NaN and move on
        if sum(recovered_rois==roi)==0
            xyz_new(roi,:) = NaN;
            continue
        end
        % get the pixels
        px_coord = recovered_xyz(recovered_rois==roi,:);
        % save the z
        xyz_new(roi,3) = px_coord(1,3);
        % find the centroid of the coordinates and save
        xyz_new(roi,1:2) = round(mean(px_coord(:,1:2),1));

    end
    % correct to the full ref stack
    xyz_new = xyz_new + property.ref_correction;
    % store in the main cell
    new_coordinates{files} = xyz_new;
end

% concatenate the cell
new_coordinates = cat(1,new_coordinates{:});
%% write in the structure and save

% copy to preserve backwards compatibility
main_str = data;

% add the info to the structure
main_str.registered_coord = new_coordinates;

% build the path
structure_path = fullfile(paths.stage3_path,strcat(data.name,'_clusters.mat'));

%save
save(structure_path,'main_str')