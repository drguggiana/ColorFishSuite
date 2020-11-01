%% Add AF info from registration

%% Clean up and load the structure

clearvars
close all
load('paths.mat')
addpath(genpath(paths(1).main_path))
registration_path = paths(1).registration_path;

cluster_path = paths(1).stage3_path;
% fig_path = strcat(paths(1).fig_path,'Classify\');

data = load_clusters(cluster_path);
%% Load the registration mask

% define the mask path
mask_path = fullfile(paths.registration_path,...
    'AF_Segmentation_VAST\full_mask\AF_mask.tif');

% get the stack info
mask_info = imfinfo(mask_path);

% allocate memory for the stack
mask_stack = zeros(mask_info(1).Height,mask_info(1).Width,size(mask_info,1));

% for all the frames
for frames = 1:size(mask_info,1)
    mask_stack(:,:,frames) = imread(mask_path,frames);
end
%% Load the labels

% define the path to the file
labels_file = fullfile(registration_path,'Labels_info','MaskDatabase.mat');
% load the labels file
labels_data = load(labels_file);

% define which dataset to load depending on the reference
if contains(data.name,{'Syn','syn'})
    field_list = {'AF8','AF4','AF5','AF9','AF6','AF7','SFGS1_2','SFGS3_4','SFGS5_6','SO','SGC','SAC_SPV'};
%     field_list = {'AF4','AF5','AF6','AF7','AF8','AF9','Tecum Neuropil','Periventriculare'};
%     coordinate_conversion = [274,0,0];
    
    label_stack = mask_stack;
    
else
    field_list = {'Tectum Stratum','Tecum Neuropil','Pretectum','Habenula','Cerebellum'};
    coordinate_conversion = [272,74,49];


    % get the number of fields
    field_number = length(field_list);
    % allocate memory for the fields
    label_stack = zeros(labels_data.height,labels_data.width,labels_data.Zs);
    % get the field numbers for these names
    for field = 1:length(field_list)
        % get the index with the first name match (so as to not get
        % subdivisions)
        idx_vector = find(contains(labels_data.MaskDatabaseNames,field_list{field}),1);
        % store the map and the name in the cell
        label_stack = label_stack + reshape(full(labels_data.MaskDatabase(:,idx_vector)),...
            labels_data.height,labels_data.width,labels_data.Zs).*field;

    end
end
%% Correct for fish of origin (p17b_syngc6s in particular)
close all

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
%% Get the anatmmy for each ROI

% get the coordinates of the ROIs
xyz_roi = data.registered_coord;

% turn them to linear
linear_roi = sub2ind(size(mask_stack),xyz_roi(:,1),xyz_roi(:,2),xyz_roi(:,3));
% exclude the NaNs
not_nan = ~isnan(linear_roi);
linear_roi = linear_roi(not_nan);

% allocate memory for the new anatomy
new_anatomy = zeros(size(xyz_roi,1),2);

% load the anatomy
new_anatomy(not_nan,1) = mask_stack(linear_roi);
% load the z of origin
new_anatomy(:,2) = xyz_roi(:,3);
%% Visualize the overlap

% create a temp stack with the ROIs
roi_stack = zeros(size(mask_stack));
roi_stack(linear_roi) = 1;

close all
% define the target z
for z_target = 1:5:138
% z_target = 150;

    figure
%     imagesc(retinotopy_mask(:,:,z_target))
    C = imfuse(roi_stack(:,:,z_target),mask_stack(:,:,z_target));
    image(C)
    hold on
    set(gca,'TickLength',[0 0])
    axis equal
    pause(0.05)
end

figure
% imagesc(max(retinotopy_mask,[],3))
C = imfuse(max(roi_stack,[],3),max(mask_stack,[],3));
image(C)
axis equal
%% Reformat the data

% define the ref brain depending on the file
if contains(data(1).name,{'Syn','syn'})
    ref_name = 'refisl2cut2.nrrd.tif';
%     full_ref_name = 'refisl2.nrrd.tif';
    full_ref_name = 'refbrain.nrrd.tif';

%     ref_name = 'refcutblursub.nrrd.tif';
%     full_ref_name = 'refbrain.nrrd.tif';
else
    ref_name = 'refcutblursub.nrrd.tif';
    full_ref_name = 'refbrain.nrrd.tif';
end
% assemble the reference path
ref_path = fullfile(registration_path,'Reference_brains',ref_name);
% load the ref info
ref_info = imfinfo(ref_path);

% reformat
fishave_cell = reformat_fish(data,fish_ori,num_files,im_info_cell,aff_cell,ref_info,'seeds');
%% Determine the region for each seed from the registration

% allocate memory for the registred anatomy
registered_anatomy = cell(num_files,1);
label_copy = zeros(size(label_stack));
corrected_coordinates = cell(num_files,1);
% for all the experiments
for files = 1:num_files
    % load the coordinates
    current_coord = fishave_cell{files,1};
    % correct the coordinates for the different references
    current_coord = round(current_coord + coordinate_conversion);
    
    % record the region for each seed
    % allocate memory for the seeds
    seed_region = zeros(size(current_coord,1),1);
    % for all the seeds
    for seeds = 1:size(current_coord,1)  

        if current_coord(seeds,1) > size(label_stack,1) || ...
                current_coord(seeds,2) > size(label_stack,2) || ...
                current_coord(seeds,3) > size(label_stack,3)
            continue
        else
%             label_copy(current_coord(seeds,2),...
%                     current_coord(seeds,1),current_coord(seeds,3)) = 1;
            seed_region(seeds) = label_stack(current_coord(seeds,1),...
                current_coord(seeds,2),current_coord(seeds,3));
        end
    end
    % store the result
    registered_anatomy{files} = seed_region;
    % also store the coordinates
    corrected_coordinates{files} = current_coord;
end

% turn the anatomy into a vector
registered_anatomy = cat(1,registered_anatomy{:});
% same with the coordinates
corrected_coordinates = cat(1,corrected_coordinates{:});

% NaN points outside the stack
nan_vector = corrected_coordinates(:,1)>size(mask_stack,1)|...
    corrected_coordinates(:,2)>size(mask_stack,2)|...
    corrected_coordinates(:,3)>size(mask_stack,3);

corrected_coordinates(nan_vector,:) = NaN;
registered_anatomy(nan_vector) = NaN;

% figure
% imagesc(cat(3,sum(label_stack,3)>0,sum(label_copy,3)>0,sum(label_stack,3)>0))
% axis equal
%% Train a kNN classifier for the non-allocated ROIs

% define the number of folds
fold_number = 5;
% get the training data and labels
train_labels = registered_anatomy(registered_anatomy>0);
train_data = corrected_coordinates(registered_anatomy>0,:);

classifier = fitcecoc(train_data,train_labels,'KFold',fold_number,'Verbose',0,...
            'Learners','kNN','Coding','onevsall');
%% Predict the assignment of the non-allocated ROIs

% get the coordinate data for those ROIs
non_allocated = corrected_coordinates(registered_anatomy==0,:);

% allocate memory for the predictions
predicted_AFs = zeros(size(non_allocated,1),fold_number);

% predict
% for all the folds
for folds = 1:fold_number
    predicted_AFs(:,folds) = predict(classifier.Trained{folds},non_allocated);
end

% take the mode across folds
prediction_mode = mode(predicted_AFs,2);
%% Assemble the complete vector

% copy the original annotation
new_anatomy = registered_anatomy;
% insert the predictions
new_anatomy(registered_anatomy==0) = prediction_mode;
%% Create a visualization stack of the allocation

% generate an empty volume
vis_stack = zeros(size(mask_stack));

% turn the seed coordinates to linear
linear_idx = sub2ind(size(mask_stack),corrected_coordinates(:,1),...
    corrected_coordinates(:,2),corrected_coordinates(:,3));
% get the unique regions
region_list = unique(new_anatomy);
% remove NaNs
region_list = region_list(~isnan(region_list));
% get the number of regions
region_number = size(region_list,1);
% populate with all the seeds
% for all the regions
for region = 1:region_number
    % get the logical vector for the region
    region_vector = new_anatomy==region_list(region);
    vis_stack(linear_idx(region_vector)) = new_anatomy(region_vector);
end
%% Plot the visualization

% generate the colormap
cmap = jet(region_number);
close all
% define the target z
for z_target = 1:5:size(mask_stack,3)
% z_target = 150;

    figure
    image(vis_stack(:,:,z_target))
    colormap(cmap)
    axis equal
    hold on
    set(gca,'TickLength',[0 0])
    pause(0.05)
end

figure
image(max(vis_stack,[],3))
colormap(cmap)
axis equal
