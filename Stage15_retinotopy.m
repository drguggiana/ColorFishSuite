%% Calculate retinotopy
%% Clean up

clearvars
close all
load('paths.mat')
addpath(genpath(paths(1).main_path))
%% Load the full ref brain

% define the path of the map stack
full_path = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\registration\Reference_brains\refbrain.nrrd.tif';

% load the stack

% get the stack info
full_info = imfinfo(full_path);

% allocate memory for the stack
full_stack = zeros(full_info(1).Height,full_info(1).Width,length(full_info));

% for all the frames
for frames = 1:length(full_info)
    full_stack(:,:,frames) = imread(full_path,frames);
end
%% Load the calcium data

load('E:\Behavioral data\Matlab\AF_proc\AFSuite\Analysis_files\Stage1\Retinotopy\pre\5432_rep2_z1_5_traces.mat')
%% OFF Save a stack to use as map for the reformatting (only run once)

% Centroids only
% % get the coordinates of the ROIs
% xy_coord = round(cat(1,seed_concat.centroid));
% z_coord = z_seed;
% % turn them into linear indexes
% linear_idx = sub2ind(size(ave_stack),xy_coord(:,2),xy_coord(:,1),z_coord);
% % generate the roi labels
% label_vector = 1:length(z_coord);

% Full ROIS
z_coord = z_seed;

% get the number of ROIs
roi_number = length(z_coord);
% get the xy full coordinates
[x_full,y_full] = ind2sub(size(ave_stack(:,:,1)),cat(1,seed_concat.pxlist));
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
    area = seed_concat(roi).area;
    % put as many copies of the z in those positions
    xyz_coord(counter:counter+area-1,3) = z_coord(roi);
    % fill the label vector
    label_vector(counter:counter+area-1) = label_counter;
    % update the counters
    counter = counter + area;
    label_counter = label_counter + 1;
end

% turn them into linear indexes
linear_idx = sub2ind(size(ave_stack),xyz_coord(:,1),xyz_coord(:,2),xyz_coord(:,3));


% allocate memory for the stack
map_stack = zeros(size(ave_stack));

% populate the stack
map_stack(linear_idx) = label_vector;
% turn it into uint16
map_stack = uint16(map_stack);

% define the saving path
tif_path = fullfile(paths.registration_path,'Pre_registration_brains','Retinotopy','5432_rep2_z1_5_anato_map.tif');
% save the stack
% for all the frames
for frames = 1:size(map_stack,3)
    if frames == 1
        oa = 'overwrite';
    else
        oa = 'append';
    end
    imwrite(map_stack(:,:,frames),tif_path,'tif','WriteMode',oa)
end
%% Load transformation for retinotopy stack

% define the path of the map stack
map_path = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\registration\Reformatted_brains\Retinotopy\5432_rep2_z1_5_anato_02_9dof.tif';

% load the stack

% get the stack info
map_info = imfinfo(map_path);

% allocate memory for the stack
map_stack = zeros(map_info(1).Height,map_info(1).Width,length(map_info));

% for all the frames
for frames = 1:length(map_info)
    map_stack(:,:,frames) = imread(map_path,frames);
end

disp(length(unique(map_stack(:))))
%% Get the ROI coordinates 

% get the remaining ROI IDs
remaining_ROIs = unique(map_stack(:));
% get rid of the 0
remaining_ROIs = remaining_ROIs(2:end);
% get the number of remaining ROIs
rem_number = length(remaining_ROIs);

% find the linear coordinates
linear_coord = find(map_stack);

% get the ROI values
ROI_values = map_stack(linear_coord);

% get the 3d coordinates
[x_reg,y_reg,z_reg] = ind2sub(size(map_stack),linear_coord);
% correct them based on the map used for registration
x_reg = x_reg + 325;
z_reg = z_reg + 29;

% get the linear coordinates for the full stack
linear_full = sub2ind(size(full_stack),x_reg,y_reg,z_reg);

% allocate memory for the full mask
retinotopy_mask = zeros(size(full_stack));

% position the ROIs in the full stack
retinotopy_mask(linear_full) = ROI_values;

% sort the indexes according to the ROI
[ROI_values,sort_idx] = sort(ROI_values);
linear_full = linear_full(sort_idx);

% get the area of each ROI
roi_area = histcounts(ROI_values,[unique(ROI_values);max(ROI_values)]);
%% OFF Visualize the mask

% close all
% % define the target z
% for z_target = 1:5:138
% % z_target = 150;
% 
%     figure
% %     imagesc(retinotopy_mask(:,:,z_target))
%     C = imfuse(retinotopy_mask(:,:,z_target),full_stack(:,:,z_target));
%     image(C)
%     hold on
%     set(gca,'TickLength',[0 0])
%     axis equal
%     pause(0.05)
% end
% 
% figure
% % imagesc(max(retinotopy_mask,[],3))
% C = imfuse(max(retinotopy_mask,[],3),max(full_stack,[],3));
% image(C)
% axis equal
%% Assign azimuth and elevation values to each ROI of the map

%0'grey OFF ON grey
%1'back right';
%2'back');
%3'backleft');
%4'LEFT');
%5'left fwd');
%6'FWD');
%7 RIGHT fwd)');
%8'RIGHT');
%9'grey ON OFF grey
%10'back right';
%11'back');
%12'backleft');
%13'LEFT');
%14'left fwd');
%15'FWD');
%16 RIGHT fwd)');
%17'RIGHT');
%18'back right';
%19'back');
%20'backleft');
%21'LEFT');
%22'left fwd');
%23'FWD');
%24 RIGHT fwd)');
%25'RIGHT');


% define the stimulation window
% time_vector = 21:63;
% time_vector = 31:53;
time_vector = 21:60;

% get the number of timepoints
t = length(time_vector);

% interpolate the diagonals
pos_diagonal = (1:0.7071:0.7071*(t+1)) - 0.7071*20.5 + 20.5;
neg_diagonal = fliplr(pos_diagonal);

% define the mapping of space to time
spacetime = cell(8,2);

spacetime{1,1} = [pos_diagonal;pos_diagonal]';
spacetime{2,1} = [nan(1,t);1:t]';
spacetime{3,1} = [neg_diagonal;pos_diagonal]';
spacetime{4,1} = [t:-1:1;nan(1,t)]';
spacetime{5,1} = [neg_diagonal;neg_diagonal]';
spacetime{6,1} = [nan(1,t);t:-1:1]';
spacetime{7,1} = [pos_diagonal;neg_diagonal]';
spacetime{8,1} = [1:t;nan(1,t)]';


spacetime{1,2} = 'back_right';
spacetime{2,2} = 'back';
spacetime{3,2} = 'back_left';
spacetime{4,2} = 'left';
spacetime{5,2} = 'left_forward';
spacetime{6,2} = 'forward';
spacetime{7,2} = 'forward_right';
spacetime{8,2} = 'right';

% define the stimuli to use
stim_vector = 19:26; %BE

% get only the BE responses
responses = mean(reps_trace(remaining_ROIs,time_vector,stim_vector,:),4);

% allocate memory for the maxima
max_response = zeros(size(responses,1),size(responses,3),2);
% for all the stimuli
for stim = 1:size(responses,3)
    % get the max coordinate for that interval
%     [max_val,ind] = max(abs(diff(responses(:,:,stim),1,2)),[],2);
    [max_val,ind] = max(abs(diff(movmean(responses(:,:,stim),5,2),1,2)),[],2);
    ind = ind + 1;
    % get a selection vector to exclude weak ROIs
%     responsive_vector = double(max_val>3*std(responses(:,:,stim),0,2));
%     responsive_vector(responsive_vector==0) = NaN;
    responsive_vector = ones(size(responses,1),1);
    
    % save the corresponding time value
    max_response(:,stim,:) = responsive_vector.*spacetime{stim-floor((stim-1)/8).*8,1}(ind,:);
end
%% Visualize the partial maps
close all
% average the azimuth values
for maps = 1:2
    switch maps
        case 1
%             azimuth_percell = squeeze(nanmean(max_response(:,6,2),2));
            azimuth_percell = squeeze(nanmean(max_response(:,4,1),2));
        case 2
%             azimuth_percell = squeeze(nanmean(max_response(:,8,1),2));
            azimuth_percell = squeeze(nanmean(max_response(:,8,1),2));
    end

    % allocate a stack
    azimuth_map = zeros(size(full_stack));

    % allocate memory for the azimuth values
    azimuth_vector = zeros(length(roi_area),1);
    % initialize a counter for the idx
    counter = 1;
    % % initialize a counter for the labels
    % label_counter = 1;
    % for all the rois
    for roi = 1:rem_number
        % get the area of the roi
        area = roi_area(roi);
        % put as many copies of the azimuth in those positions
        azimuth_vector(counter:counter+area-1) = azimuth_percell(roi);
    %     % fill the label vector
    %     label_vector(counter:counter+area-1) = label_counter;
        % update the counters
        counter = counter + area;
    %     label_counter = label_counter + 1;
    end

    % apply the values to the map
    azimuth_map(linear_full) = azimuth_vector;
    %% Generate an interpolated map
    
%     % prepare the grid
%     [X,Y,Z] = meshgrid(400:1:800,1:1:size(full_stack,2),75:3:105);
%     pts = [X(:),Y(:),Z(:)];
%     % also get the coord in linear
%     linear_grid = sub2ind(size(full_stack),pts(:,1),pts(:,2),pts(:,3));
%     
%     % get the coordinates of the points in 3d
%     [x,y,z] = ind2sub(size(full_stack),linear_full);
%     
%     % get a vector to exclude nans
%     in_vector = ~isnan(azimuth_vector);
%     interpolant = scatteredInterpolant(x(in_vector),y(in_vector),z(in_vector),azimuth_vector(in_vector),...
%         'linear','none');
%     
%     vq = interpolant(pts(:,1),pts(:,2),pts(:,3));
%     % allocate a new map
%     new_map = zeros(size(azimuth_map));
%     % load the new values
%     new_map(linear_grid) = vq;
    %% Visualize the map
    
    new_map = azimuth_map;
    new_map(new_map==0) = NaN;
    
    % define interval
    interval = 10;
    figure
    z_vector = 75:interval:105;
    z_num = length(z_vector);
    count = 1;
    cmap = jet(256);
    cmap(1,:) = [0 0 0];
    % define the target z
    for z_target = z_vector
    % z_target = 150;
    
%         figure
        subplot(round(sqrt(z_num)),ceil(sqrt(z_num)),count)
        count = count + 1;
        imagesc(imgaussfilt(nanmean(new_map(400:800,:,z_target:z_target+interval-1),3)))
    %     C = imfuse(retinotopy_mask(:,:,z_target),full_stack(:,:,z_target));
    %     image(C)
        hold on
        set(gca,'TickLength',[0 0])
        axis equal
        colormap(cmap)
        set(gca,'visible','off')    
%         pause(0.05)
    end
    figure
    imagesc(nanmean(new_map(400:800,:,:),3))
    
    colormap(cmap)
    % C = imfuse(max(retinotopy_mask,[],3),max(full_stack,[],3));
    % image(C)
    axis equal
end
%% Interpolate for the whole map


%% Generate visual space maps for each stimulus

