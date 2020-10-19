%% Load obj file and translate regions to ref brain
clearvars
close all
load('paths.mat')
addpath(genpath(paths(1).main_path))
%% Get the path to the template stack


% get the target stack
target_path = uipickfiles('FilterSpec',fullfile(paths.registration_path,'Pre_registration_brains','*tif'));
target_path = target_path{1};

% define the zero threshold
zero_threshold = 200;
%% Load the template stack

% get the stack size
stack_info = imfinfo(target_path,'tif');
im_size = [stack_info(1).Height,stack_info(1).Width,size(stack_info,1)];

% allocate memory for the stack
template_stack = zeros(im_size);
% load the stack
% for all the frames
for frame = 1:im_size(3)
    template_stack(:,:,frame) = imread(target_path,frame);
end
%% Load the map images

% define the folder and file name
[folder,name] = fileparts(target_path);
folder = strsplit(folder,'\');
folder = folder{end};

% allocate memory for the maps
map_stacks = cell(3,1);

% define the dimension order
dim_order = [3,2,1];
% load the maps
% for all the maps
for maps = 1:3
    % define the map paths
    map_path = dir(fullfile(paths.registration_path,'Map_stacks',folder,'Reformatted_tif',...
        strcat('*0',num2str(dim_order(maps)+1),'*.tif')));
    map_path = fullfile(map_path.folder,map_path.name);
    % if it's the first map
    if maps == 1
        % load the stack info
        map_info = imfinfo(map_path);
        % define the image size
        im_size = [map_info(1).Height,map_info(1).Width,size(map_info,1)];
    end
    

    % allocate memory for the stack
    temp = zeros(im_size);
    % load the stack
    % for all the frames
    for frame = 1:im_size(3)
        temp(:,:,frame) = imread(map_path,frame);
    end
    % store in the cell
    map_stacks{maps} = temp;
end

% concatenate the stacks
map_stacks = cat(4,map_stacks{:});
% flip z
map_stacks(:,:,:,3)
%% Plot map stacks

close all

for dims = 1:3
    figure
    imagesc(max(map_stacks(:,:,:,dims),[],3))
end
%% Load the reference brain

% define the path
reference_path = fullfile(paths.registration_path,'Reference_brains\ZBB_isl2b-GFP_cut.tif');

% load the stack info
ref_info = imfinfo(reference_path);
% define the image size
ref_size = [ref_info(1).Height,ref_info(1).Width,size(ref_info,1)];


% allocate memory for the stack
ref_stack = zeros(ref_size);
% load the stack
% for all the frames
for frame = 1:ref_size(3)
    ref_stack(:,:,frame) = imread(reference_path,frame);
end
%% Load the full ref stack

% define the path
full_reference_path = fullfile(paths.registration_path,'Reference_brains\ZBB_isl2b-GFP.tif');

% load the stack info
fullref_info = imfinfo(full_reference_path);
% define the image size
fullref_size = [fullref_info(1).Height,fullref_info(1).Width,size(fullref_info,1)];


% allocate memory for the stack
fullref_stack = zeros(fullref_size);
% load the stack
% for all the frames
for frame = 1:fullref_size(3)
    fullref_stack(:,:,frame) = imread(full_reference_path,frame);
end
%% Invert the map to index directly

% allocate memory for the inverted map
inverted_map_ori = zeros([size(template_stack),3]);

% get the nonzero points
ref_idx = find(map_stacks(:,:,:,1));
% allocate memory for the template coordinates
ref_coord = zeros(size(ref_idx,1),3);
% for all the coordinates
for coord = 1:3
    % create a temp with ref dimensions
    temp = map_stacks(:,:,:,coord);
    % store the linear index of the values
    ref_coord(:,coord) = temp(ref_idx);
end

% get the linear indexes for the template
template_idx = sub2ind(size(template_stack),ref_coord(:,1),ref_coord(:,2),ref_coord(:,3));
% allocate memory for the ref coordinates in the template
template_coord = zeros(size(template_idx,1),3);
% get the subindexes for the ref
[template_coord(:,1),template_coord(:,2),template_coord(:,3)] = ...
    ind2sub(size(map_stacks(:,:,:,1)),ref_idx);

% for all the dimensions
for dims = 1:3
    % create a temp variable
    temp = zeros(size(template_stack));
    % load using the linear indexes
    temp(template_idx) = template_coord(:,dims);
    % load in the final map
    inverted_map_ori(:,:,:,dims) = temp;
end
%% OFF Plot stack as a scatter

% close all
% 
% ind = find(template_stack>200);
% [x,y,z] = ind2sub(size(template_stack),ind);
% 
% figure
% scatter3(x,y,z)
%% OFF Plot obj as scatter

% close all
% 
% 
% scatter3(obj_struct.v(1,:)',obj_struct.v(2,:)',obj_struct.v(3,:)')
% % scatter3(obj_struct.f3(1,:)',obj_struct.f3(2,:)',obj_struct.f3(3,:)')
%% Load the obj file

% define the path of the obj file
% obj_path = 'C:\Users\Drago\Downloads\Segment__0005_AF9.obj';
% obj_path = 'C:\Users\Drago\Downloads\Segment__0023_AF10.SO1-2.obj';
% obj_path = 'C:\Users\Drago\Downloads\Segment__0004_AF5.obj';

obj_path = fullfile(paths.registration_path,'AF_Segmentation_VAST','Segment__0023_AF10.SO1-2.obj');
tic
% read the file
obj_struct = loadawobj(obj_path);
toc
%% Convert obj coordinates to  stack coordinates

% define the template pixels sizes
x_size = 0.1758;
y_size = 0.1758;
z_size = 1;

% define the raw sizes (pre cropping and flips)
raw_size_x = 1064;
raw_size_y = 1064;
raw_size_z = 243;

% get the size of the stack in microns
total_x = raw_size_x*x_size;
total_y = raw_size_y*y_size;
total_z = raw_size_z*z_size;

% get the obj coordinates
x_obj = obj_struct.v(1,:)';
y_obj = obj_struct.v(2,:)';
z_obj = obj_struct.v(3,:)';

% rescale the object coordinates to pixels from um
x_rescale = round(x_obj.*raw_size_x./total_x);
y_rescale = round(y_obj.*raw_size_y./total_y);
z_rescale = round(z_obj.*raw_size_z./total_z);

% convert the coordinates according to the stack processing before
% registration (change between raw brain and preregistration brain). This
% includes cropping, flipping z and flipping x
x_rescale = size(template_stack,1)-(x_rescale-39);
y_rescale = y_rescale-39;
z_rescale = size(template_stack,3)+z_rescale;
% x_rescale = round(size(template_stack,2)-x_obj.*size(template_stack,1)./total_x);
% y_rescale = round(y_obj.*size(template_stack,2)./total_y);
% z_rescale = round(size(template_stack,3)+z_obj.*size(template_stack,3)./total_z);

% assemble a single matrix (flipping x and y due to the image/matrix
% coordinates)
rescaled_px = cat(2,y_rescale,x_rescale,z_rescale);
% exclude all points in the 0 and negative range (due to cropping for
% registration
rescaled_px = rescaled_px(sum(rescaled_px<1,2)==0,:);
%% OFF Plot an overlay per slice

% close all
% % define the target z
% for z_target = 100:10:240
% % z_target = 150;
% 
%     figure
%     imagesc(template_stack(:,:,z_target))
%     hold on
%     set(gca,'TickLength',[0 0])
%     scatter(rescaled_px(rescaled_px(:,3)==z_target,1),rescaled_px(rescaled_px(:,3)==z_target,2),0.1,'k')
%     pause(0.05)
% end
%% OFF Plot a collapsed overlay
% close all
% figure
% 
% imagesc((imrotate(max(template_stack,[],3),0)))
% hold on
% scatter(y_rescale,x_rescale,0.1,'k')
%% Interpolate the inverted maps
tic

% create a copy of the inverted map to edit
inverted_map = inverted_map_ori;
inverted_idx = find(inverted_map(:,:,:,1)>0);
[invX,invY,invZ] = ind2sub(size(template_stack),inverted_idx);


% zero_idx = find(inverted_map(:,:,:,1)==0);
% [zeroX,zeroY,zeroZ] = ind2sub(size(template_stack),zero_idx);

% find the missing points from the template
rescaled_ind = sub2ind(size(template_stack),rescaled_px(:,2),rescaled_px(:,1),rescaled_px(:,3));
rescaled_idx = inverted_map(:,:,:,1);
rescaled_idx = find(rescaled_idx(rescaled_ind)==0);
[zeroX,zeroY,zeroZ] = ind2sub(size(template_stack),rescaled_idx);

% for all dims
for dims = 1:3
    % get the nonzeros values
    temp = inverted_map(:,:,:,dims);
    inverted_vals = temp(inverted_idx);
    % create the interpolant
%     interpolated_points = griddatan([invX,invY,invZ],inverted_vals,[zeroX,zeroY,zeroZ]);
    interpolated_points = griddata(invX,invY,invZ,inverted_vals,zeroX,zeroY,zeroZ);
    
    
%     % query the interpolant
%     interpolated_points = interpolant(zeroX,zeroY,zeroZ);
    
    % save the values
    temp(rescaled_idx) = round(interpolated_points);
    % overwrite the maps
    inverted_map(:,:,:,dims) = temp;
    
end
toc
%% OFF Plot the inverted map

% close all
% 
% for dims = 1:3
%     figure
%     imagesc(max(inverted_map(:,:,:,dims),[],3))
% end
%% OFF Plot the inverted map in slices

% close all
% % define the target dim
% target_dim = 3;
% % define the target z
% for z_target = 100:10:240
% % z_target = 150;
% 
%     figure
%     imagesc(squeeze(inverted_map(:,:,z_target,target_dim)))
%     hold on
%     set(gca,'TickLength',[0 0])
%     colorbar
%     pause(0.05)
% end
%% Reformat the region coordinates

% allocate memory for the new points
reformat_px = zeros(size(rescaled_px));

% for all the points
for points = 1:size(reformat_px,1)

    
%     if rescaled_px(points,1)>size(inverted_map,1) || ...
%             rescaled_px(points,2)>size(inverted_map,2) || ...
%             rescaled_px(points,3)>size(inverted_map,3)
%         reformat_px(points,:) = [NaN,NaN,NaN];
%         continue
%     end
    reformat_px(points,1) = inverted_map(rescaled_px(points,2),rescaled_px(points,1),rescaled_px(points,3),1);
    reformat_px(points,2) = inverted_map(rescaled_px(points,2),rescaled_px(points,1),rescaled_px(points,3),2);
    reformat_px(points,3) = inverted_map(rescaled_px(points,2),rescaled_px(points,1),rescaled_px(points,3),3);
end
%% Plot the overlay on the reference brain

close all
figure

imagesc(max(ref_stack,[],3))
hold on
scatter(reformat_px(:,2),reformat_px(:,1),1,'k')
%% Plot per section

close all
% define the target z
for z_target = 1:5:125
% z_target = 150;

    figure
    imagesc(ref_stack(:,:,z_target))
    hold on
    set(gca,'TickLength',[0 0])
    scatter(reformat_px(reformat_px(:,3)==z_target,2),reformat_px(reformat_px(:,3)==z_target,1),1,'k')
    pause(0.05)
end
%% Correct the coordinates to fit the full stack

full_coordinates = reformat_px;
full_coordinates(:,2) = full_coordinates(:,2) + 288;
full_coordinates(:,3) = full_coordinates(:,3) + 10;
%% Create a mask from the reformatted area

% allocate memory for the mask
mask = zeros(fullref_stack);





