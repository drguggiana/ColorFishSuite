%% Load obj file, turn into reformat-ready stack and save
%% Clean up
clearvars
close all
load('paths.mat')

addpath(genpath(paths(1).main_path))

tic
%% Load the template stack

% get the target stack path
target_path = fullfile(paths.registration_path,'Pre_registration_brains\Anatomy','7dpf_top_down_nocrop.tif');

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
%% Get the obj paths and loop through

% get the paths
obj_list = dir(fullfile(paths.registration_path,'AF_Segmentation_VAST','*.obj'));

% get the number of files
obj_num = size(obj_list,1);

% for all the files
for files = 1:obj_num
    %% Load the obj file
    
    % assemble the obj path
    obj_path = fullfile(obj_list(files).folder,obj_list(files).name);
    % define the path of the obj file
    % obj_path = 'C:\Users\Drago\Downloads\Segment__0005_AF9.obj';
    % obj_path = 'C:\Users\Drago\Downloads\Segment__0023_AF10.SO1-2.obj';
    % obj_path = 'C:\Users\Drago\Downloads\Segment__0004_AF5.obj';
    
    % obj_path = fullfile(paths.registration_path,'AF_Segmentation_VAST','Segment__0023_AF10.SO1-2.obj');
    % read the file
    obj_struct = loadawobj(obj_path);
    %% Convert obj coordinates to  stack coordinates
    
    % define the template pixels sizes
    x_size = 0.1758;
    y_size = 0.1758;
    z_size = 1;
    
    % define the raw sizes (pre cropping and flips)
%     raw_size_x = 1064;
%     raw_size_y = 1064;
%     raw_size_z = 243;
    raw_size_x = size(template_stack,1);
    raw_size_y = size(template_stack,2);
    raw_size_z = size(template_stack,3);
    
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
    x_rescale = raw_size_x-x_rescale;
    z_rescale = raw_size_z+z_rescale;
    
    % assemble a single matrix (flipping x and y due to the image/matrix
    % coordinates)
    rescaled_px = cat(2,y_rescale,x_rescale,z_rescale);
    % exclude all points in the 0 and negative range (due to cropping for
    % registration
    rescaled_px = rescaled_px(sum(rescaled_px<1,2)==0,:);
    %% Generate a stack with the outline
    
    % generate the stack
    output_stack = zeros(size(template_stack));
    
    % convert the coordinates to linear
    rescaled_idx = sub2ind(size(template_stack),rescaled_px(:,2),rescaled_px(:,1),rescaled_px(:,3));
    
    % populate it
    output_stack(rescaled_idx) = 1;
    %% Fill in the outline
    
    % for all the frames
    for frames = 1:size(output_stack,3)
        % if it's empty, skip
        if sum(output_stack(:,:,frames),'all') == 0
            continue
        end
        % fill in the hole
        output_stack(:,:,frames) = imfill(imclose(output_stack(:,:,frames),[1 1 1;1 1 1;1 1 1]),4,'holes');
    end
    
    % convert to 16 bit
    output_stack = uint16(output_stack);
    %% Combine with the mixed stack
    
    % if it's the first iteration, create it
    if files == 1
        mixed_stack = output_stack;
    else
        % zero the overlaps
        and_matrix = mixed_stack&output_stack;
        disp(num2str(sum(and_matrix(:))))
        output_stack(and_matrix) = 0;
        % combine the matrices
        mixed_stack = mixed_stack + output_stack.*files;
    end
    
end
    %% OFF Plot the overlap
    
%     close all
%     % define the target z
%     for z_target = 100:10:240
%         % z_target = 150;
%         
%         figure
%         %     imagesc(template_stack(:,:,z_target))
%         C = imfuse(template_stack(:,:,z_target),output_stack(:,:,z_target));
%         image(C)
%         %     hold on
%         set(gca,'TickLength',[0 0])
%         %     scatter(rescaled_px(rescaled_px(:,3)==z_target,1),rescaled_px(rescaled_px(:,3)==z_target,2),0.1,'k')
%         pause(0.05)
%     end
%% Save for reformatting

% assemble the saving path
[~,name] = fileparts(obj_path);
save_path = fullfile(paths.registration_path,'AF_Segmentation_VAST\preregistration_stacks','mixed_masks.tif');

% save the file
for frames = 1:size(mixed_stack,3)
    if frames == 1
        oa = 'overwrite';
    else
        oa = 'append';
    end
    imwrite(mixed_stack(:,:,frames),save_path,'tif','WriteMode',oa)
    
end
toc