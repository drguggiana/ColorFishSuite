%% Put the reformatted cut section in the main reference brain
%% Clean up
clearvars
close all
load('paths.mat')
addpath(genpath(paths(1).main_path))
%% Load the reformatted cut stack

% define the mask path
mask_path = fullfile(paths.registration_path,...
    'AF_Segmentation_VAST\reformatted\7dpf_top_down_nocrop_02_warp_m0g15c4e1e-1x52r1.tif');

% get the stack info
mask_info = imfinfo(mask_path);

% allocate memory for the stack
mask_stack = zeros(mask_info(1).Height,mask_info(1).Width,size(mask_info,1));

% for all the frames
for frames = 1:size(mask_info,1)
    mask_stack(:,:,frames) = imread(mask_path,frames);
end
%% Load the ZBB-isl2-GFP ref brain

% define the mask path
ref_path = fullfile(paths.registration_path,...
    'Reference_brains','ZBB_isl2b-GFP.tif');

% get the stack info
ref_info = imfinfo(ref_path);

% allocate memory for the stack
ref_stack = zeros(ref_info(1).Height,ref_info(1).Width,size(ref_info,1));

% for all the frames
for frames = 1:size(ref_info,1)
    ref_stack(:,:,frames) = imread(ref_path,frames);
end
%% Insert into the full reference brain

% create the empty ref brain stack
full_mask = zeros(size(ref_stack));

% define the shifts between the cut and full brains
y_shift = 0;
x_shift = 288;
z_shift = 10;

% insert the masks in the empty stack
full_mask(x_shift:x_shift+size(mask_stack,1)-1,1:size(mask_stack,2),z_shift:end) = mask_stack;

% convert to uint16 for saving
full_mask = uint16(full_mask);
%% Plot an overlay of the stacks
close all

% define the target z
for z_target = 1:10:size(ref_stack,3)
    % z_target = 150;
    
    figure
    %     imagesc(template_stack(:,:,z_target))
    C = imfuse(full_mask(:,:,z_target),ref_stack(:,:,z_target));
    image(C)
    %     hold on
    set(gca,'TickLength',[0 0])
    %     scatter(rescaled_px(rescaled_px(:,3)==z_target,1),rescaled_px(rescaled_px(:,3)==z_target,2),0.1,'k')
    pause(0.05)
end
%% Save the mask

% define path
save_path = fullfile(paths.registration_path,'AF_Segmentation_VAST','full_mask','AF_mask.tif');

% save the stack
for frames = 1:size(ref_stack,3)
    if frames == 1
        oa = 'overwrite';
    else
        oa = 'append';
    end
    
    imwrite(full_mask(:,:,frames),save_path,'tif','WriteMode',oa)
    
end