%% Create a map image for reformatting registration
%% Clean up
clearvars
close all
load('paths.mat')
addpath(genpath(paths(1).main_path))
%% Define the size of the image and other constants

% get the target stack
target_path = uipickfiles('FilterSpec',fullfile(paths.registration_path,'Pre_registration_brains','*tif'));

% define the zero threshold
zero_threshold = 200;
%% Load the template stack

% get the stack size
stack_info = imfinfo(target_path{1},'tif');
im_size = [stack_info(1).Height,stack_info(1).Width,size(stack_info,1)];

% allocate memory for the stack
template_stack = zeros(im_size);
% load the stack
% for all the frames
for frame = 1:im_size(3)
    template_stack(:,:,frame) = imread(target_path{1},frame);
end
%% Create the stack

% map_stack_r = zeros(im_size,'uint16');
% map_stack_g = zeros(im_size,'uint16');
% map_stack_b = zeros(im_size,'uint16');
% 
% 
% % for all the pixels in the stack
% for pix = 1:numel(map_stack_r)
%     % if the pixel is above the 0 threshold
%     if template_stack(pix) > zero_threshold
%     
%         % get the coordinates of the point
%         [x,y,z] = ind2sub(im_size,pix);
%         map_stack_r(pix) = uint16(x);
%         map_stack_g(pix) = uint16(y);
%         map_stack_b(pix) = uint16(z);
%     end
% end

% generate the maps
map_stack_02 = repmat(uint16(1:im_size(1))',[1 im_size(2) im_size(3)]);
map_stack_03 = repmat(uint16(1:im_size(2)),[im_size(1) 1 im_size(3)]);
map_stack_04 = repmat(uint16(reshape(1:im_size(3),1,1,[])),[im_size(1) im_size(2) 1]);
%% Save the stack

% get the file name
[~,name] = fileparts(target_path{1});
% define the output path
output_path_02 = fullfile(paths.registration_path,'Map_stacks',strcat(name,'_02.tif'));
output_path_03 = fullfile(paths.registration_path,'Map_stacks',strcat(name,'_03.tif'));
output_path_04 = fullfile(paths.registration_path,'Map_stacks',strcat(name,'_04.tif'));


% save the tif

% for all the frames
for frame = 1:im_size(3)
    if frame == 1
        imwrite(map_stack_02(:,:,1),output_path_02,'tif','WriteMode','overwrite')
        imwrite(map_stack_03(:,:,1),output_path_03,'tif','WriteMode','overwrite')
        imwrite(map_stack_04(:,:,1),output_path_04,'tif','WriteMode','overwrite')

    else
        imwrite(map_stack_02(:,:,frame),output_path_02,'tif','WriteMode','append')
        imwrite(map_stack_03(:,:,frame),output_path_03,'tif','WriteMode','append')
        imwrite(map_stack_04(:,:,frame),output_path_04,'tif','WriteMode','append')
    end
end