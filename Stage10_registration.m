clearvars
close all

% addpath(genpath('E:\Behavioral data\Matlab'))
% addpath(genpath('C:\cygwin64\home\'))

% define the figure path
load('paths.mat')
addpath(genpath(paths(1).main_path))

cluster_path = paths(1).stage3_path;
fig_path = strcat(paths(1).fig_path,'Registration\');

data = load_clusters(cluster_path);

%get the folder where the image files are
tar_path_reg = uipickfiles('FilterSpec',paths(1).registration_path);
%% Load reformatted average stack
%get the number of experiments selected
num_files = length(tar_path_reg);

%allocate memory to store the maps
temp_tiff = uint32(loadtiff(tar_path_reg{1}));
all_maps = zeros(size(temp_tiff,1),size(temp_tiff,2),size(temp_tiff,3),num_files,'uint32');
%for all the experiments
for files = 1:num_files
    %load the tiff file
    temp_tiff = uint32(loadtiff(tar_path_reg{files}));
    
    all_maps(:,:,:,files) = temp_tiff;
    
end
%% Load transformation matrix

%translate
%rotate
%scale
%shear
%center

%define the main loading path
% tmat_path = 'C:\cygwin64\home\Drago\Registration\affine\C\Users\Drago\Desktop\registrationtrial\refh2bcut_C\Users\Drago\Desktop\registrationtrial\procbrains';
% tmat_path = 'C:\cygwin64\home\Drago\Registration\affine\C\Users\Drago\Desktop\registrationtrial\refh2bwrp2_C\Users\Drago\Desktop\registrationtrial\procbrains';
% tmat_path = 'C:\cygwin64\home\Drago\Registration\affine\C\Users\Drago\Desktop\registrationtrial\refh2bRESCALE_C\Users\Drago\Desktop\registrationtrial\procbrains';
tmat_path = 'C:\cygwin64\home\Drago\Registration\affine\C\Users\Drago\Desktop\registrationtrial\refh2bFINAL_C\Users\Drago\Desktop\registrationtrial\procbrains';

%allocate memory for the pa and final matrices
pa_reg_cell = cell(2,1);

 %allocate memory for the affine objects
 % aff_mat = zeros(4,4,num_exp);
 aff_cell = cell(num_files,1);

%for both types of regs
for regs = 1:2
    %go file by file extracting the matrix
    %allocate memory for the matrices
    x_mat = zeros(5,3,num_files);
    %for all the files
    for files = 1:num_files
        %get the file name
        [~,f_name] = fileparts(tar_path_reg{files});
        switch regs
            case 1
                %get rid of the extension and add the "list" extension to define the
                %target path
                f_name = strcat(f_name(1:end-9),'pa.list');
                %load the actual matrix
                temp_cell = importdata(fullfile(tmat_path,f_name,'registration'));
                %allocate memory for the parsed data
                parse_mat = zeros(5,2);
                
                parse_mat(1,:) = [str2double(temp_cell.textdata{3,2}) str2double(temp_cell.textdata{4,1})]';
                parse_mat(2,:) = [str2double(temp_cell.textdata{5,2}) str2double(temp_cell.textdata{6,1})]';
                parse_mat(3,:) = [str2double(temp_cell.textdata{7,2}) str2double(temp_cell.textdata{8,1})]';
                parse_mat(4,:) = [str2double(temp_cell.textdata{9,2}) str2double(temp_cell.textdata{10,1})]';
                parse_mat(5,:) = [str2double(temp_cell.textdata{11,2}) str2double(temp_cell.textdata{12,1})]';
                %store the parse matrix in the main matrix
                x_mat(:,:,files) = cat(2,parse_mat,[0;0;1;0;0]);
                
            case 2
                %get rid of the extension and add the "list" extension to define the
                %target path
                f_name = strcat(f_name(1:end-4),'list');
        
                %load the actual matrix
                temp_cell = importdata(fullfile(tmat_path,f_name,'registration'));
                %and parse it
                %allocate memory for the parsed data
                parse_mat = zeros(5,3);
                %for all the relevant lines
                for lines = 6:10
                    %split the string via spaces
                    temp_str = strsplit(temp_cell{lines},' ');
                    %convert the last 3 to numbers and store
                    parse_mat(lines-5,:) = [str2double(temp_str{2});str2double(temp_str{3});str2double(temp_str{4})];
                end
                %store the parse matrix in the main matrix
                x_mat(:,:,files) = parse_mat;
        end
        
    end   
    %% Turn the transformation matrix into an affine matrix

    %allocate memory for the transformation matrices
    trans_mats = cell(num_files,1);

    %for all the files
    for files = 1:num_files
        %turn the cmtk parameter matrix into an affine matrix and store
        trans_mats{files} = cmtkparams2affmat(x_mat(:,:,files));
        
    end
    %store the matrices
    pa_reg_cell{regs} = trans_mats;
end


%get the corresponding affine object
%for all the files
for files = 1:num_files
%     aff_cell{files} = affine3d(pa_reg_cell{2}{files}*pa_reg_cell{1}{files});
    aff_cell{files} = affine3d(pa_reg_cell{2}{files});
end
%% Load the reference brain

%define the path
ref_path = 'E:\Behavioral data\Matlab\AF_proc\Analysis\Refbrain_tif\';
%define the filename
% ref_name = 'refh2bcut.tif';
% ref_name = 'refh2bwrp2.tif';
% ref_name = 'refh2bRESCALE.tif';
ref_name = 'refh2bFINAL.tif';

%allocate memory for the stack
im_info = imfinfo(fullfile(ref_path,ref_name));
ref_stack = zeros(im_info(1).Height,im_info(1).Width,size(im_info,1));
%load the stack
%for all the z sections
for z = 1:size(im_info,1)
    ref_stack(:,:,z) = imread(fullfile(ref_path,ref_name),z);
end

%get the dimensions of the ref brain to allocate for the registered brains
ref_dim = size(ref_stack);
% ref_dim = [573 496 274];
% ref_dim([1 2]) = round(ref_dim([1 2]).*310/320);
% ref_dim = [621 718 138];
%% Load the pre-reg tif (from pre-reg NRRD file)

%define the path
prereg_path = 'E:\Behavioral data\Matlab\AF_proc\Analysis\Refbrain_tif\Stack_tif';
%get the tif files in the directory
tif_list = dir(strcat(prereg_path,'\*.tif'));
%get the number of tifs
tif_num = length(tif_list);
%allocate memory to store the prereg coords
prereg_coord = cell(tif_num,2);
%for all the files
for tifs = 1:tif_num
    
    %define the filename
%     prereg_name = '20160914_h2b6s_p17b_1_anato.tif_01.tif';
    prereg_name = tif_list(tifs).name;

    %allocate memory for the stack
    im_info = imfinfo(fullfile(prereg_path,prereg_name));
    prereg_stack = zeros(im_info(1).Height,im_info(1).Width,size(im_info,1));
    %load the stack
    %for all the z sections
    for z = 1:size(im_info,1)
        prereg_stack(:,:,z) = imread(fullfile(prereg_path,prereg_name),z);
    end

    %also extract its coordinates
    ind_coord = find(prereg_stack);
    [x,y,z] = ind2sub(size(prereg_stack),ind_coord);
    prereg_coord{tifs,1} = prereg_stack;
    prereg_coord{tifs,2} = cat(2,x,y,z);
end
%% Register the fish data
close all

% get the fish of origin
fish_ori = data.fish_ori;

% get the seed coordinates
cat_seed_all = data.xy_seed;
cat_z_all = data.z_seed;

% get the gains
delta_norm = data.delta_norm;

%select the target file
tar_file = ([1 2 3 4 5 6 7 8 9 10 11 12]);

% %define the new number of files
% new_num_files = length(tar_file);


%allocate memory for the registered fish
fishave_cell = cell(num_files,2);

%get the number of fish in the coordinate file
num_fish_coord = length(unique(fish_ori(:,1)));

% %leave only the selected fish in cat_seed_all
% cat_seed_all = cat_seed_all(include_vec);
 
%allocate memory for the final structure
%get the sizes of the structure
%         seeds_sizes = vertcat(cat_seed_all(:).area);
seeds_px = vertcat(cat_seed_all(:).pxlist);
fish_xyz = zeros(size(seeds_px,1),3);

%initialize a seed counter
seed_c = 1;
%and also allocate memory for a vector to replace fish_ori
new_ori = zeros(size(seeds_px,1),1);
%and the stimulus vector
new_stim = zeros(size(seeds_px,1),4);
%for all the seeds
for seeds = 1:size(cat_seed_all,1)
    %get the number of pixels in this seeds
    num_pix = cat_seed_all(seeds).area;
    %and the actual pixels
    list_pix = cat_seed_all(seeds).pxlist;
    %transform the indexes to subindexes
    [sub_y,sub_x] = ind2sub([320 320],list_pix);
    %flip the x indexes
    sub_y = 320 - sub_y + 1;
    %load the transformed indexes and the z into the fish_xyz
    %matrix
    fish_xyz(seed_c:seed_c+num_pix-1,:) = ...
        [sub_x,sub_y,ones(num_pix,1).*cat_z_all(seeds)];
    %load the fish of origin and stimulus into the new vector
    new_ori(seed_c:seed_c+num_pix-1) = fish_ori(seeds,1);
    %do the same with the stimulus vector
    new_stim(seed_c:seed_c+num_pix-1,:) = repmat(delta_norm(seeds,:),num_pix,1);
    
    %update the index
    seed_c = seed_c + num_pix;
    
end

%invert z
fish_xyz(:,3) = max(fish_xyz(:,3)) - fish_xyz(:,3)+1;

%renumber z
fish_xyz(:,3) = mod(fish_xyz(:,3)-1,max(cat_z_all)/num_fish_coord)+1;

% %use the prereg coord instead of the fish
% fish_xyz = prereg_coord{tar_file,2};

%invert x and y
fish_xyz(:,[2 1]) = fish_xyz(:,[1 2]);

%rescale dimensions based on the micron/pixel divergence
fish_xyz(:,[1 2]) = (fish_xyz(:,[1 2])-5)*0.906;
fish_xyz(:,3) = fish_xyz(:,3)*5;
    
%initialize a counter for the fish
fish_c = 1;
%for all the files
for files = tar_file%:num_files
    %get the coordinates for the current fish
    curr_fish = fish_xyz(new_ori==files,:);
    
%     %use the prereg coord instead of the fish    
%     curr_fish = fish_xyz;

    %transform the volume
    new_fishave_xyz = transformPointsInverse(aff_cell{fish_c},curr_fish);
    
%     %don't transform
%     new_fishave_xyz = curr_fish;
    %bring the coordinates back to pixels
    new_fishave_xyz(:,[1 2]) = (new_fishave_xyz(:,[1 2]))/(495.56/621);
    new_fishave_xyz(:,3) = new_fishave_xyz(:,3)/2;

    %store the data in the fish cell
    fishave_cell{fish_c,1} = new_fishave_xyz; 
    fishave_cell{fish_c,2} = new_stim(new_ori==files,:);
    
    %update the counter
    fish_c = fish_c + 1;
    
end

%plot the average
%allocate memory for the matrix
new_volume = zeros(horzcat(ref_dim,num_files));
%initialize a fish counter
fish_c = 1;
%for all fish
for files  = tar_file%:num_files
    %round the new coordinates
    curr_seeds = round(fishave_cell{fish_c,1});
    
    %clip the seeds that fall outside the volume in z, and show a flag
    over_seeds = curr_seeds(:,3)>ref_dim(3) | (curr_seeds(:,3)<1) | ...
        curr_seeds(:,1)>ref_dim(2) | curr_seeds(:,1)<1 |...
        curr_seeds(:,2)>ref_dim(1) | curr_seeds(:,2)<1;
    %if there were seeds outside the boundaries, exclude them
    if sum(over_seeds) > 0
        curr_seeds = curr_seeds(~over_seeds,:);
        fprintf(strcat('file:',num2str(files),'\r\n',...
            'excluded over lim seeds:',num2str(sum(over_seeds)),'\r\n'))
    end
    %load the data in the common volume
    new_volume(sub2ind(size(new_volume),curr_seeds(:,2),curr_seeds(:,1),...
        curr_seeds(:,3),ones(length(curr_seeds),1).*fish_c)) = 1;
    %update the counter
    fish_c = fish_c + 1;
    
end
%% Plot the fish overlap with the registered average volume
close all
% %define the volume/slice to plot
% tar_slice = 1;

%for all the fish

for files = tar_file
    figure
    %for all three dimensions to slice
    for dims = 1:3
        ave_stack = squeeze(sum(all_maps(:,:,:,files),dims));
    %     ave_stack = squeeze(sum(all_maps(:,:,tar_slice,tar_file),dims));
    %     ave_stack = squeeze(sum(ref_stack,dims));
    %     reg_stack = squeeze(sum(ref_stack,dims));
        reg_stack = squeeze(sum(permute(new_volume(:,:,:,files),[1 2 3]),dims));
    %     reg_stack = squeeze(sum(permute(new_volume(:,:,tar_slice,tar_file),[2 1 3]),dims));
    %     reg_stack = squeeze(sum(prereg_stack,dims));
    %     %flip reg_stack
    %     reg_stack = reg_stack(:,size(reg_stack,2):-1:1,:);


%         %plot as fused image
%         comb_2 = imfuse(ave_stack,reg_stack,'falsecolor');
        %plot as binary
        blank = zeros(size(ave_stack));
        ave_stack = cat(3,blank,ave_stack,blank);
        reg_stack = repmat(reg_stack,1,1,3);
        comb_2 = normr_1(ave_stack,1) + reg_stack;
        %select the subplot (projection)
        switch dims
            case 1
                subplot(2,5,1:3)
                imagesc(permute(comb_2,[2 1 3]))
            case 2
                subplot(2,5,6:8)
                imagesc(permute(comb_2,[2 1 3]))
            case 3
                subplot(2,5,[4:5 9:10])
                imagesc(comb_2)
        end
        axis equal
    end
end
autoArrangeFigures
%% Plot an average projection

close all

%for all three dimensions to slice
for dims = 1:3
%     ave_stack = squeeze(mean(mean(all_maps,4),dims));
    %     ave_stack = squeeze(sum(all_maps(:,:,tar_slice,tar_file),dims));
        ave_stack = squeeze(mean(ref_stack,dims));
%         reg_stack = squeeze(sum(ref_stack,dims));
    reg_stack = squeeze(mean(mean(new_volume,4),dims));
    %     reg_stack = squeeze(sum(permute(new_volume(:,:,tar_slice,tar_file),[2 1 3]),dims));
    %     reg_stack = squeeze(sum(prereg_stack,dims));
    %     %flip reg_stack
    %     reg_stack = reg_stack(:,size(reg_stack,2):-1:1,:);
    
    
%     %plot as fused image
%     comb_2 = imfuse(ave_stack,reg_stack,'falsecolor');
    %plot as binary
    blank = zeros(size(ave_stack));
    ave_stack = cat(3,blank,ave_stack,blank);
    reg_stack = cat(3,reg_stack,blank,blank);
    comb_2 = normr_1(ave_stack,1) + normr_1(reg_stack,1);
    %select the subplot (projection)
    switch dims
        case 1
            subplot(2,5,1:3)
            imagesc(permute(comb_2,[2 1 3]))
        case 2
            subplot(2,5,6:8)
            imagesc(permute(comb_2,[2 1 3]))
        case 3
            subplot(2,5,[4:5 9:10])
            imagesc(comb_2)
    end
end
%% Calculate maps based on gains
close all
%define the target percentile
perc = 90;
%allocate memory for the matrix
new_volume = zeros(horzcat(ref_dim,num_files,4));
%for all fish
for files  = tar_file%:num_files
    %round the new coordinates
    curr_seeds = round(fishave_cell{files,1});
    curr_act = fishave_cell{files,2};
    
    %clip the seeds that fall outside the volume in z, and show a flag
    over_seeds = curr_seeds(:,3)>ref_dim(3) | (curr_seeds(:,3)<1) | ...
        curr_seeds(:,1)>ref_dim(2) | curr_seeds(:,1)<1 |...
        curr_seeds(:,2)>ref_dim(1) | curr_seeds(:,2)<1;
    if sum(over_seeds) > 0
        curr_seeds = curr_seeds(~over_seeds,:);
        curr_act = curr_act(~over_seeds,:);
        fprintf(strcat('file:',num2str(files),'\r\n',...
            'excluded over lim seeds:',num2str(sum(over_seeds)),'\r\n'))
    end
    %for the 4 cone types
    for cones = 1:4
        %threshold for gains only in the highest percentile of the
        %distributions
        perc_thres = prctile(curr_act(curr_act(:,cones)>0,cones),perc);
        curr_act(curr_act(:,cones)<perc_thres,cones) = 0;
        new_volume(sub2ind(size(new_volume),curr_seeds(:,2),curr_seeds(:,1),...
            curr_seeds(:,3),ones(length(curr_seeds),1).*files,...
            ones(length(curr_seeds),1).*cones)) = curr_act(:,cones);
    end
end
%% Plot activity maps
close all
%take the average across animals 
animal_ave = squeeze(sum(new_volume,4));
% %make the elements that are zero across the map NaNs
% %get the indexes
% nan_idx = repmat(sum(animal_ave,4),1,1,1,4);
% animal_ave(nan_idx==0) = NaN;

figure
%plot the average across animals
for cones = 1:4
    subplot(2,2,cones)
    imagesc(squeeze(mean(animal_ave(:,:,:,cones),3)))
end

%calculate a map of the pixels with high gain
high_gain = squeeze(sum(animal_ave>0,4));
projection_plot(high_gain);

%leave only the max gain as non-zero pixel
[~,max_ind] = max(animal_ave,[],4);
%correct the indexing for the zero pixels
max_ind = max_ind.*(sum(animal_ave,4)~=0);

figure
histogram(max_ind(:))

% figure
%allocate memory for the map
max_map = zeros(horzcat(size(max_ind),3));
%for all the cones
for cones = 1:4
    %load the current color
    curr_color = animal_ave(:,:,:,cones);
%     %get a matrix to store the 3D color data
%     temp_color = zeros(size(curr_color));
    %get the indexing matrix
    idx_mat = max_ind==cones;
    switch cones
        case 1
            max_map(:,:,:,1) = curr_color.*idx_mat + max_map(:,:,:,1);
        case 2
            max_map(:,:,:,2) = curr_color.*idx_mat + max_map(:,:,:,2);
        case 3
            max_map(:,:,:,3) = curr_color.*idx_mat + max_map(:,:,:,3);
        case 4
            max_map(:,:,:,1) = curr_color.*idx_mat + max_map(:,:,:,1);
            max_map(:,:,:,3) = curr_color.*idx_mat + max_map(:,:,:,3);
    end
%     subplot(2,2,cones)
%     imagesc(sum(curr_color.*idx_mat,3))
%     axis equal
projection_plot(curr_color.*idx_mat);
end

projection_plot(max_map);
clear('max_map','animal_ave','max_ind','high_gain','curr_color')
autoArrangeFigures
%% Plot standard deviation maps

close all
%take the std across animals
animal_std = squeeze(sum(zscore(new_volume,0,4),4));
% animal_std = squeeze(std(new_volume,0,4));
% %make the elements that are zero across the map NaNs
% %get the indexes
% nan_idx = repmat(sum(animal_ave,4),1,1,1,4);
% animal_ave(nan_idx==0) = NaN;

figure
%plot the average across animals
for cones = 1:4
    subplot(2,2,cones)
    imagesc(squeeze(mean(animal_std(:,:,:,cones),3)))
end

%calculate a map of the pixels with high std across cones
high_std = squeeze(sum(animal_std>0,4));
projection_plot(high_std);

%leave only the max std as non-zero pixel
[~,max_ind] = max(animal_std,[],4);
%correct the indexing for the zero pixels
max_ind = max_ind.*(sum(animal_std,4)~=0);

% figure
% histogram(max_ind(:))

figure
%allocate memory for the map
max_map = zeros(horzcat(size(max_ind),3));
%for all the cones
for cones = 1:4
    %load the current color
    curr_color = normr_1(animal_std(:,:,:,cones),1);
%     %get a matrix to store the 3D color data
%     temp_color = zeros(size(curr_color));
    %get the indexing matrix
    idx_mat = max_ind==cones;
    switch cones
        case {1,2,3}
            max_map(:,:,:,cones) = curr_color.*idx_mat + max_map(:,:,:,cones);
        case 4
            max_map(:,:,:,1) = curr_color.*idx_mat + max_map(:,:,:,1);
            max_map(:,:,:,3) = curr_color.*idx_mat + max_map(:,:,:,3);
    end
    subplot(2,2,cones)
    imagesc(sum(curr_color.*idx_mat,3))
    axis equal
end

projection_plot(max_map);
clear('max_map','animal_std','max_ind','high_std','curr_color')
autoArrangeFigures
%% Calculate variation in tuning for each pixel
close all
%get the standard deviation across colors and average across animals
col_std = squeeze(mean(std(new_volume,0,5),4));

%plot the projection combined with the ref stack
projection_plot(col_std,ref_stack);
%clear memory
clear('col_std')
%% calculate z projection
z_projection = squeeze(sum(new_volume,3));

% also collapse data into a single matrix
new_fish_xyz = new_fishave_xyz;

% %for all the files
% for files = 1:num_files
%     imagesc(z_projection(:,:,files))
% %     hold('on')
% end
%% Save a tiff stack with the desired data

%get the target volume
print_stack = max_map;
%add across fish too
% r_img = squeeze(mean(print_stack,4));
% r_img(r_img<0) = 0;
r_img = uint8(255.*normr_1(print_stack,1));

%and store as tif stack
fig_path = 'E:\Behavioral data\Matlab\AF_proc\Analysis\Figures';
fig_name = strcat('Max_Intensity');
fig_full = fullfile(fig_path,strcat(fig_name,'.tif'));

imwrite(squeeze(r_img(:,:,1,:)),fig_full,'tif','WriteMode','overwrite')

for z = 2:size(r_img,3)
    
    
    imwrite(squeeze(r_img(:,:,z,:)),fig_full,'tif','Resolution',[size(r_img,1),size(r_img,2)],'WriteMode','append')
end
%% Compare the ref brain to the registered data

close all
%get the max projection of the ref brain
max_ref = max(ref_stack,[],3);

% %rescale the ref
% max_ref = imresize(max_ref,ref_dim(1:2));
%create a combined image
figure
% comb_im = imfuse(max_ref,fish_projection,'falsecolor');
comb_im = imfuse(max_ref,z_projection(:,:,1),'falsecolor');
image(comb_im)
axis equal
%% Calculate intensity maps for all the stimuli

close all

%round the new coordinates
work_fish_xyz = round(new_fish_xyz);


% load the traces
conc_trace = data.conc_trace;
% load stim num
stim_num2 = data.stim_num;

%create a new matrix capable of containing the registered data

% %get the bounds x y and z of the registered data
% max_all = max(work_fish_xyz);
% min_all = min(work_fish_xyz);
% 
% %bring the coordinates to the positive
% work_fish_xyz = work_fish_xyz - min_all + 1;

%reshape the response matrix
reshape_trace = reshape(conc_trace,size(conc_trace,1),stim_num2,[]);
figure
%allocate memory to store the projections
% project_mat = zeros(horzcat(max_all(1:2)-min_all(1:2)+1,stim_num2));
project_mat = zeros(horzcat(ref_dim([1 2]),stim_num2));
%also allocate memory to store the whole volumes
volume_cell = cell(stim_num2,1);
%for all the stimuli
for stim = 1:stim_num2
    
%     %calculate the average response during the stimulation period for this
%     %stim
%     stim_vec = squeeze(mean(abs(reshape_trace(:,stim,21:60)),3));
    %or just extract the gains for the stimulus
    stim_vec = delta_norm(:,stim);
    %and center their distributions
    stim_vec = ((stim_vec - mean(stim_vec))./std(stim_vec));
    subplot(2,2,stim)
    histogram(stim_vec,100)
    %allocate memory for the matrix
%     new_volume = zeros(horzcat(max_all-min_all + 1,num_files));
    new_volume = zeros(horzcat(ref_dim,num_files));

    %for all fish
    for files  = 1:num_files
        %load the corresponding seeds in the target volume
        curr_seeds = work_fish_xyz(fish_ori(:,1)==files,:);
        %and the corresponding responses
        curr_stim = stim_vec(fish_ori(:,1)==files);
        %clip the seeds that fall outside the volume in z, and show a flag
        over_seeds = (curr_seeds(:,3)>ref_dim(3)) | (curr_seeds(:,3)<1);
        if sum(over_seeds) > 0
            curr_seeds = curr_seeds(~over_seeds,:);
            curr_stim = curr_stim(~over_seeds);
            fprintf(strcat('file:',num2str(files),'\r\n',...
                'excluded over z seeds:',num2str(sum(over_seeds)),'\r\n'))
        end
        
        new_volume(sub2ind(size(new_volume),curr_seeds(:,1),curr_seeds(:,2),...
            curr_seeds(:,3),ones(length(curr_seeds),1).*files)) = curr_stim;
        
    end
     
    %save the volume
    volume_cell{stim} = new_volume;

    %plot the fish overlapping each other
    %calculate z projection
%     z_projection = squeeze(max(new_volume,[],3));

    z_projection = squeeze(max(new_volume,[],3));
 
    % %for all the files
    % for files = 1:num_files
    %     imagesc(z_projection(:,:,files))
    % %     hold('on')
    % end

    %add across fish too
    fish_projection = squeeze(max(z_projection,[],3));

%     subplot(2,2,stim)
%     imagesc(fish_projection)
    
    %store the normalized projection
    project_mat(:,:,stim) = normr_1(fish_projection,1);
end
error('Stop Here')
%% Save average tiff stacks

close all

%for all the stimuli
for stim = 1:stim_num2

    %get the target volume
    new_volume = volume_cell{stim};
    %add across fish too
    r_img = squeeze(mean(new_volume,4));
    r_img(r_img<0) = 0;
    r_img = uint16(255.*normr_1(r_img,1));

    %and store as tif stack
     fig_path = 'E:\Behavioral data\Matlab\AF_proc\Figures\20180103_maps';
    fig_name = strcat('AverageStimulus_',num2str(stim));
    fig_full = fullfile(fig_path,strcat(fig_name,'.tif'));

    imwrite(squeeze(r_img(:,:,1)),fig_full,'tif','WriteMode','overwrite')

    for z = 2:size(r_img,3)


        imwrite(squeeze(r_img(:,:,z)),fig_full,'tif','Resolution',[size(r_img,1),size(r_img,2)],'WriteMode','append')
    end
end
%% Plot the projections in 4 channels
close all
figure

%allocate memory for the combined map
comb_map = zeros(size(project_mat,1),size(project_mat,2),3);

%load the colors in the different channels
comb_map(:,:,1:3) = project_mat(:,:,1:3);

imagesc((comb_map))
axis equal
figure
%for all the stimuli
for stim = 1:stim_num2
   subplot(2,2,stim)
   imagesc(project_mat(:,:,stim))
   axis equal
end
%% Threshold the values per stimulus

%allocate memory for the thresholded matrix
thres_mat = zeros(size(project_mat));
%for all the stimuli
for stim = 1:stim_num2
    %get the current stim
    curr_stim = project_mat(:,:,stim);
    
    %get the values only (to get the threshold
    vals = curr_stim(curr_stim~=0);

    %get the lowest percentile
    thres = prctile(vals,50);
    %threshold the matrix
    temp_mat = curr_stim;
    temp_mat(temp_mat<thres) = 0;
    
    %save the results
    thres_mat(:,:,stim) = temp_mat;
end
%% Plot only the highest signal out of the 4 channels
close all

%get the value for each color
[max_val,max_ind] =  max(thres_mat,[],3);
%get the filled pixels
empty_ind = sum(thres_mat,3)~=0;

%combine the 2 maps
max_filled = max_ind.*empty_ind;

figure
imagesc(max_filled)
axis equal
col_map = [0 0 0;1 0 0;0 1 0;0 0 1;1 0 1];
colormap(col_map)
%% Plot pairwise distance distributions for each stimulus

close all

%allocate memory for the distributions
dist_cell = cell(stim_num2,1);
%allocate memory to store the fits
fit_cell = cell(stim_num2,1);
%for each stimulus
for stim = 1:stim_num2
    %extract the gains for the stimulus
    stim_vec = delta_norm(:,stim);
    %and center their distributions
    stim_vec = ((stim_vec - mean(stim_vec))./std(stim_vec));
    %get the lowest percentile to filter out
    thres = prctile(stim_vec,50);
    thres_vec = stim_vec>thres;
    %allocate memory to store the distances per fish
    stim_dist = cell(num_files,1);
    
    %for all fish
    for files  = 1:num_files
        %load the corresponding seeds in the target volume
%         curr_seeds = work_fish_xyz(fish_ori(:,1)==files&thres_vec,:);
        curr_seeds = fishave_cell{files,1}(thres_vec,:);
        %get the pairwise distances
        stim_dist{files} = pdist(curr_seeds)';
        
        
    end
    %concatenate across fish and save
    dist_cell{stim} = vertcat(stim_dist{:});
    
    %also plot
    subplot(2,2,stim)
%     histogram(dist_cell{stim},100)
    temp_vec = dist_cell{stim};
    temp_vec = temp_vec(temp_vec~=0);
    params = histfit(temp_vec,100,'kernel');

    %save the datapoints from the non parametric fit
    fit_cell{stim} = get(params(2),'YData');
end
%% Generate null pairwise distance distributions

%define the number of iterations
num_iter = 10;

%allocate memory for the fits
null_cell = cell(stim_num2,num_iter);


%for each stimulus
for stim = 1:stim_num2
    %extract the gains for the stimulus
    stim_vec = delta_norm(:,stim);
    %and center their distributions
    stim_vec = ((stim_vec - mean(stim_vec))./std(stim_vec));
    %get the lowest percentile to filter out
    thres = prctile(stim_vec,50);
    thres_vec = stim_vec>thres;
    
    %for all the iterations
    for iter = 1:num_iter
        %allocate memory to store the distances per fish
        stim_dist = cell(num_files,1);
        
        %for all fish
        for files  = 1:num_files
            %get the number of items to draw
            num_items = sum(fish_ori(:,1)==files&thres_vec);
            %generate a set of random indices for the required number of
            %seeds
%             rand_ind = randperm(prod(ref_dim),num_items);
            rand_ind = randperm(size(work_fish_xyz,1),num_items);
            %load the corresponding seeds in the target volume
            curr_seeds = work_fish_xyz(rand_ind,:);
            %get the pairwise distances
            stim_dist{files} = pdist(curr_seeds)';


        end
        %concatenate across fish and save
        temp_vec = vertcat(stim_dist{:});

%         %also plot
%         subplot(2,2,stim)
%     %     histogram(dist_cell{stim},100)
%         temp_vec = dist_cell{stim};
        temp_vec = temp_vec(temp_vec~=0);
        params = histfit(temp_vec,100,'kernel');

        %save the datapoints from the non parametric fit
        null_cell{stim,iter} = get(params(2),'YData');
    end
end

%allocate memory for the average distributions
null_ave = cell(stim_num2,1);
%average the null distribution reps
%for all the stimuli
for stim = 1:stim_num2
    %concatenate and average
    null_ave{stim} = mean(cat(1,null_cell{stim,:}),1);
end
%% Plot and fit the distributions
close all
figure
%for all the stimuli
for stim = 1:stim_num2
    plot(fit_cell{stim})
    hold('on')
    plot(null_ave{stim})
end

% %combine the data and null cells
% comb_cell = {fit_cell{:},null_ave{:}};
% figure
% comb_list = combnk(1:stim_num2*2,2);
% comb_im = zeros(stim_num2*2);
% for combs = 1:length(comb_list)
%     comb_im(comb_list(combs,1),comb_list(combs,2)) = immse(comb_cell{comb_list(combs,1)}',comb_cell{comb_list(combs,2)}');
% end
% imagesc(comb_im)
% set(gca,'XTick',[1 2 3 4],'XTickLabels',{'R','G','B','U'},'YTick',[1 2 3 4],'YTickLabels',{'R','G','B','U'})

%compare immse to null

figure
%allocate memory to store the values
stim_immse = zeros(stim_num2,1);
%for all the stimuli
for stim = 1:stim_num2
    %calculate the immse between real and null
    stim_immse(stim) = immse(fit_cell{stim}',null_ave{stim}');
end

plot(stim_immse,'*')
set(gca,'XTick',[1 2 3 4],'XTickLabels',{'R','G','B','U'})
%% More plots

%variability maps
%difference between space across animals
%compare pairwise distance distributions to null
%% Load the clusters

%get the folder where the image files are
tar_path_clu = uipickfiles('FilterSpec','E:\Behavioral data\Matlab\AF_proc\Analysis\Stage3_cluster\*.mat');

%load the cluster number
clu_num = load(tar_path_clu{1},'clu_num');
clu_num = clu_num.clu_num;
%and the cluster index
idx_clu = load(tar_path_clu{1},'idx_clu');
idx_clu = idx_clu.idx_clu;
%% Color seeds based on cluster and save tiff stack

close all

%round the new coordinates
work_fish_xyz = round(new_fish_xyz);

% %reshape the response matrix
% reshape_trace = reshape(conc_trace,size(conc_trace,1),stim_num2,[]);
% figure
%allocate memory to store the projections
% project_mat = zeros(horzcat(max_all(1:2)-min_all(1:2)+1,stim_num2));
% project_mat = zeros(horzcat(ref_dim([1 3]),stim_num2));
% %for all the stimuli
% for stim = 1:stim_num2
    
%     %calculate the average response during the stimulation period for this
%     %stim
%     stim_vec = squeeze(mean(abs(reshape_trace(:,stim,21:60)),3));
    %or just extract the gains for the stimulus
    clu_vec = idx_clu;
%     %and center their distributions
%     stim_vec = ((stim_vec - mean(stim_vec))./std(stim_vec));
%     subplot(2,2,stim)
%     histogram(stim_vec,100)
    %allocate memory for the matrix
%     new_volume = zeros(horzcat(max_all-min_all + 1,num_files));
    new_volume = zeros(horzcat(ref_dim,num_files));

    %for all fish
    for files  = 3%:num_files
        %load the corresponding seeds in the target volume
        curr_seeds = work_fish_xyz(fish_ori(:,1)==files,:);
        %and the corresponding responses
        curr_clu = clu_vec(fish_ori(:,1)==files);
        %clip the seeds that fall outside the volume in z, and show a flag
        over_seeds = curr_seeds(:,3)>ref_dim(3);
        if sum(over_seeds) > 0
            curr_seeds = curr_seeds(~over_seeds,:);
            curr_clu = curr_clu(~over_seeds);
            fprintf(strcat('file:',num2str(files),'\r\n',...
                'excluded over z seeds:',num2str(sum(over_seeds)),'\r\n'))
        end
        
        new_volume(sub2ind(size(new_volume),curr_seeds(:,1),curr_seeds(:,2),...
            curr_seeds(:,3),ones(length(curr_seeds),1).*files)) = 1;

    end

    %plot the fish overlapping each other
    %calculate z projection
%     z_projection = squeeze(max(new_volume,[],3));

%     z_projection = squeeze(max(new_volume,[],2));
 
    % %for all the files
    % for files = 1:num_files
    %     imagesc(z_projection(:,:,files))
    % %     hold('on')
    % end

    %add across fish too
    r_img = squeeze(mean(new_volume,4));
    r_img = uint16(255.*normr_1(r_img,1));

    %and store as tif stack
     fig_path = 'E:\Behavioral data\Matlab\AF_proc\Figures\20180103_maps';
    fig_name = 'test3';
    fig_full = fullfile(fig_path,strcat(fig_name,'.tif'));

    imwrite(squeeze(r_img(:,:,1)),fig_full,'tif','WriteMode','overwrite')

for z = 2:size(r_img,3)


    imwrite(squeeze(r_img(:,:,z)),fig_full,'tif','Resolution',[size(r_img,1),size(r_img,2)],'WriteMode','append')
end

%     subplot(2,2,stim)
%     imagesc(fish_projection)
    
%     %store the normalized projection
%     project_mat(:,:,stim) = normr_1(fish_projection,1);
% end


%% calculate mode and std deviation
temp_tiff = std(double(all_maps),0,4);
% temp_tiff = uint32(temp_tiff);
% r_img = zeros(size(temp_tiff,1),size(temp_tiff,2),size(temp_tiff,3),3,'uint32');
%     r_img(:,:,:,3) = mod(temp_tiff, 256);
%     t = (temp_tiff - r_img(:,:,:,3)) / 256;
%     r_img(:,:,:,2) = mod(t, 256);
%     t = (t - r_img(:,:,:,2)) / 256;
%     r_img(:,:,:,1) = t;
%     max(r_img(:))
%     
%     r_img = double(r_img);

r_img = temp_tiff;

r_img = uint16(255.*normr_1(r_img,1));



    %%
%     [fig_name,fig_path] = uiputfile('E:\Behavioral data\Matlab\AF_proc\Figures\*.*');
%         fig_path = 'E:\Behavioral data\Matlab\AF_proc\Figures\';
    fig_path = 'E:\Behavioral data\Matlab\AF_proc\Figures\20171111_maps';
    fig_name = 'test';
    fig_full = fullfile(fig_path,strcat(fig_name,'_.tif'));
%     temp_name = strsplit(reg_fold{fish},'.');
%     fig_full = fullfile(fig_path,strcat(temp_name{1},'.tif'));
    
    %         %if a grayscale image is desired
%     C = uint32(C(:,:,:,1)*2^16+C(:,:,:,2)*2^8+C(:,:,:,3));
%     options.overwrite = true;
%     options.append = false;
%     saveastiff(C(:,:,1),fig_full,options)
%     options.overwrite = false;
%     options.append = true;
%             imwrite(squeeze(r_img(:,:,1,:)),fig_full,'tif','WriteMode','overwrite')
                        imwrite(squeeze(r_img(:,:,1)),fig_full,'tif','WriteMode','overwrite')

for z = 2:size(r_img,3)
%         saveastiff(C(:,:,z),fig_full,options)
%         imwrite(squeeze(r_img(:,:,z,:)),fig_full,'tif','Resolution',[size(r_img,1),size(r_img,2)],'WriteMode','append')

            imwrite(squeeze(r_img(:,:,z)),fig_full,'tif','Resolution',[size(r_img,1),size(r_img,2)],'WriteMode','append')
end