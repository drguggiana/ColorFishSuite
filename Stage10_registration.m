clearvars
close all

% define the figure path
load('paths.mat')
addpath(genpath(paths(1).main_path))

cluster_path = paths(1).stage3_path;
fig_path = strcat(paths(1).fig_path,'Registration\');
registration_path = paths(1).registration_path;

data = load_clusters(cluster_path);

%get the folder where the image files are
reformatted_maps = dir(fullfile(registration_path,'Reformatted_brains',data.name,'*.tif'));

% turn the structure into just the paths
%get the number of experiments selected
num_files = length(reformatted_maps);
% allocate memory to store the paths
reformatted_cell = cell(num_files,1);
for files = 1:num_files
    reformatted_cell{files} = fullfile(reformatted_maps(files).folder,...
        reformatted_maps(files).name);
end
%% Load reformatted average stack


%allocate memory to store the maps
temp_tiff = uint32(loadtiff(reformatted_cell{1}));
all_maps = zeros(size(temp_tiff,1),size(temp_tiff,2),size(temp_tiff,3),num_files,'uint32');
%for all the experiments
for files = 1:num_files
    %load the tiff file
    temp_tiff = uint32(loadtiff(reformatted_cell{files}));
    
    all_maps(:,:,:,files) = temp_tiff;
    
end
%% Load transformation matrix

%translate
%rotate
%scale
%shear
%center

%define the main loading path
tmat_path = fullfile(registration_path,'Registration_info',data.name);

%allocate memory for the affine objects
aff_cell = cell(num_files,1);

%go file by file extracting the matrix
%allocate memory for the matrices
x_mat = zeros(5,3,num_files);
%for all the files
for files = 1:num_files
    %get the file name
    [~,f_name,~] = fileparts(reformatted_cell{files});
    [~,f_name] = fileparts(f_name);

    %get rid of the extension and add the "list" extension to define the
    %target path
    f_name = strcat(f_name,'.list');

    %load the actual matrix
    temp_cell = importdata(fullfile(tmat_path,f_name,'registration'));
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

end   
%% Turn the transformation matrix into an affine matrix

%allocate memory for the transformation matrices
trans_mats = cell(num_files,1);

%for all the files
for files = 1:num_files
    %turn the cmtk parameter matrix into an affine matrix and store
    trans_mats{files} = cmtkparams2affmat(x_mat(:,:,files));

end

%get the corresponding affine object
%for all the files
for files = 1:num_files
    aff_cell{files} = affine3d(trans_mats{files});
end
%% Load the reference brain

% define the ref brain depending on the file
if contains(data(1).name,{'Syn','syn'})
    ref_name = 'refisl2cut2.nrrd.tif';
    full_ref_name = 'refisl2.nrrd.tif';
%     ref_name = 'refcutblursub.nrrd.tif';
%     full_ref_name = 'refbrain.nrrd.tif';
else
    ref_name = 'refcutblursub.nrrd.tif';
    full_ref_name = 'refbrain.nrrd.tif';
end
% assemble the reference path
ref_path = fullfile(registration_path,'Reference_brains',ref_name);
full_ref_path = fullfile(registration_path,'Reference_brains',full_ref_name);
%allocate memory for the stack
ref_info = imfinfo(ref_path);
ref_stack = zeros(ref_info(1).Height,ref_info(1).Width,size(ref_info,1));

full_ref_info = imfinfo(full_ref_path);
full_ref_stack = zeros(full_ref_info(1).Height,full_ref_info(1).Width,size(full_ref_info,1));
%load the stack
%for all the z sections
for z = 1:size(ref_info,1)
%     ref_stack(:,:,z) = imread(fullfile(ref_path,ref_name),z);
    ref_stack(:,:,z) = imread(ref_path,z);
end

%for all the z sections
for z = 1:size(full_ref_info,1)
%     ref_stack(:,:,z) = imread(fullfile(ref_path,ref_name),z);
    full_ref_stack(:,:,z) = imread(full_ref_path,z);
end

%get the dimensions of the ref brain to allocate for the registered brains
ref_dim = size(ref_stack);
full_ref_dim = size(full_ref_stack);
%% Load the pre-reg tif (from pre-reg NRRD file)

%define the path
prereg_search = fullfile(registration_path,'Pre_registration_brains',data.name,'*.tif');
%get the tif files in the directory
tif_list = dir(prereg_search);
%get the number of tifs
tif_num = length(tif_list);
%allocate memory to store the prereg coords
prereg_coord = cell(tif_num,2);
% also for the im_info
im_info_cell = cell(tif_num,1);
% get the path
prereg_path = tif_list(1).folder;
%for all the files
for tifs = 1:tif_num
    
    %define the filename
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
    % store the stack information
    im_info_cell{tifs} = im_info;
end
%% Load the labels

% define the path to the file
labels_file = fullfile(registration_path,'Labels_info','MaskDatabase.mat');
% load the labels file
labels_data = load(labels_file);

% define which dataset to load depending on the reference
if contains(data.name,{'Syn','syn'})
    field_list = {'AF4','AF5','AF6','AF7','AF8','AF9','Tecum Neuropil','Periventriculare'};
    coordinate_conversion = [274,0,0];
else
    field_list = {'Tectum Stratum','Tecum Neuropil','Pretectum','Habenula','Cerebellum'};
    coordinate_conversion = [272,74,49];
end

% get the number of fields
field_number = length(field_list);
% allocate memory for the fields
label_cell = cell(field_number,2);
label_stack = zeros(labels_data.height,labels_data.width,labels_data.Zs);
% get the field numbers for these names
for field = 1:length(field_list)
    % get the index with the first name match (so as to not get
    % subdivisions)
    idx_vector = find(contains(labels_data.MaskDatabaseNames,field_list{field}),1);
    % store the map and the name in the cell
    label_stack = label_stack + reshape(full(labels_data.MaskDatabase(:,idx_vector)),...
        labels_data.height,labels_data.width,labels_data.Zs).*field;
%     label_cell{field,1} = reshape(full(labels_data.MaskDatabase(:,idx_vector)),...
%         labels_data.height,labels_data.width,labels_data.Zs);
    label_cell{field,2} = labels_data.MaskDatabaseNames{idx_vector};
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
%% Reformat the data

fishave_cell = reformat_fish(data,fish_ori,num_files,im_info_cell,aff_cell,ref_info,'seeds');
%% Determine the region for each seed from the registration

% allocate memory for the registred anatomy
registered_anatomy = cell(num_files,1);
label_copy = zeros(size(label_stack));
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
end

% turn the anatomy into a vector
registered_anatomy = cat(1,registered_anatomy{:});
registered_anatomy = ones(size(registered_anatomy));
% % plot stack
% figure
% % imagesc(sum(label_copy,3))
% for el = 1:15
%     subplot(3,5,el)
%     imagesc(cat(3,label_stack(:,:,el*5+63)>0,label_copy(:,:,el*5+63)>0,label_stack(:,:,el*5+63)>0))
% %     imagesc(label_copy(:,:,el*5+63)>0)
% end

figure
imagesc(cat(3,sum(label_stack,3)>0,sum(label_copy,3)>0,sum(label_stack,3)>0))
axis equal
%% Get the coordinate vector for each map

% combine the coordinate vectors
coord = round(cat(1,fishave_cell{:,1}))+coordinate_conversion;
% take the negative seeds out (also in the index)

% get a vector to keep the seeds only within the volume
selection_vector = coord(:,1)<full_ref_dim(1) & coord(:,2)<full_ref_dim(2) &...
       coord(:,3)<full_ref_dim(3);
selection_vector = selection_vector & sum(coord<1,2)==0;

% get the seed indexes corresponding to each pixel
indexes = cat(1,fishave_cell{:,2});
indexes = indexes(selection_vector);
coord = coord(selection_vector,:);
registered_anatomy = registered_anatomy(selection_vector);

% turn the coordinates into linear indexes
coord = sub2ind(full_ref_dim,coord(:,1),coord(:,2),coord(:,3));
%% Plot maps of the gains
close all
% full_ref_dim = ref_dim;
% full_ref_stack = ref_stack;
% only do this if it's p17b
if contains(data.name,'p17b')
%     % get the gains
%     delta_norm = data.delta_norm;
    
    % get the reshaped activity
    delta_norm = reshape(data.conc_trace,[],data.time_num,data.stim_num);
    % take only the stimulation time
    delta_norm = delta_norm(:,21:60,:);
    % take the absolute average
    delta_norm = squeeze(mean(abs(delta_norm),2));
    % get the max gain for each seed
    [max_gain,max_idx] = max(delta_norm,[],2);
    
    % allocate memory to store each map
    gain_maps = cell(4,1);
    for color = 1:4
        gain_maps{color} = zeros(full_ref_dim);
    end
    % run through all the seeds
    for seeds = 1:size(delta_norm,1)
        index_vector = coord(indexes==seeds&registered_anatomy>0);
%         % accumulate the gain for the pixels of each seed
%         gain_maps{max_idx(seeds)}(index_vector) = ...
%             gain_maps{max_idx(seeds)}(index_vector) + max_gain(seeds);
        % for each color
        for color = 1:4
            gain_maps{color}(index_vector) = ...
                gain_maps{color}(index_vector) + abs(delta_norm(seeds,color));
%             prc = prctile(gain_maps{color}(:),90);
%             gain_maps{color}(gain_maps{color}>prc) = prc;
        end
        
    end
    % allocate memory to save the projections
    color_cell = cell(4,1);
    % generate the maps
    for color = 1:4
        % allocate memory for the stack
%         target_stack = zeros(full_ref_dim(1),full_ref_dim(2),full_ref_dim(3),3);
%         target_stack = normr_1(repmat(full_ref_stack,1,1,1,3),1);
        curr_stack = normr_1(full_ref_stack,1);
        curr_stack(gain_maps{color}~=0) = normr_1(gain_maps{color}(gain_maps{color}~=0),1).*5;
        curr_blank = curr_stack;
        curr_blank(gain_maps{color}~=0) = 0;
        if color == 1
            target_stack = cat(4,curr_stack,curr_blank,curr_blank);
%             target_stack = cat(4,curr_blank,1-curr_stack,1-curr_stack);
        elseif color == 2
            target_stack = cat(4,curr_blank,curr_stack,curr_blank);
        elseif color == 3
            target_stack = cat(4,curr_blank,curr_blank,curr_stack);
        else
            target_stack = cat(4,curr_stack,curr_blank,curr_stack);
        end
        
        % normalize again cause of the sum
        target_stack = normr_1(target_stack,1);
        
        % take the anterior half
        half_stack = target_stack(140:760,:,:,:);
        
%         if contains(data.name,{'Syn','syn'})
%             max_projection = permute(max(half_stack,[],2),[1 3 2 4]);
%         else
            % produce a max intensity projection
            max_projection = max(half_stack,[],3);
          color_cell{color} = max_projection;
%         end
%         % turn zeros into ones
%         black_points = sum(max_projection,4)==0;
%         % for all three channels
%         for chan = 1:3
%             curr_slice = max_projection(:,:,1,chan);
%             curr_slice(black_points) = 1;
%             max_projection(:,:,1,chan) = curr_slice;
%         end
        % equalize histogram
%         max_projection = histeq(max_projection);
        % assemble the path for the output stack
%         save_path = fullfile(fig_path,strjoin({'MaxProjGain',data.name,'Color',num2str(color),'.tif'},'_'));
%         save_stack(save_path,max_projection,full_ref_info)

    end
    %% Save the projections

    % for all the clusters
    for color = 1:4
        close all
        h = fig('units','centimeters','height',4,'width',4);

        %     subplot(round(sqrt(clu_num)),ceil(sqrt(clu_num)),clu)
        %     subplot('Position',[(j-1)*1/numRecsAcross (numRecsDown-i)*1/numRecsDown 1/numRecsAcross 1/numRecsDown])
        % create the image
        I = squeeze(color_cell{color});
        %     I = histeq(I);
        I = imadjust(I,[0 0.3]);
        imagesc(I)
%         brighten(1)
        set(gca,'XTick',[],'YTick',[])
        axis square
        box off
        % assemble the save path
%         save_path = fullfile(fig_path,'clusters',strjoin({'Cluster',data.name,'Number',num2str(clu),'.tif'},'_'));
        save_path = fullfile(fig_path,strjoin({'MaxProjGain',data.name,'Color',num2str(color),'.tif'},'_'));

        % save it
        saveas(h,save_path,'tif')
    end
    error('stop')
    %% save a stack with all 4 colors
    % allocate memory for the stack
    % normalize the full ref stack
    norm_stack = normr_1(full_ref_stack,1);
    r_channel = norm_stack;
    g_channel = norm_stack;
    b_channel = norm_stack;
    a_channel = norm_stack;
    % normalize the gain
    norm_gain = normr_1(max_gain,1);
    % run through all the seeds
    for seeds = 1:size(delta_norm,1)
        index_vector = coord(indexes==seeds&registered_anatomy>0);
        switch max_idx(seeds)
            case 1
                r_channel(index_vector) = 1;
            case 2
                g_channel(index_vector) = 1;
            case 3
                b_channel(index_vector) = 1;
            case 4
                r_channel(index_vector) = 1;
                b_channel(index_vector) = 1;
        end
    end
    
    % put the final stack together
    gain_stack = cat(4,r_channel,g_channel,b_channel);
    % assemble the save path
    save_path = fullfile(fig_path,strjoin({'GainMap',data.name,'AllColors','.tif'},'_'));
    % save it
%     save_stack(save_path,gain_stack,full_ref_info)
    %% Plot max projections of the gains
    
    % take the anterior half
    half_stack = gain_stack(190:700,:,:,:);

    if contains(data.name,{'Syn','syn'})
        max_projection = permute(max(half_stack,[],2),[1 3 2 4]);
    else
        % produce a max intensity projection
        max_projection = max(half_stack,[],3);
    end

%     % save it
%     % assemble the save path
%     save_path = fullfile(fig_path,strjoin({'MaxProjGain',data.name,'.tif'},'_'));
%     % save it
%     save_stack(save_path,max_projection)
end
%% Plot clusters
close all

% get the clusters
idx_clu = data.idx_clu;
% get the number of clusters
clu_num = data.clu_num;

% get a color map for the clusters
cluster_color = lines(clu_num);

% allocate memory for the stack
% normalize the full ref stack
norm_stack = normr_1(full_ref_stack,1);
r_channel = norm_stack;
g_channel = norm_stack;
b_channel = norm_stack;
a_channel = norm_stack.*0.1;
% for all the seeds
for seeds = 1:size(data.conc_trace,1)
    % get the seed cluster
    curr_cluster = idx_clu(seeds);
    % if it's zero, skip
    if curr_cluster == 0
        continue
    end
    % get the color for the seed
    seed_color = cluster_color(curr_cluster,:);
    % get the index vector
    index_vector = coord(indexes==seeds&registered_anatomy>0);
    % color the corresponding pixels
    r_channel(index_vector) = seed_color(1);
    g_channel(index_vector) = seed_color(2);
    b_channel(index_vector) = seed_color(3);
end

% put the final stack together
full_stack = cat(4,r_channel,g_channel,b_channel);
% assemble the save path
save_path = fullfile(fig_path,strjoin({'ClusterMap',data.name,'.tif'},'_'));
% save it
% save_stack(save_path,full_stack,full_ref_info)
%% Plot max intensity projections of the clusters

% take the anterior half
half_stack = full_stack(190:700,:,:,:);

if contains(data.name,{'Syn','syn'})
    max_projection = permute(max(half_stack,[],2),[1 3 2 4]);
else
    % produce a max intensity projection
    max_projection = max(half_stack,[],3);
end

% save it
% assemble the save path
save_path = fullfile(fig_path,strjoin({'MaxProj',data.name,'.tif'},'_'));
% save it
save_stack(save_path,max_projection)
%% Plot the clusters one by one
close all
% get the clusters
idx_clu = data.idx_clu;
% get the number of clusters
clu_num = data.clu_num;

% get a color map for the clusters
cluster_color = distinguishable_colors(clu_num,[0 0 0;1 1 1]);
% normalize the full ref stack
norm_stack = cat(4,repmat(normr_1(full_ref_stack,1),1,1,1,3));
% allocate memory for the max projections
cluster_cell = cell(clu_num,1);
% for all the clusters
for clu = 1:clu_num
    fprintf(strjoin({'Current cluster',num2str(clu),'\r\n'},'_'))
    % allocate memory for the stack
    temp_stack = norm_stack;
%     r_channel = norm_stack;
%     g_channel = norm_stack;
%     b_channel = norm_stack;
%     a_channel = norm_stack.*0.1;
%     % for all the seeds
%     for seeds = 1:size(data.conc_trace,1)
    % get the seeds for this cluster
    seed_list = find(idx_clu==clu);
%         % get the seed cluster
%         curr_cluster = idx_clu(seeds);
%         % if it's zero, skip
%         if curr_cluster == 0
%             continue
%         end
%         % get the color for the seed
%         seed_color = cluster_color(curr_cluster,:);
    % get the color of the cluster
    seed_color = cluster_color(clu,:);
    
    index_vector = coord(ismember(indexes,seed_list)&(registered_anatomy>0));
    
    [x,y,z] = ind2sub(full_ref_dim,index_vector);
    index_vector_r = sub2ind([full_ref_dim,3],x,y,z,ones(size(x)));
    index_vector_g = sub2ind([full_ref_dim,3],x,y,z,ones(size(x)).*2);
    index_vector_b = sub2ind([full_ref_dim,3],x,y,z,ones(size(x)).*3);
    
    temp_stack(index_vector_r) = seed_color(1);
    temp_stack(index_vector_g) = seed_color(2);
    temp_stack(index_vector_b) = seed_color(3);
    
%     % for all the seeds in the cluster
%     for seeds = 1:length(seed_list)
%         % get the index vector
%         index_vector = coord(indexes==seed_list(seeds)&registered_anatomy>0);
%     % color the corresponding pixels
%     r_channel(index_vector) = seed_color(1);
%     g_channel(index_vector) = seed_color(2);
%     b_channel(index_vector) = seed_color(3);
%     end
    
    
    % put the final stack together
%     full_stack = cat(4,r_channel,g_channel,b_channel);
%     % assemble the save path
%     save_path = fullfile(fig_path,strjoin({'ClusterMap',data.name,'.tif'},'_'));
    % save it
    % save_stack(save_path,full_stack,full_ref_info)

    % take the anterior half
%     half_stack = full_stack(190:700,:,:,:);
    half_stack = temp_stack(140:760,:,:,:);

%     if contains(data.name,{'Syn','syn'})
%         max_projection = permute(max(half_stack,[],2),[1 3 2 4]);
%     else
        % produce a max intensity projection
        max_projection = max(half_stack,[],3);
%     end
    % store the projection
    cluster_cell{clu} = max_projection;
end
%% Plot the projections 
close all

% for all the clusters
for clu = 1:clu_num
    close all
    h = fig('units','centimeters','height',4,'width',4);

%     subplot(round(sqrt(clu_num)),ceil(sqrt(clu_num)),clu)
%     subplot('Position',[(j-1)*1/numRecsAcross (numRecsDown-i)*1/numRecsDown 1/numRecsAcross 1/numRecsDown])
    % create the image
    I = squeeze(cluster_cell{clu});
%     I = histeq(I,50000);
    I = imadjust(I,[0 0.7]);
    imagesc(I)
%     brighten(h,-0.9)
    set(gca,'XTick',[],'YTick',[])
    axis square
    box off
    % assemble the save path
    save_path = fullfile(fig_path,'clusters',strjoin({'Cluster',data.name,'Number',num2str(clu),'.tif'},'_'));
    % save it
    saveas(h,save_path,'tif')
end
%%
error('stop')
%% Prepare the average for plotting

%allocate memory for the matrix
new_volume = zeros(horzcat(ref_dim,num_files));
%initialize a fish counter
fish_c = 1;
%for all fish
for files  = 1%:num_files
    %round the new coordinates
%     curr_seeds = round(fishave_cell{fish_c,1});
    curr_seeds = round(cat(1,fishave_cell{:,1}));
    
    %clip the seeds that fall outside the volume in z, and show a flag
    over_seeds = curr_seeds(:,3)>ref_dim(3) | (curr_seeds(:,3)<1) | ...
        curr_seeds(:,1)>ref_dim(1) | curr_seeds(:,1)<1 |...
        curr_seeds(:,2)>ref_dim(2) | curr_seeds(:,2)<1;
    %if there were seeds outside the boundaries, exclude them
    if sum(over_seeds) > 0
        curr_seeds = curr_seeds(~over_seeds,:);
        fprintf(strcat('file:',num2str(files),'\r\n',...
            'excluded over lim seeds:',num2str(sum(over_seeds)),'\r\n'))
    end
    %load the data in the common volume
    new_volume(sub2ind(size(new_volume),curr_seeds(:,1),curr_seeds(:,2),...
        curr_seeds(:,3),ones(length(curr_seeds),1).*fish_c)) = 1;
    %update the counter
    fish_c = fish_c + 1;
    
end
%% Plot the fish overlap with the registered average volume
close all
% %define the volume/slice to plot
% tar_slice = 1;

%for all the fish

for files = 1
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