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
% define the stimulus labels
% get the dataset name
stim_name = data.name;
%scan for the p17b
if contains(stim_name, 'p17b')
    %if it's p17b
    stim_labels = {'Red','Green','Blue','UV'};
else %if it's p6p8 instead
    %define the stim labels (for p6p8 data)
    stim_labels = {'Red CK','UV CK','Red GR','UV GR','Red FL','UV FL'};
end
% get the stimulus number
stim_num = data.stim_num;
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
%% Plot maps of the max response
close all
% full_ref_dim = ref_dim;
% full_ref_stack = ref_stack;
% only do this if it's p17b
% if contains(data.name,'p17b')
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
% get the number of stimuli
stim_num = data.stim_num;
% allocate memory to store each map
gain_maps = cell(stim_num,1);
for color = 1:stim_num
    gain_maps{color} = zeros(full_ref_dim);
end
% run through all the seeds
for seeds = 1:size(delta_norm,1)
    index_vector = coord(indexes==seeds&registered_anatomy>0);
    %         % accumulate the gain for the pixels of each seed
    %         gain_maps{max_idx(seeds)}(index_vector) = ...
    %             gain_maps{max_idx(seeds)}(index_vector) + max_gain(seeds);
    % for each color
    for color = 1:stim_num
        gain_maps{color}(index_vector) = ...
            gain_maps{color}(index_vector) + abs(delta_norm(seeds,color));
        %             prc = prctile(gain_maps{color}(:),90);
        %             gain_maps{color}(gain_maps{color}>prc) = prc;
    end
    
end
% allocate memory to save the projections
color_cell = cell(stim_num,1);
% generate the maps
for color = 1:stim_num
    % allocate memory for the stack
    %         target_stack = zeros(full_ref_dim(1),full_ref_dim(2),full_ref_dim(3),3);
    %         target_stack = normr_1(repmat(full_ref_stack,1,1,1,3),1);
    curr_stack = normr_1(full_ref_stack,1);
    curr_stack(gain_maps{color}~=0) = normr_1(gain_maps{color}(gain_maps{color}~=0),1).*5;
    curr_blank = curr_stack;
    curr_blank(gain_maps{color}~=0) = 0;
    
    if contains(data.name,'p17b')
        switch color
            case 1
                target_stack = cat(4,curr_stack,curr_blank,curr_blank);
                %             target_stack = cat(4,curr_blank,1-curr_stack,1-curr_stack);
            case 2
                target_stack = cat(4,curr_blank,curr_stack,curr_blank);
            case 3
                target_stack = cat(4,curr_blank,curr_blank,curr_stack);
            case 4
                target_stack = cat(4,curr_stack,curr_blank,curr_stack);
        end
    else
        switch color
            case {1,3,5}
                target_stack = cat(4,curr_stack,curr_blank,curr_blank);
            case {2,4,6}
                target_stack = cat(4,curr_stack,curr_blank,curr_stack);
        end
    end
    % normalize again cause of the sum
    target_stack = normr_1(target_stack,1);
    
    % take the anterior half
    if contains(data.name,'p17b')
        half_stack = target_stack(140:760,:,:,:);
    else
        half_stack = target_stack(300:700,110:511,:,:);
    end
    
    % produce a max intensity projection
    max_projection = max(half_stack,[],3);
    color_cell{color} = max_projection;
    
    
end
%% Save the max response projections

% for all the clusters
for color = 1:stim_num
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
    axis equal
    box off
    % assemble the save path
    %         save_path = fullfile(fig_path,'clusters',strjoin({'Cluster',data.name,'Number',num2str(clu),'.tif'},'_'));
    save_path = fullfile(fig_path,strjoin({'MaxProjGain',data.name,'Color',num2str(color),'.tif'},'_'));
    
    % save it
    saveas(h,save_path,'tif')
end
%% OFF Save a stack with all 4 colors
% % allocate memory for the stack
% % normalize the full ref stack
% norm_stack = normr_1(full_ref_stack,1);
% r_channel = norm_stack;
% g_channel = norm_stack;
% b_channel = norm_stack;
% a_channel = norm_stack;
% % normalize the gain
% norm_gain = normr_1(max_gain,1);
% % run through all the seeds
% for seeds = 1:size(delta_norm,1)
%     index_vector = coord(indexes==seeds&registered_anatomy>0);
%     switch max_idx(seeds)
%         case 1
%             r_channel(index_vector) = 1;
%         case 2
%             g_channel(index_vector) = 1;
%         case 3
%             b_channel(index_vector) = 1;
%         case 4
%             r_channel(index_vector) = 1;
%             b_channel(index_vector) = 1;
%     end
% end
% 
% % put the final stack together
% gain_stack = cat(4,r_channel,g_channel,b_channel);
% % assemble the save path
% save_path = fullfile(fig_path,strjoin({'GainMap',data.name,'AllColors','.tif'},'_'));
% % save it
% %     save_stack(save_path,gain_stack,full_ref_info)
% %% Plot max projections of the gains
% 
% % take the anterior half
% half_stack = gain_stack(190:700,:,:,:);
% 
% if contains(data.name,{'Syn','syn'})
%     max_projection = permute(max(half_stack,[],2),[1 3 2 4]);
% else
%     % produce a max intensity projection
%     max_projection = max(half_stack,[],3);
% end
% 
% %     % save it
% %     % assemble the save path
% %     save_path = fullfile(fig_path,strjoin({'MaxProjGain',data.name,'.tif'},'_'));
% %     % save it
% %     save_stack(save_path,max_projection)
% % end
%% OFF Plot clusters all in the same map
% close all
% 
% % get the clusters
% idx_clu = data.idx_clu;
% % get the number of clusters
% clu_num = data.clu_num;
% 
% % get a color map for the clusters
% cluster_color = distinguishable_colors(clu_num);
% 
% % allocate memory for the stack
% % normalize the full ref stack
% norm_stack = normr_1(full_ref_stack,1);
% r_channel = norm_stack;
% g_channel = norm_stack;
% b_channel = norm_stack;
% a_channel = norm_stack.*0.1;
% % for all the seeds
% for seeds = 1:size(data.conc_trace,1)
%     % get the seed cluster
%     curr_cluster = idx_clu(seeds);
%     % if it's zero, skip
%     if curr_cluster == 0
%         continue
%     end
%     % get the color for the seed
%     seed_color = cluster_color(curr_cluster,:);
%     % get the index vector
%     index_vector = coord(indexes==seeds&registered_anatomy>0);
%     % color the corresponding pixels
%     r_channel(index_vector) = seed_color(1);
%     g_channel(index_vector) = seed_color(2);
%     b_channel(index_vector) = seed_color(3);
% end
% 
% % put the final stack together
% full_stack = cat(4,r_channel,g_channel,b_channel);
% % assemble the save path
% save_path = fullfile(fig_path,strjoin({'ClusterMap',data.name,'.tif'},'_'));
% % save it
% % save_stack(save_path,full_stack,full_ref_info)
% %% Plot max intensity projections of the clusters all in the same map
% 
% % take the anterior half
% half_stack = full_stack(190:700,:,:,:);
% 
% if contains(data.name,{'Syn','syn'})
%     max_projection = permute(max(half_stack,[],2),[1 3 2 4]);
% else
%     % produce a max intensity projection
%     max_projection = max(half_stack,[],3);
% end
% 
% % save it
% % assemble the save path
% save_path = fullfile(fig_path,strjoin({'MaxProj',data.name,'.tif'},'_'));
% % save it
% save_stack(save_path,max_projection)
%% Calculate the clusters one by one
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

    % get the seeds for this cluster
    seed_list = find(idx_clu==clu);

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

    % take the anterior half (protocol dependent)
    if contains(data.name,'p17b')
        half_stack = temp_stack(140:760,:,:,:);
    else
        half_stack = temp_stack(300:700,:,:,:);
    end


    % produce a max intensity projection
    max_projection = max(half_stack,[],3);
    % store the projection
    cluster_cell{clu} = max_projection;
end
%% Plot the cluster projections one by one 
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
%% OFF Calculate the overlap between stimulus responses
% % determine the number of stimulus pairs
% stim_combinations = nchoosek(1:stim_num,2);
% number_combinations = size(stim_combinations,1);
% % determine the filtering threshold
% conc_trace = abs(reshape(data.conc_trace,[],data.time_num,stim_num));
% baseline = conc_trace(:,1:20,:);
% mean_baseline = squeeze(mean(conc_trace,2));
% cell_ave = mean(mean_baseline,1);
% cell_std = 2*std(mean_baseline,1);
% ROI_threshold = (cell_ave+cell_std);
% % concat_baseline = reshape(baseline,[],4);
% % ROI_threshold = prctile(concat_baseline,8,1);
% % allocate memory for the overlaps
% overlap_matrix = zeros(stim_num);
% % for all the combinations
% for combs = 1:number_combinations
%     % get the maps
%     map1 = abs(gain_maps{stim_combinations(combs,1)});
%     map2 = abs(gain_maps{stim_combinations(combs,2)});
%     
%     % gaussian blur them
% %     map1_gauss = imgaussfilt(map1,3);
% %     map1_gauss(map1>0) = 0;
% %     map1 = map1_gauss(:)>0;
%     map1 = map1(:)>ROI_threshold(stim_combinations(combs,1));
% %     map2_gauss = imgaussfilt(map2,3);
% %     map2_gauss(map2>0) = 0;
% %     map2 = map2_gauss(:)>0;
%     map2 = map2(:)>ROI_threshold(stim_combinations(combs,2));
%     
%     % quantify the overlap
%     overlap_matrix(stim_combinations(combs,1),stim_combinations(combs,2)) = sum(sum([map1,map2],2)>1);
%     
% end
%% Calculate the overlap between stimuli and compare to reps

% determine the filtering threshold
all_reps = abs(reshape(data.single_reps,size(data.single_reps,1),data.time_num,stim_num,[]));
% get the number of reps
rep_num = size(all_reps,4);
% allocate memory to store the maps per rep
rep_maps = cell(rep_num,1);
% fir all the reps
for reps = 1:rep_num
    % get the current rep
    conc_trace = all_reps(:,:,:,reps);
    
    % take only the stimulation time
    delta_norm = conc_trace(:,21:60,:);
    % take the absolute average
    delta_norm = squeeze(mean(abs(delta_norm),2));
    % calculate the maps
    maps_cell = property_map(coord,indexes,registered_anatomy,delta_norm,stim_num,full_ref_dim);

    baseline = conc_trace(:,1:20,:);
    mean_baseline = squeeze(mean(conc_trace,2));
    
    cell_ave = mean(mean_baseline,1);
    cell_std = 2*std(mean_baseline,1);
    ROI_threshold = (cell_ave+cell_std);
    
    % threshold them
    
    % for all the maps
    for stim = 1:stim_num
        % apply the threshold and save
        maps_cell{stim} = abs(maps_cell{stim});
        maps_cell{stim} = maps_cell{stim}(:)>ROI_threshold(stim);
    end
    % save in the rep map concatenated by stimulus
    rep_maps{reps} = cat(2,maps_cell{:});
end

% concatenate the rep maps into a single matrix
rep_maps = cat(3,rep_maps{:});
%% Calculate the overlap matrices

% determine the number of stimulus pairs
stim_combinations = nchoosek(1:stim_num,2);
number_combinations = size(stim_combinations,1);
% concat_baseline = reshape(baseline,[],4);
% ROI_threshold = prctile(concat_baseline,8,1);
% allocate memory for the overlaps
overlap_matrix = zeros(stim_num,stim_num,rep_num);
% for all the reps
for reps = 1:rep_num
    % for all the combinations
    for combs = 1:number_combinations
        % get the maps
        map1 = rep_maps(:,stim_combinations(combs,1),reps);
        map2 = rep_maps(:,stim_combinations(combs,2),reps);

        % quantify the overlap
        overlap_matrix(stim_combinations(combs,1),stim_combinations(combs,2),reps) = sum(sum([map1,map2],2)>1);

    end
end
%% Calculate the control matrix with the reps for each stimulus

% allocate memory to store the overlaps between reps
control_overlap = zeros(stim_num,1);
% get the combinations for the reps
control_combs = nchoosek(1:rep_num,2);
control_num = size(control_combs,1);
% for all the stimuli
for stim = 1:stim_num
    % get the maps for this stim
    stim_maps = squeeze(rep_maps(:,stim,:));
    % for all the combinations
    for combs = 1:control_num
        % select the corresponding maps
        map1 = stim_maps(:,control_combs(combs,1));
        map2 = stim_maps(:,control_combs(combs,2));
        % calculate the overlap
        control_overlap(stim) = control_overlap(stim) + sum(sum([map1,map2],2)>1)./control_num;
    end
end
%% Plot the matrix
close all
% define the fontsize
fontsize = 15;
% normalize the matrix by the average of the rep distances
% allocate memory for the normalized matrix
norm_overlap = mean(overlap_matrix,3);
% for all the combinations
for combs = 1:number_combinations
    % normalize
    norm_overlap(stim_combinations(combs,1),stim_combinations(combs,2)) = ...
        1-(norm_overlap(stim_combinations(combs,1),stim_combinations(combs,2))./mean(...
        control_overlap([stim_combinations(combs,1),stim_combinations(combs,2)])));
end

figure
imagesc(norm_overlap)
axis square
% set the color scale to the max number of trials per category
title(data.figure_name)
set(gca,'CLim',[0, 1])
set(gca,'TickLength',[0 0])
set(gca,'XTick',1:stim_num,'XTickLabels',stim_labels,'FontSize',fontsize,...
    'XTickLabelRotation',45)
set(gca,'YTick',1:stim_num,'YTickLabels',stim_labels,'FontSize',fontsize)
set(gca,'FontSize',20)

cba = colorbar;
set(cba,'TickLength',0)
ylabel(cba,'Dissimilarity Index')
% assemble the figure path
file_path = strjoin({'anatomicalOverlap',stim_name},'_');
% saveas(gcf, fullfile(fig_path,file_path), 'png')
print(fullfile(fig_path,file_path),'-dpng','-r600')
