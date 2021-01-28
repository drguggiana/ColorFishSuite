%% Calculate the convolved responses to the different RGC types from Zhou et al based on my stimuli

% load the data
clearvars
close all
load('paths.mat')
addpath(genpath(paths(1).main_path))

cluster_path = paths(1).stage3_path;
fig_path = strcat(paths(1).fig_path,'Convolution\');

% define the stimulus time
stim_time = 21:60;

data = load_clusters(cluster_path);
% define the color scheme depending on the stimulus type
if contains(data(1).name,'p17')
    color_scheme = [1 0 0;0 1 0;0 0 1;1 0 1];
    cone_color_scheme = [0.5 0 0;0 0.5 0;0 0 0.5;0.5 0 0.5];
    stim_labels = {'Red','Green','Blue','UV'};
else
%     color_scheme = distinguishable_colors(6);
    color_scheme = [1 0 0;1 0 1;1 0 0;1 0 1;1 0 0;1 0 1];
    cone_color_scheme = [];
    stim_labels = [];
end
% get the number of data sets
num_data = size(data,2);
%% Define the stimulus oscillations

close all

% define the number of stimulus frames
stimulus_frames = 40;
% define the framerate of the Zhou, Bear imaging
imaging_framerate = 15.6;

% get the time vector for the stim period
time_vector = 1:1/imaging_framerate:stimulus_frames;
% create the sine
% filtered = sin(2*pi*0.125*(6:43))';
stim_sine = 0.5+0.5.*(sin(2*pi.*0.125.*(time_vector)));

% concatenate with the flat edges
full_stim = cat(2,0.5.*ones(1,20*imaging_framerate),stim_sine,0.5.*ones(1,20*imaging_framerate));

% update the time vector to include the edges
time_vector = 1:1/imaging_framerate:80;

% normalize to 1-256
% full_stim = round(full_stim.*255+1);

figure
plot(full_stim)
%% Load the projector power

%define the number of power levels
pow_num = 256;

% define the number of LEDs
led_num = 4;
%define the LED colors for plots
led_colors = 'rbmk';

%define the path to the power files
p_power = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\20160729_powermeter';

%get the files in the selected path
p_list = dir(p_power);
%eliminate the . and ..
p_list = p_list(3:end);
%get the names from the cell array
p_names = {p_list(:).name}';
%for each color, extract the corresponding file name
p_red = p_names{contains(p_names,'red')};
p_green = p_names{contains(p_names,'green')};
p_blue = p_names{contains(p_names,'blue')};
p_uv = p_names{contains(p_names,'UV')};

%define a cell array with the spectra file names
power_files = {p_red,p_green,p_blue,p_uv};

%allocate memory for the power measurements
power_read = cell(led_num,1);

%allocate memory for the power fits
power_fits = cell(led_num,1);

%allocate memory to store the values per intensity level
power_vals = zeros(256,led_num);

%for all the spectra
for led = 1:led_num
    %open the target power file
    FID = fopen(fullfile(p_power,power_files{led}),'r');
    %read the text from the file
    power_raw = textscan(FID,'%s %s %s %f %s',...
        'HeaderLines',1,'MultipleDelimsAsOne',1);
    %close the file
    fclose(FID);
    %load the actual measurements on the allocated cell
    power_read{led} = power_raw{4};
    %low pass the curve
    wind_size = 10;
    power_read{led} = filtfilt(ones(1,wind_size)/wind_size,1,power_read{led});
    
    %fit the data to a polynomial

    %define the vector with the x values to fit
    x_vec = (255:-256/length(power_read{led}):0)';
    power_fits{led} = fit(x_vec,power_read{led}...
    ,'smoothingspline');%,'exclude',[1:5,end-3:end]);
    
    %for all 256 intensity values
    for pow = 1:pow_num
        %evaluate the fit at the desired values
        power_vals(pow,led) = feval(power_fits{led},pow-1);
    end
end
% normalize the power
power_vals = normr_1(power_vals,1);
%% Load the Zhou et al. data

% load the Zhou data
ref_data = load(paths.reference_path);
%% Get the dorsal retina kernels
% isolate the kernels for the color channels
kernel_matrix = cat(3,ref_data.AK_R_Mat,ref_data.AK_G_Mat,ref_data.AK_B_Mat,ref_data.AK_UV_Mat);

% get only the dorsal retina ones (ventral FOV)
% dorsal_bool = ref_data.Pos_Mat(:,1)>=1.5&ref_data.Pos_Mat(:,1)<=2.5;
dorsal_bool = ref_data.Pos_Mat(:,1)>=1&ref_data.Pos_Mat(:,1)<=3;
% dorsal_bool = ones(size(ref_data.Pos_Mat(:,1),1),1)==1;
kernel_matrix = kernel_matrix(dorsal_bool,:,:);
%% Convolve the stimulus with the kernels

% allocate memory for the convolved ROIs
convolved_rois = zeros(size(kernel_matrix,1),size(full_stim,2),size(kernel_matrix,3));


stim_power = [0.6397 1 0.9397 0.1727];
% get the number of ROIs
roi_number = size(kernel_matrix,1);
% for all the ROIs
for roi = 1:roi_number
    % for all the stimuli
    for stim = 1:data(1).stim_num
        convolved_rois(roi,:,stim) = conv(full_stim.*stim_power(stim),kernel_matrix(roi,:,stim),'same');
    end
end
%% Downsample the traces

% get the time points of the original stimulation (basically the frames)
original_times = 1:80;

% allocate memory for the downsampled data
downsampled_traces = zeros(roi_number,length(original_times),data(1).stim_num);

% for all the stimuli
for stim = 1:data(1).stim_num
    
    % interpolate the convolved ROIs
    downsampled_traces(:,:,stim) = interp1(time_vector,convolved_rois(:,:,stim)',original_times)';
end

% reshape for sPCA
downsampled_traces = reshape(downsampled_traces,roi_number,[]);
%% Plot the convolved, downsampled traces

close all

figure

imagesc(sort_traces(normr_1(downsampled_traces,0)))
colormap(magma)
colorbar
%% Get the traces in clustering-ready form for allocation with the GMM

% save the full traces
downsample_full = downsampled_traces;
% get only the stim period
downsampled_traces = data(1).conc_trace;
downsampled_traces = reshape(downsampled_traces,[],data(1).time_num,data(1).stim_num);
downsampled_traces = reshape(downsampled_traces(:,stim_time,:),[],length(stim_time)*data(1).stim_num);

%define the sPCA parameters to use
bounds_top = 1:length(stim_time):size(downsampled_traces,2);
bounds_bottom = [bounds_top(2:end)-1,size(downsampled_traces,2)];
bounds = [bounds_top;bounds_bottom];
K = ones(1,data(1).stim_num).*4;
t_bins = ones(data(1).stim_num,1).*10;
pca_vec = ones(data(1).stim_num,1).*1;

%define the vector of cluster numbers to try
clu_vec = [];

replicates = 20;

[~,~,~,~,~,f_data] = sPCA_GMM_cluster_Color(downsampled_traces,bounds...
    ,K,t_bins,pca_vec,[],clu_vec,replicates,[],2);

% use the GMM to allocate the Zhou traces
convolved_clusters = cluster(data(1).GMModel,f_data);
%% Infer the valid clusters by comparing original and cleaned up idx

% get only the stimulus period
conc_trace = data.conc_trace;
conc_trace = reshape(conc_trace,[],data(1).time_num,data(1).stim_num);
conc_trace = reshape(conc_trace(:,stim_time,:),[],length(stim_time)*data(1).stim_num);

[~,~,~,~,~,f_data] = sPCA_GMM_cluster_Color(conc_trace,bounds...
    ,K,t_bins,pca_vec,[],clu_vec,replicates,[],2);

% get the original idx
original_idx = cluster(data.GMModel,f_data);
% get the final idx
final_idx = data.idx_clu;
% get a vector with the cluster number
original_clunum = data.GMModel.NumComponents;

% allocate memory for the LUT
idx_LUT = zeros(original_clunum,2);
% for all the original clusters
for clu = 1:original_clunum
    idx_LUT(clu,1) = clu;
    idx_LUT(clu,2) = mode(final_idx(original_idx==clu));
end
%% Correct the idx for the kernels

for roi = 1:roi_number
    convolved_clusters(roi) = idx_LUT(convolved_clusters(roi)==idx_LUT(:,1),2);
end
%% Load the Zhou clusters

kernels_data = load(paths.kernel_cluster_path);
kernels_idx = kernels_data.clusterOpt_BIC;
% get the valid clusters
kernels_valid = kernels_data.TrueClusIndex_Vec; 
% get only the dorsal retina indexes
kernels_idx = kernels_idx(dorsal_bool);
% get the number of clusters
kernel_clunum = length(kernels_valid);
%% Compare the clusters allocation between Zhou and the present data

close all

% assemble a comparison matrix

% allocate memory for the matrix
cluster_matrix = zeros(data(1).clu_num,kernel_clunum);

% for all the traces
for roi = 1:roi_number
    % if the Zhou cluster number is not there, skip
    if any(kernels_valid,kernels_idx(roi)) == 0
        continue
    end
    % get the row coordinate (local clusters)
    x = convolved_clusters(roi);
    % get the Zhou cluster number
    y = find(kernels_valid==kernels_idx(roi));
    % add the position of the cluster as a counter
    cluster_matrix(x,y) = cluster_matrix(x,y) + 1;
end

% plot the matrix
imagesc(sortrows(cluster_matrix,'descend'))
cmap = magma(256);
cmap(1,:) = [1 1 1];
colormap(cmap)
% xlabel('Zhou et al. clusters')
% ylabel('Convolved clusters')
set(gca,'TickLength',[0 0],'FontSize',15)
axis equal
axis tight
set(gcf,'Color','w')

% create the settings
fig_set = struct([]);

fig_set(1).fig_path = fig_path;
fig_set(1).fig_name = strjoin({'Convolution_clusters',data(1).name,'.eps'},'_');
fig_set(1).fig_size = 4;
% fig_set(1).fig_size = [5 4.7];
% fig_set(2).fig_size = [2 4.7];
% fig_set(3).fig_size = [5 2];
fig_set(1).colorbar = 1;
fig_set(1).colorbar_label = '# of ROIs';
fig_set(1).box = 'on';
fig_set(1).cmap = cmap;

h = style_figure(gcf,fig_set);
%% Get the gains

data_copy = data;
data_copy(1).conc_trace = downsample_full;
% get the gains
[delta_norm, qual_res, cross_res] = gain_analysis(data_copy(1),stim_time);
%% Get the types

% get the 10th percentile
zero_threshold = prctile(abs(delta_norm),10,1);
% zero_threshold = prctile(abs(delta_norm),5,1);
% zero the values below a threshold
delta_norm(abs(delta_norm)<zero_threshold&abs(delta_norm)>0) = 0;
% turn negatives into -1 and positives into 1
delta_norm(delta_norm>0) = 1;
delta_norm(delta_norm<0) = -1;

% quantify the occurrence of each pattern
[pattern,ia,ic_conv] = unique(delta_norm,'rows');

% get the number of patterns
pattern_num = length(ia);

% allocate vector for the number
pattern_counts_conv = zeros(pattern_num,1);
% count the occurrences
% for all the patterns
for pat = 1:pattern_num
    pattern_counts_conv(pat) = sum(ic_conv==pat);
end

% sort by abundance
[pattern_counts_conv,sort_idx_conv] = sort(pattern_counts_conv,'descend');

pattern = pattern(sort_idx_conv,:);

% allocate memory for the colors
pattern_full_conv = zeros(size(pattern,1),4,3);
% transform the indexes into colors
for channel = 1:3
    pattern_full_conv(pattern(:,channel)==1,channel,channel) = 1;
    pattern_full_conv(pattern(:,channel)==0,channel,:) = 1;
    if channel == 1
        pattern_full_conv(pattern(:,4)==1,4,[1 3]) = 1;
        pattern_full_conv(pattern(:,4)==0,4,:) = 1;
    end
    
end

% % store the matrix
% type_cell{datas,1} = pattern;
% type_cell{datas,2} = pattern_counts./sum(pattern_counts);
% type_cell{datas,3} = pattern_full;

% % eliminate the patterns with only 1 instance
% elim_vector = pattern_counts_conv<2;
% pattern_counts_conv = pattern_counts_conv(~elim_vector);
% pattern_full_conv = pattern_full_conv(~elim_vector,:,:);
%% Get the types with the original traces

% calculate the 0 kernels based on the 10SD criterion used in the ref
% calculate the SD
sd_matrix = squeeze(std(kernel_matrix(:,1:150,:),0,2));

% classify in on and off

% get the maxima and minima
[max_val,max_idx] = max(kernel_matrix,[],2);
[min_val,min_idx] = min(kernel_matrix,[],2);

% calculate the null kernels
null_kernels = squeeze(max_val-min_val)<(5.*sd_matrix);
% null_kernels = sd_matrix<10;

% allocate memory for the allocation
on_off_matrix = zeros(size(max_idx,1),size(max_idx,3));
on_off_matrix(squeeze(max_idx)>squeeze(min_idx)) = 1;
on_off_matrix(squeeze(max_idx)<squeeze(min_idx)) = -1;
on_off_matrix(null_kernels) = 0;

[pattern_ref,ia,ic] = unique(on_off_matrix,'rows');
        
% get the number of patterns
pattern_num = length(ia);

% allocate vector for the number
pattern_counts = zeros(pattern_num,1);
% count the occurrences
% for all the patterns
for pat = 1:pattern_num
    pattern_counts(pat) = sum(ic==pat);
end

% sort by abundance
[pattern_counts,sort_idx] = sort(pattern_counts,'descend');
pattern_ref = pattern_ref(sort_idx,:);

% allocate memory for the colors
pattern_full = zeros(size(pattern_ref,1),4,3);
% transform the indexes into colors
for channel = 1:3
    pattern_full(pattern_ref(:,channel)==1,channel,channel) = 1;
    pattern_full(pattern_ref(:,channel)==0,channel,:) = 1;
    if channel == 1
        pattern_full(pattern_ref(:,4)==1,4,[1 3]) = 1;
        pattern_full(pattern_ref(:,4)==0,4,:) = 1;
    end
    
end
%% Match the types based on Zhou et al

% determine the total number of patterns
total_number = max(size(pattern_full_conv,1),size(pattern_full,1));

% determine th unique patterns in my data if any
[unique_target,ia] = setdiff(pattern,pattern_ref,'stable','rows');

% get the color patterns
total_patterns = vertcat(pattern_full,pattern_full_conv(ia,:,:));
total_list = vertcat(pattern_ref,pattern(ia,:));

% get the counts
ref_counts = vertcat(pattern_counts,zeros(size(unique_target,1),1));
conv_counts = zeros(size(ref_counts));
% allocate memory also for the indexes
new_ref_idx = zeros(size(ic));
new_conv_idx = zeros(size(ic_conv));
% for all the types
for types = 1:size(conv_counts,1)
    
    % get the vector for this type
    target_type = sum(on_off_matrix==total_list(types,:),2)==4;
    new_ref_idx(target_type) = types;
    
    target_type_conv = sum(delta_norm==total_list(types,:),2)==4;
    new_conv_idx(target_type_conv) = types;
    
    % find the matching pattern ad get the counts
    % if the pattern is not present, skip (leave 0)
    target = sum(pattern==total_list(types,:),2)==4;
    if sum(target)==0
        continue
    end
    conv_counts(types) = pattern_counts_conv(target);
end
%% Compare types
close all

% assemble the comparison matrix

% allocate memory for the matrix
% comparison_matrix = zeros(size(pattern_full_conv,1),size(pattern_full,1));
comparison_matrix = zeros(size(total_list,1));

% for all the traces
for roi = 1:roi_number
    % get the row and column
%     row = find(sort_idx_conv==ic_conv(roi));
%     col = find(sort_idx==ic(roi));
    row = new_conv_idx(roi);
    col = new_ref_idx(roi);
    comparison_matrix(row,col) = comparison_matrix(row,col) + 1;
end
figure
subplot('Position',[0 0.2 0.2 0.8])
% image(permute(pattern_full_conv,[1 2 3]))
image(permute(total_patterns,[1 2 3]))
set(gca,'TickLength',[0 0],'XTick',[],'YTick',[])

subplot('Position',[0.2 0.2 0.8 0.8])
imagesc(((comparison_matrix)))
cmap = magma(256);
cmap(1,:) = [1 1 1];
colormap(gca,cmap)
set(gca,'TickLength',[0 0],'XTick',[],'YTick',[])
subplot('Position',[0.2 0 0.8 0.2])
% image(permute(pattern_full,[2 1 3]))
image(permute(total_patterns,[2 1 3]))
set(gca,'TickLength',[0 0],'XTick',[],'YTick',[])

set(gcf,'Color','w')


% create the settings
fig_set = struct([]);

fig_set(1).fig_path = fig_path;
fig_set(1).fig_name = strjoin({'Convolution_types',data(1).name,'.eps'},'_');
fig_set(1).fig_size = [5 4.8];
fig_set(2).fig_size = [2 4.8];
fig_set(3).fig_size = [5 2];
% fig_set(3).colorbar = 1;
% fig_set(3).colorbar_label = '# of ROIs';
fig_set(1).box = 'on';
fig_set(3).cmap = cmap;

h = style_figure(gcf,fig_set);

% generate an independent colorbar

figure
imagesc((1:256)')
set(gca,'XTick',[],'YTick',[])

colormap(cmap)
% create the settings
fig_set = struct([]);

fig_set(1).fig_path = fig_path;
fig_set(1).fig_name = strjoin({'Convolution_types_cbar',data(1).name,'.eps'},'_');
fig_set(1).fig_size = [1 10];
fig_set(1).box = 'on';

h = style_figure(gcf,fig_set);
