function [loaded_structure] = load_clusters(in_path)

% initialize the structure
loaded_structure = struct([]);

% if it's a cell, assume it's full paths
if iscell(in_path)
    tar_path_all = in_path;
else
    % detect whether the path has an extension or not
    [~,file_name] = fileparts(in_path);
    
    % if no extension, use the path to open uipickfiles
    if isempty(file_name)
        %get the folder where the image files are
        tar_path_all = uipickfiles('FilterSpec',in_path);
    else
        % encapsulate in a cell and use as full path
        tar_path_all = {in_path};
    end
end


%get the number of experiments selected
num_exp = length(tar_path_all);

% %define the list of labels to sort the files
% label_list = {'_clusters.mat'};

% determine the label_list based on the extensions present
[~,~,label_list] = cellfun(@(x) fileparts(x),tar_path_all,'UniformOutput',0);
label_list = unique(label_list);

%get the number of each type of data (the round is to avoid the nonscalar
%warning for the for loop)
num_data = round(num_exp./length(label_list));

%allocate memory for the different types of files
name_cell = cell(num_data,length(label_list));

%for the types of files
for f_type = 1:length(label_list)
    
    %get the coordinates of the file names
    name_map = ~cellfun(@isempty,strfind(tar_path_all,label_list{f_type}));
 
    %store them in the corresponding layer of the name cell
    name_cell(:,f_type) = tar_path_all(name_map);
end

% for all the files
for files = 1:num_data
    % get the file name
    file_name = name_cell{files,1};
    [~,file_name] = fileparts(file_name);
    file_name = file_name(1:end-9);
    %% Define/load constants

    %get the number of stimuli
    stim_num2 = load(name_cell{files,1},'stim_num2');
    stim_num2 = stim_num2.stim_num2;

    %get the colors
    col_out = load(name_cell{files,1},'col_out');
    col_out = col_out.col_out;

    %also load the raw traces
    conc_trace = load(name_cell{files,1},'conc_trace');
    conc_trace = conc_trace.conc_trace;

    %also load the raw traces
    fish_ori = load(name_cell{files,1},'fish_ori');
    fish_ori = fish_ori.fish_ori;

    time_num = size(conc_trace,2)/stim_num2;
    if contains(label_list,'_clusters')
        pcs = load(name_cell{files,1},'pcs');
        pcs = pcs.pcs;
    else
        pcs = [];
    end

    if contains(label_list,'_gains')
        delta_norm = load(name_cell{files,1},'delta_norm');
        delta_norm = delta_norm.delta_norm;
    else
        delta_norm = [];
    end
    %% Load the cluster info
    %load the cluster indexes for this fish
    idx_clu = load(name_cell{files,1},'idx_clu');
    idx_clu = idx_clu.idx_clu;

    %using the indexes, calculate the average traces

    %get the number of clusters in this fish
    % clu_num = load(name_cell{1},'clu_num');
    % clu_num = clu_num.clu_num;
    clu_num = length(unique(idx_clu(idx_clu~=0)));

    %allocate memory for the averages
    clu_ave = zeros(clu_num,size(conc_trace,2));
    %and for the trace number
    clu_number = zeros(clu_num,1);
    %for all the clusters
    for clu = 1:clu_num
        %calculate the cluster average
        clu_ave(clu,:) = mean(conc_trace(idx_clu==clu,:),1);
        %and store the number of traces going into each average
        clu_number(clu) = sum(idx_clu==clu);
    end
    %% Package the outputs

    loaded_structure(files).file_name = file_name;
    loaded_structure(files).stim_num = stim_num2;
    loaded_structure(files).col_out = col_out;
    loaded_structure(files).conc_trace = conc_trace;
    loaded_structure(files).clu_ave = clu_ave;
    loaded_structure(files).fish_ori = fish_ori;
    loaded_structure(files).time_num = time_num;
    loaded_structure(files).delta_norm = delta_norm;
    loaded_structure(files).idx_clu = idx_clu;
    loaded_structure(files).clu_number = clu_number;
    loaded_structure(files).pcs = pcs;
end