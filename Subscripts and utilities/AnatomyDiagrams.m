%% Plot anatomical diagrams
clearvars
close all

% define the figure path
load('paths.mat')
addpath(genpath(paths(1).main_path))

cluster_path = paths(1).stage3_path;
fig_path = strcat(paths(1).fig_path,'Setup\');
registration_path = paths(1).registration_path;

% define which map to target
target_map = 'syn';
%% Load the reference brain

% define the ref brain depending on the file
switch target_map
    case {'syn','AF10'}
        ref_name = 'refisl2cut2.nrrd.tif';
%         full_ref_name = 'refisl2.nrrd.tif';
        full_ref_name = 'refbrain.nrrd.tif';

    %     ref_name = 'refcutblursub.nrrd.tif';
    %     full_ref_name = 'refbrain.nrrd.tif';
    case {'huc','Tectum'}
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
%% Load the labels

% define the path to the file
labels_file = fullfile(registration_path,'Labels_info','MaskDatabase.mat');
% load the labels file
labels_data = load(labels_file);

% define which dataset to load depending on the reference
switch target_map
    case 'syn'
        field_list = {'AF4','AF5','AF6','AF7','AF8','AF9','Tecum Neuropil'};
        field_names = {'AF4','AF5','AF6','AF7','AF8','AF9','AF10'};
        bar_ratio = 3;
        bar_length = 4.3;
        coordinate_conversion = [274,0,0];
        front_limit = 300;
        side_limit = 621;
%         cmap = lines(length(field_names));
        cmap = distinguishable_colors(length(field_names),{'k','w','g','r','b','m','y'});
        im_range = 0.7;

    case 'huc'
        field_list = {'Tectum Stratum','Tecum Neuropil','Pretectum','Habenula','Cerebellum'};
        field_names = {'TcN','TcP','Pt','Hb','Cb'};
        bar_ratio = 5;
        bar_length = 5.4;
        coordinate_conversion = [272,74,49];
        front_limit = 200;
        side_limit = 621;
        cmap = distinguishable_colors(length(field_names),{'k','y','w','g','r','b','m'});
        im_range = 0.7;
    case 'AF10'
        field_list = {'Tecum Neuropil'};
        field_names = {'AF10'};
        coordinate_conversion = [274,0,0];
        front_limit = 300;
        side_limit = 310;
        cmap = distinguishable_colors(7,{'k','w','g','r','b','m','y'});
        cmap = cmap(end,:);
        im_range = 0.3;
    case 'Tectum'
        field_list = {'Tectum Stratum','Tecum Neuropil'};
        field_names = {'TcN','TcP'};
        coordinate_conversion = [272,74,49];
        front_limit = 200;
        side_limit = 310;
        cmap = distinguishable_colors(length(field_names),{'k','y','w','g','r','b','m'});
        cmap = cmap(2,:);
%         cmap = magma(3);
%         cmap = cmap([2 2],:);
        im_range = 0.3;
end

% get the number of fields
field_number = length(field_list);
% allocate memory for the fields
label_cell = cell(field_number,2);
label_stack = zeros(labels_data.height,labels_data.width,labels_data.Zs);
% get the field numbers for these names
for field = 1:field_number
    % get the index with the first name match (so as to not get
    % subdivisions)
    idx_vector = find(contains(labels_data.MaskDatabaseNames,field_list{field}),1);
    % store the map and the name in the cell
    label_stack = label_stack + reshape(full(labels_data.MaskDatabase(:,idx_vector)),...
        labels_data.height,labels_data.width,labels_data.Zs).*field;
    label_cell{field,1} = reshape(full(labels_data.MaskDatabase(:,idx_vector)),...
        labels_data.height,labels_data.width,labels_data.Zs);
    label_cell{field,2} = labels_data.MaskDatabaseNames{idx_vector};
end
% if it's tectum, combine the two parts of tectum
if strcmp(target_map,'Tectum')
    label_cell = {label_cell{1,1}|label_cell{2,1},'Tectum'};
    field_number = 1;
end
%% Plot a max projection with the labels

close all

% for all projections
for proj = 1:3
    % max project the ref stack
    max_ref = 1-repmat(normr_1(squeeze(max(full_ref_stack,[],proj)),1),1,1,3);
    
    % allocate memory to store the max projections
    max_cell = cell(field_number,1);
    % allocate memory to store the colored pics
    color_cell = cell(field_number,1);
    
    % get the color map
    % cmap = distinguishable_colors(field_number+3,{'w','k'});
    % cmap = cmap(4:end,:);
    % max project and plot all the maps, also turn them into full color stacks
    for field = 1:field_number
        % calculate the projection
        max_cell{field} = squeeze(max(label_cell{field,1},[],proj));
        
        %     figure
        %     imagesc(max_cell{field})
        %     axis equal
        %     axis tight
        % create the color version
        %     color_temp = zeros(size(max_cell{field},1),size(max_cell{field},2),3);
        color_temp = 1-double(repmat(normr_1(max_cell{field},1),1,1,3));
        % for all 3 channels
        for channel = 1:3
            temp = color_temp(:,:,channel);
            temp(temp==0) = cmap(field,channel);
            color_temp(:,:,channel) = temp;
%             color_temp(:,:,channel) = color_temp(:,:,channel).*cmap(field,channel);
        end
        % store the image
        color_cell{field} = color_temp;
        %     figure
        %     imagesc(color_temp)
        %     axis equal
        %     axis tight
        
        
    end
    % add the projections with alpha blending
    max_all = cat(4,max_ref,color_cell{:});
    switch proj
        case 1
            permute_vector = [2 1 3];
            x_range = 1:138;
            y_range = 1:side_limit;
            y_dir = 'normal';
            aspect_vector = [0.798*side_limit 2*138 1];
        case 2
            permute_vector = [2 1 3];
            x_range = 1:138;
            y_range = front_limit:700;
            y_dir = 'normal';
            aspect_vector = [0.798*(700-front_limit) 2*138 1];
        case 3
            permute_vector = [1 2 3];
            x_range = front_limit:700;
            y_range = 1:side_limit;
            y_dir = 'reverse';
            aspect_vector = [side_limit 700-front_limit 1];
    end
    figure
    I = permute(mean(max_all,4),permute_vector);
    I = imadjust(I(x_range,y_range,:),[im_range 1]);
    imagesc(I)
    set(gca,'YDir',y_dir,'TickLength',[0 0])
    set(gca,'XTick',[],'YTick',[],'LineWidth',2)
%     axis equal
    pbaspect(aspect_vector)
    axis tight
    % save the figure
    file_path = strjoin({'Anatomy',target_map,num2str(proj),'.png'},'_');
    export_fig(fullfile(fig_path,file_path),'-r600')
    if proj == 3
        %%
        figure
        imagesc((1:field_number)')
%         colormap(gca,cmap)
        set(gca,'YAxisLocation','right','TickLength',[0 0],'FontSize',15,'LineWidth',2)
%         set(gca,'XTick',[],'YTick',1:field_number,'YTickLabels',field_names)
        set(gca,'XTick',[],'YTick',[])
%         set(gcf,'Color','w')
%         axis equal
        daspect([bar_ratio 1 1])
        axis tight
%         % save the figure
%         file_path = strjoin({'Anatomy',target_map,'legend','.png'},'_');
%         export_fig(fullfile(fig_path,file_path),'-r600')
                % create the settings
        fig_set = struct([]);
        
        fig_set(1).fig_path = fig_path;
        fig_set(1).fig_name = strjoin({'Anatomy',target_map,'legend','.png'},'_');
        fig_set(1).fig_size = [1 bar_length];
        fig_set(1).box = 'on';
        fig_set(1).font_size = 'medium';
        fig_set(1).cmap = cmap;
        
        h = style_figure(gcf,fig_set);
    end
end
autoArrangeFigures