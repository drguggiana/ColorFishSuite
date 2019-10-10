% plot clustering
clearvars
close all force
%% Load the files and define paths

%get the folder where the image files are
tar_path_all = uipickfiles('FilterSpec','E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\Stage3\*.mat');

%get the number of experiments selected
num_exp = length(tar_path_all);

%define the list of labels to sort the files
label_list = {'_clusters.mat'};

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
%% Define/load constants

%define whether to load the raw traces
raw_var = 0;
%define whether to cluster traces
clu_var = 0;

%get the number of stimuli
stim_num2 = load(name_cell{1},'stim_num2');
stim_num2 = stim_num2.stim_num2;

%get the colors
col_out = load(name_cell{1},'col_out');
col_out = col_out.col_out;

%calculate the number of time bins
temp_conc = load(name_cell{1},'conc_trace');
temp_conc = temp_conc.conc_trace;

time_num = size(temp_conc,2)/stim_num2;

clear('temp_conc')

%path for the clustering BIC results
bic_name = strcat('E:\Behavioral data\Matlab\AF_proc\Clemens_suite\Results\bic_vec\'...
    ,'GroupFile_',num2str(num_data),'_fish');
%% Load the clusters and format

%allocate memory to store all the groups of average clusters and the number
%of traces in each cluster
clu_all = cell(num_data,2);
%also allocate memory for the index vectors
idx_all = cell(num_data,1);

%for all the fish
for fish = 1:num_data
    
    %show the current fish
    fprintf(strcat('Current fish:',num2str(fish),'\r\n'))
    %load the cluster indexes for this fish
    idx_clu = load(name_cell{fish},'idx_clu');
    idx_clu = idx_clu.idx_clu;
    
    %also load the raw traces
    conc_trace = load(name_cell{fish},'conc_trace');
    conc_trace = conc_trace.conc_trace;
        
    %using the indexes, calculate the average traces
    
    %get the number of clusters in this fish
    clu_num = load(name_cell{fish},'clu_num');
    clu_num = clu_num.clu_num;
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
    %store the average in the storage cell
    clu_all{fish,1} = clu_ave;
    %and the number of traces
    clu_all{fish,2} = clu_number;
    %store the idx_clu vector
    idx_all{fish} = idx_clu;
end
%% Plot the raw clustering results
% close all
%sort the indexes and the fish IDs
[~,s_ind] = sort(idx_clu);
figure
% subplot(1,8,1:7)
imagesc(normr_1(conc_trace(s_ind,:),0))
ylabel('Trace #')
set(gca,'XTick',[],'YTick',[],'FontSize',25)
colormap(parula)
colorbar
% subplot(1,8,8)
% imagesc(fish_ori(s_ind))

%average within each cluster
%allocate memory for the averages
clu_ave = zeros(clu_num,size(conc_trace,2));
%and for the number of items in each cluster
clu_count = zeros(clu_num,1);
%for all the clusters
for clu = 1:clu_num
    %average
    clu_ave(clu,:) = mean(conc_trace(idx_clu==clu,:),1);
    %and count
    clu_count(clu) = sum(idx_clu==clu);
end

%order the traces by abundance
[sort_count,ind_count] = sort(clu_count,'descend');

% figure
% imagesc(clu_ave)
%plot them
clu_plot = clu_ave(ind_count,:);
image_stim_1(normr_1(clu_plot,0),time_num,stim_num2,col_out,0)
plot_handle = get(gcf,'Children');
set(plot_handle(2),'YTick',1:length(ind_count),'YTickLabel',mat2cell(ind_count,ones(length(ind_count),1)),...
    'TickLength',[0 0],'FontSize',25)
set(gcf,'Colormap',parula)
colorbar(plot_handle(2))
set(plot_handle(1),'FontSize',25)
% plot_pos = get(plot_handle(2),'Position');
% stim_pos = get(plot_handle(1),'Position');
% stim_pos(3) = plot_pos(3);
% set(plot_handle(1),'Units','normalized','Position',stim_pos)

figure
bar(sort_count)
set(gca,'XTick',1:length(ind_count),'XTickLabel',mat2cell(ind_count,ones(length(ind_count),1)),...
    'TickLength',[0 0],'FontSize',20)
xlabel('Cluster','FontSize',20)
ylabel('Number of traces','FontSize',20)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])


% imagesc_cluster_beh(conc_trace(s_ind,:),sort(idx_clu),bout_signal)
% set(gca,'XTick',40:80:size(conc_trace,2),'XTickLabel',stim_labels2,'XTickLabelRotation',90)

%also plot a dendrogram of the cluster averages 
clustergram(normr_1(clu_ave,0),'Colormap',parula,'Cluster','column','Symmetric',false)

figure
%for all the clusters
for clu = 1:clu_num
    plot(normr_1(clu_ave(clu,:),0)+clu*0.7)
    hold('on')
end


error('Stop here')

%% Get the anatomy of a cluster of interest
close all

plot_save = 0;
%define the target cluster
tar_clu = 4;
%get the seed indexes of the cluster
post_seed_index = find(idx_clu==tar_clu);

%also get the indexes of the fish and original seed
ori_seed = fish_ori(post_seed_index,:);

%get the number of fish with this cluster
clu_fish = unique(ori_seed(:,1));

%load the seed xy information
seed_concat = load(name_cell{1},'cat_seed_all');
seed_concat = seed_concat.cat_seed_all;
%load the z position of each seed
z_seed = load(name_cell{1},'cat_z_all');
z_seed = z_seed.cat_z_all;
%load the average stack
ave_stack = load(name_cell{1},'cat_stack_all');
ave_stack = ave_stack.cat_stack_all;

%for all the fish
for fish = clu_fish'
    
    
    seed_index = ori_seed(ori_seed(:,1)==fish,2);
    
    %get the seeds of interest
    seed_inter = seed_concat(seed_index);
    %and their zs
    seed_inter_z = z_seed(seed_index);
    
    %get the number of seeds
    seed_clunum = length(seed_index);
    %allocate memory to store the pixels and the z
    pix_cell = cell(seed_clunum,1);
    %for all the seeds of interest
    for seeds = 1:seed_clunum
        %retrieve the pixel list of the seed
        pix_cell{seeds} = seed_inter(seeds).pxlist;
    end
    
    %using the average stack, plot the seeds
    % figure
    %allocate memory to calculate a maximum intensity projection
    max_pro = zeros(size(ave_stack,1),size(ave_stack,2),length(unique(z_seed)));
    %for all the zs
    for z = 1:length(unique(z_seed))
%             subplot(round(sqrt(z_num)),ceil(sqrt(z_num)),z)
        %allocate memory for the seed locations
        temp_im = zeros(size(ave_stack,1),size(ave_stack,2));
        %if there are seeds in this z
        if ismember(z,seed_inter_z)
            figure
            %get the seeds in this z
            temp_seed = pix_cell(z==seed_inter_z);
            
            %for all the seeds in this z
            for seeds = 1:length(temp_seed)
                %put the locations of this seed in
                temp_im(temp_seed{seeds}) = 1;
            end
            %calculate a fusion of the seed image and the average image
            fused_im = imfuse(temp_im,ave_stack(:,:,z),'ColorChannels',[1 2 2]);
%             fused_im = temp_im.*ave_stack(:,:,z);
            %store the image for max projection
            max_pro(:,:,z) = temp_im;
            %plot it
            image(fused_im)
%             colormap(hsv)
%             colormap(parula)
            title(strcat('Fish: ',num2str(fish)))
        else%otherwise, just plot the average image
%                     imagesc(ave_stack(:,:,z))
            fused_im = imfuse(temp_im,ave_stack(:,:,z),'ColorChannels',[1 2 2]);
        end
        
        %if the plotting variable is on
        if plot_save == 1
            C = fused_im;
            if fish == 1 && z == 1
                [fig_name,fig_path] = uiputfile('E:\Behavioral data\Matlab\AF_proc\Figures\*.*');
            end
            fig_full = fullfile(fig_path,strcat(fig_name,'_',num2str(fish),'.tif'));
            if z ==1
                imwrite(C,fig_full,'tif','WriteMode','overwrite')
            else
                imwrite(C,fig_full,'tif','Resolution',size(C),'WriteMode','append')
            end
        end
    end
%     %calculate the max projection and plot
%     figure
%     imagesc(max(max_pro,[],3))
end
%% Plot traces
close all

%for all the clusters
for clu = [22]%clu_num
%define the target cluster
tar_clu = clu;
%get the seed indexes of the cluster
seed_index = find(idx_clu==tar_clu);
%get the traces
tar_traces = conc_trace(seed_index,:);

%plot them
image_stim_1(tar_traces,time_num,stim_num2,col_out,2)
plot_handle = get(gcf,'Children');
set(plot_handle(2),'FontSize',25)
set(plot_handle(1),'FontSize',25)
end
%% Plot the anatomy with the dominant LED colored on top of the seeds
close all

save_color = 0;

%define the save path
fig_path = 'E:\Behavioral data\Matlab\AF_proc\Analysis\Stage3_cluster\';

%create a color vector with the number of clusters
cluster_color = hsv(clu_num);

% tar_clu = 1;
% clu_erase = zeros(clu_num,3);
% clu_erase(tar_clu,:) = [1 1 1];
% cluster_color = cluster_color.*clu_erase;

%fish progress counter
fish_c = 1;
%for all the fish
for fish = 1:num_data
    
    figure
    
    %get the name of the fish
    [~,ori_name,~] = fileparts(name_cell{fish,1});
    ori_name = ori_name(1:end-9);
    %define the new file name
    fig_name = strcat(ori_name,'_colortif.tif');
    %define the path
    fig_full = fullfile(fig_path,fig_name);
    
    %load the average stack
    ave_stack = load(name_cell{fish,1},'cat_stack_all');
    ave_stack = ave_stack.cat_stack_all;
    %get the number of z sections
    z_num = size(ave_stack,3);
    %and also the height and width
    im_width = size(ave_stack,2);
    im_height = size(ave_stack,1);
    %allocate memory for the colored stack
%     r_stack = zeros(im_height,im_width,z_num);
%     g_stack = zeros(im_height,im_width,z_num);
%     b_stack = zeros(im_height,im_width,z_num);
    ave_copy = ave_stack;
    ave_copy(isnan(ave_copy)) = mode(ave_copy(:));
    ave_copy = (ave_copy-min(ave_copy(:)))./(max(ave_copy(:))-min(ave_copy(:)));
    %for all the z
    for z = 1:z_num
        ave_copy(:,:,z) = imadjust(ave_copy(:,:,z));
    end
    r_stack = ave_copy;
    g_stack = ave_copy;
    b_stack = ave_copy;

    %load the seed xy information
    seed_concat = load(name_cell{fish,1},'cat_seed_all');
    seed_concat = seed_concat.cat_seed_all;
    %load the z position of each seed
    z_seed = load(name_cell{fish,1},'cat_z_all');
    z_seed = z_seed.cat_z_all;
    %get the number of seeds
    seed_num = size(seed_concat,1);
    %create a matrix indicating the color to use for each seed
    seed_color = zeros(seed_num,3);
%     %FOURIER INPUT
%     %load the Fourier peaks from the target fish
%     four_local = four_cat(fish_c:seed_num+fish_c-1,:);
%     %get the index of the highest fourier power
%     [~,four_power] = max(four_local(:,5:8),[],2);
%     %actually write the colors on the matrix (FLIPPING GREEN AND BLUE CAUSE
%     %PROJECTOR)
%     seed_color(four_power==1|four_power==4,1) = 1;
%     seed_color(four_power==2|four_power==4,3) = 1;
%     seed_color(four_power==3,2) = 1;

    %CLUSTER INPUT
    %get the subset of cluster ids for this fish
    idx_fish = idx_clu(fish_c:seed_num+fish_c-1);
    
    %for a desired subset of clusters
    for clu = 1:clu_num
        %define the clusters to leave
        if any(clu == [4])
            continue
        end
        %eliminate the cluster in question
        idx_fish(idx_fish==clu) = 0;
    end
    %for all the seeds
    for seed = 1:seed_num
        %if it's not a 0
        if idx_fish(seed)~=0
            %fill up the color vector with the corresponding cluster colors
            seed_color(seed,:) = cluster_color(idx_fish(seed),:);
        else
            %fill up the color vector with transparent dark
            seed_color(seed,:) = [0 0 0];
        end
    end    

%     %GAIN REGRESSOR INPUT
%     %color the points based on the regressor subsets calcualted above
%     seed_color(reg_subind(:,1),1) = 1;
%     seed_color(reg_subind(:,2),2) = 1;
%     seed_color(reg_subind(:,3),3) = 1;
%     seed_color(reg_subind(:,4),1) = 1;
%     seed_color(reg_subind(:,4),2) = 1;
%     seed_color(reg_subind(:,4),3) = 1;
    
    %update the fish counter
    fish_c = fish_c + seed_num;

    %for all the seeds
    for seed = 1:seed_num
        %get the xy indices of the seed
        [y,x] = ind2sub([im_height,im_width],seed_concat(seed).pxlist);
        %and now the linear indices in 3d
        lin_ind = sub2ind([im_height,im_width,z_num],y,x,ones(length(x),1).*z_seed(seed));
        %if the pixel is dark, skip the iteration (i.e. not a seed of
        %interest)
        if all(seed_color(seed,:)==0)
            continue
        end
        %write the color layers for this seed, one channel at a time
        r_stack(lin_ind) = seed_color(seed,1);
        g_stack(lin_ind) = seed_color(seed,2);
        b_stack(lin_ind) = seed_color(seed,3);
    end
    
    rgb_stack = cat(4,r_stack,g_stack,b_stack);
    %for all the z sections
    for z = 1:z_num
        
        figure
%         subplot(round(sqrt(z_num)),ceil(sqrt(z_num)),z)
%         C = squeeze(rgb_stack(:,:,z,:));
%         imagesc(squeeze(rgb_stack(:,:,z,:)))
%         new_ave = ave_stack(:,:,z);
%         new_ave(isnan(new_ave)) = mode(new_ave(:));
        rgb_im = squeeze(rgb_stack(:,:,z,:));
%         new_ave = ave_stack(:,:,[z z z]);
%         new_ave(isnan(new_ave)) = 0;
%         new_ave(new_ave&rgb_im) = 0;
%         new_ave = histeq(new_ave);
%         C = imfuse(squeeze(rgb_stack(:,:,z,:)),new_ave,'method','blend','Scaling','independent');
%         C = imadd(squeeze(rgb_stack(:,:,z,:)),new_ave);
%         C = new_ave + rgb_im;
%         C = imfuse(ave_stack(:,:,z),squeeze(rgb_stack(:,:,z,:)),'method','blend','Scaling','joint');
        C = rgb_im;
        imagesc(C)
        set(gca,'XTick',[],'YTick',[])
        if save_color ==1 
            if z ==1
                imwrite(C,fig_full,'tif','WriteMode','overwrite')
            else
                imwrite(C,fig_full,'tif','Resolution',size(C),'WriteMode','append')
            end
        end
        set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
%         set(gca,'DataAspectRatio',[1 1 1])
    end
end
%% Plot the anatomy stack and seeds

close all

save_anato = 0;
%define the save path
fig_path = 'E:\Behavioral data\Matlab\AF_proc\Figures';
%for all of the fish
for fish  = 1:num_data
    %get the name of the fish
    [~,ori_name,~] = fileparts(name_cell{fish,1});
    ori_name = ori_name(1:end-9);
    %define the new file name
    fig_name = strcat(ori_name,'_anato.tif');
    %define the path
    fig_full = fullfile(fig_path,fig_name);
    
    %load the average stack
    ave_stack = load(name_cell{fish,2},'ave_stack');
    ave_stack = ave_stack.ave_stack;
    %get the number of z sections
    z_num = size(ave_stack,3);
    
    figure
    %normalize the entire stack for 16 bit
    plot_stack = 255.*(ave_stack-min(ave_stack(:)))./(max(ave_stack(:))-min(ave_stack(:)));
    %for all the z slices
    for z = 1:z_num
        subplot(round(sqrt(z_num)),ceil(sqrt(z_num)),z)
        C = uint8(plot_stack(:,:,z));
%         C = ave_stack(:,:,z);
        image(C);
%         C = imfuse(ave_stack(:,:,z),im_cell{z});
%         imagesc(C)
        if save_anato == 1
            if z ==1
                imwrite(C,fig_full,'tif','WriteMode','overwrite')
            else
                imwrite(C,fig_full,'tif','Resolution',size(C),'WriteMode','append')
            end
        end
    end
        set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

end