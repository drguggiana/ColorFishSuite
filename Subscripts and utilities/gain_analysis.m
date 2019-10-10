% gain analysis
clearvars
close all force
%% Load clusters

[clu_ave, idx_clu, clu_number, conc_trace, stim_num2, time_num, col_out] = load_clusters();
%% Extract the Fourier power at the stimulus frequency

[four_cat, four_qual] = fourier_extraction(conc_trace,time_num,stim_num2);
%% Quantify the cone gains based on the cone isolating stimuli

close all
%gotta define it
cone_num = 4;
% define the stimulus time
stim_time = 21:60;
%allocate memory for the delta signals and delta excitations
delta_signals = zeros(size(conc_trace,1),stim_num2);
%initialize a frame counter
frame_c = 1;
%for all the stimuli
for stim = 1:stim_num2
    %         %load the traces for this stim (just the stimulus portion)
    %         stim_traces = conc_trace(:,frame_c:frame_c+time_num-1);
    %         stim_traces = stim_traces(:,stim_time);
    %         %find the delta max-min intensity swings for each trace
    %         delta_signals(:,stim) = (max(stim_traces,[],2)-min(stim_traces,[],2))./min(stim_traces,[],2);
    
    %use the Fourier assignments as the value to solve for (they are
    %calculated below, should fix the code post DAC)
    delta_signals(:,stim) = four_cat(:,stim+4);
    
    %update the frame counter
    frame_c = frame_c + time_num;
end

%append the fourier quality for clustering (both in rbgu)
delta_clu = cat(2,delta_signals,four_qual(:,5:8));
%fix the color order to rgbu
delta_signals = delta_signals(:,[1 3 2 4]);
figure
imagesc(delta_signals)
figure
%for all the cones
for cones = 1:cone_num
    subplot(2,2,cones)
    histogram(delta_signals(:,cones),100)
end
figure
errorbar(1:cone_num,mean(delta_signals,1),std(delta_signals,0,1)./sqrt(11))
%% Cluster the patterns using GMM
%define the vector of cluster numbers to try
clu_vec = [2 3 4 5 10 20 30 40 50];
[idx_clu,GMModel,clu_num] = sPCA_GMM_cluster_2(delta_clu,0,[],[],[],[],bic_name,clu_vec);

%% Determine the gain each trace gets from each cone

%the idea is that each stimulus presents a different equation with 4
%variables. Since there are 4 LEDs (and hence 4 stimuli) that are linearly
%independent then the system of equations should have a solution
close all
%define common constants
cone_num = 4;
led_num = 4;

%load the cone excitation matrix

%define the path
cone_exc_path = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\201609019_coneexc';
% %get the file
% [cexc_name,cexc_fpath] = uigetfile(cexc_path);

%if the name is known
cone_exc_fpath = cone_exc_path;
cone_exc_name = '20160925T215300_ConeExc.mat';
%load it (dimensions of cone_exc: CONE, LED, POWER)
load(fullfile(cone_exc_fpath,cone_exc_name))
%     %scale the variable to a more intuitive range
%     cone_exc = cone_exc*1e9;

%     %flip the middle columnts in cexc_mat since these are switched in the
%     %stimuli coming from the scope
%     cexc_mat = cexc_mat(:,[1 3 2 4],:);
%calculate the delta excitations for the stimuli

%% threshold the fourier components
%extract the relevant portion of the fourier components
curr_four = four_cat(:,5:8);

%     %plot the overall histogram of the components
%     figure
%     histogram(curr_four(:),100)

%make the signals in the bottom 5th percentile 0
curr_four(curr_four<prctile(curr_four(:),10)) = 0;
%replace the values in the original fourier matrix
four_cat(:,5:8) = curr_four;
%% Perform the gain analysis
%allocate memory for the delta signals and delta excitations
delta_signals = zeros(size(conc_trace,1),stim_num2);
delta_exc = zeros(cone_num,led_num);
%and for the phase shift value
xcorr_signals = zeros(size(conc_trace,1),stim_num2);
%initialize a frame counter
frame_c = 1;
%for all the stimuli
for stim = 1:stim_num2
    
    %load the color info for this stimulus
    col_info = squeeze(col_out(stim,:,:));
    %reshape to put the two sides in different dimensions
    col_info = reshape(col_info,size(col_info,1),4,2);
    %get the delta excitations
    max_col = max(col_info,[],3);
    max_col = max(max_col,[],1)+1;
    min_col = min(col_info,[],3);
    min_col = min(min_col,[],1)+1;
    
    %using that information, extract the excitations for this LED
    max_exc = squeeze(cone_exc(:,max_col>1,max(max_col)));
    min_exc = squeeze(cone_exc(:,max_col>1,min(min_col)));
    %calculate the delta exc (STIM/LED in rows, CONES in columns)
    delta_exc(stim,:) = (max_exc - min_exc);
    
    %load the traces for this stim (just the stimulus portion)
    stim_traces = conc_trace(:,frame_c:frame_c+time_num-1);
    stim_traces = stim_traces(:,stim_time);
    %         %find the delta max-min intensity swings for each trace
    %         delta_signals(:,stim) = (max(stim_traces,[],2)-min(stim_traces,[],2));
    
    % use the Fourier assignments as the value to solve for
    %         delta_signals(:,stim) = four_cat(:,stim+4);
    
    %use the full swing of the signal
    delta_signals(:,stim) = signal_swing(:,stim);
    
    %turn the stimulus to linear form
    %         lin_stim = reshape(col_info(:,:,1),320,1);
    lin_stim = col_info(:,stim,1);
    lin_stim = lin_stim(stim_time);
    %calculate the cross correlation between stimulus and trace to
    %determine the amount of phase shift with respect to the stimulus
    %for all the traces
    for traces = 1:size(conc_trace,1)
        [~,xcorr_signals(traces,stim)] = max(xcorr(stim_traces(traces,:),lin_stim));
    end
    
    %update the frame counter
    frame_c = frame_c + time_num;
end

%solve the system of equations for all the traces

%allocate memory for the results
delta_res = zeros(size(conc_trace,1),cone_num);
qual_res = zeros(size(conc_trace,1),cone_num);
cross_res = zeros(size(conc_trace,1),cone_num);

%for all the traces
for traces = 1:size(conc_trace,1)
    %calculate the gains for this trace
    delta_res(traces,:) = delta_exc\(delta_signals(traces,:))';
    qual_res(traces,:) = delta_exc\(qual_cat(traces,5:8))';
    cross_res(traces,:) = delta_exc\(xcorr_signals(traces,:))';
end

%normalize for clustering
delta_norm = delta_res./max(delta_res(:));
delta_clu = delta_norm;

%     %z score the gain values
%     delta_clu = zscore(delta_res,0,2);

%     %append the fourier quality for clustering (changing the four_qual to
%     %rgbu also
%     delta_clu = cat(2,delta_res./max(delta_res(:)),...
%         four_qual(:,:,2)./max(four_qual(:)),xcorr_signals./max(xcorr_signals(:)));
%      delta_clu = cat(2,delta_res./max(delta_res(:)),...
%         qual_res./max(qual_res(:)),xcorr_signals./max(xcorr_signals(:)));
%     delta_clu = normr_1(cat(2,round(normr_1(log(normr_1(zscore(delta_res,0,1),1)+1),1).*20),...
%         round(normr_1(log(normr_1(xcorr_signals,1)+1),1).*20)),1);
%     delta_clu = normr_1(cat(2,round(normr_1(log(normr_1(zscore(delta_res,0,1),1)+1),1).*20),...
%         xcorr_signals),1);
% %LAST USED
%     delta_clu = cat(2,normr_1(log(normr_1(zscore(delta_res,0,1),1)+1),1),...
%         normr_1(log(normr_1(xcorr_signals,1)+1),1));

%     delta_clu = normr_1(round(normr_1(log(normr_1(zscore(delta_res,0,1),1)+1),1).*20),1);

figure
imagesc(delta_clu)

figure
%for all the cones
for cones = 1:cone_num
    subplot(2,2,cones)
    histogram(delta_res(:,cones),100)
end
%% Plot the average gain for each cone type
figure
h = errorbar(1:cone_num,mean(delta_res(:,[4 3 2 1]),1),std(delta_res(:,[4 3 2 1]),0,1)./sqrt(fish_num));
set(h,'LineStyle','none','Marker','o','MarkerSize',10)
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'UV','Blue','Green','Red'},'Fontsize',20)
xlabel('Cone','Fontsize',20)
ylabel('Gain population average (a.u.)','Fontsize',20)
%     hold('on')
%     %for all the points
%     for points = 1:size(delta_res,1)
%         %generate a random x
%         randx = rand*0.4-0.2;
%         %generate a vector with the mod x
%         mod_x = (1:4) + randx;
%         %plot the points
%         plot(mod_x,delta_res(po ints,:),'*')
%     end
error('Gains')
%% Cluster the patterns using GMM
%define the vector of cluster numbers to try
%     clu_vec = [5 10 20 30 40 50 70 100];
%     clu_vec = [80 90 100 150 200];
%     clu_vec = [32 33 34 35 36 37 38 40];
clu_vec = [5 10 20 30 40];
%     clu_num = 0;
replicates = 1;
[idx_clu,GMModel,clu_num] = sPCA_GMM_cluster_4(delta_clu,[],[],[],[],bic_name,clu_vec,replicates);
%% OFF Cluster the patters using clustergram

%     cobj = clustergram(delta_clu,'Cluster','column');
%% Cluster the patterns using k-means
%define the number of clusters
clu_num = 20;
%perform the clustering
idx_clu = kmeans(delta_clu,clu_num);
%% OFF tSNE the combinations of gains
%     no_dims = 2;
%     perplex = 30;
%     in_dims = 4;
%     t_label = [];
%     mapped_g = tsne(delta_res,t_label,no_dims,in_dims,perplex);
%% Plot tSNE results
close all
figure
gscatter(mapped_g(:,1),mapped_g(:,2),idx_clu)
%% OFF redundant with corrplot Generate a correlation plot matrix
%     close all
%     [h,ax] = plotmatrix(delta_res);
%     labels = {'Red','Green','Blue','UV'};
%
%     for i = 1:4                                       % label the plots
%         xlabel(ax(4,i), labels{i})
%         ylabel(ax(i,1), labels{i})
%     end
%% Use a 3D scatter
close all
figure
scatter3(delta_res(:,1),delta_res(:,2),delta_res(:,3),40,delta_res(:,4),'filled')
xlabel('Red')
ylabel('Green')
zlabel('Blue')

colorbar
%% Plot the correlation matrix

close all

labels = {'Red','Green','Blue','UV'};
%     labels = {'Red','Green','Blue','UV','q1','q2','q3','q4','s1','s2','s3','s4'};

corrplot(delta_res,'varNames',labels)
%% Plot the gain histograms

close all
figure
labels = {'Red','Green','Blue','UV'};
%for all the stimuli
for stim = 1:stim_num2
    subplot(2,2,stim)
    histogram(delta_res(:,stim),100)
    set(gca,'FontSize',25)
    xlabel(labels{stim},'FontSize',25)
end
%% Perform a shuffle analysis with the fourier components

%get the relevant fourier assignments
curr_four = four_cat(:,5:8);

%define the number of shuffles
shuff_num = 100;

%allocate memory for the shuffle results
shuff_res = zeros(shuff_num,4);

%for all the shuffles
for shuff = 1:shuff_num
    %get the sizes of the matrix
    [M,N] = size(curr_four);
    % Preserve the row indices
    rowIndex = repmat((1:M)',[1 N]);
    % Get randomized column indices by sorting a second random array
    [~,randomizedColIndex] = sort(rand(M,N),2);
    % Need to use linear indexing to create B
    newLinearIndex = sub2ind([M,N],rowIndex,randomizedColIndex);
    B = curr_four(newLinearIndex);
    %now compute the averages and store
    shuff_res(shuff,:) = mean(B,1);
end

%calculate the shuffle average and plot
shuff_ave = mean(shuff_res,1);
shuff_std = std(shuff_res,0,1);

figure
errorbar(4:-1:1,shuff_ave,shuff_std)
%% Plot the clustered features
close all
%sort the indexes and the fish IDs
[~,s_ind] = sort(idx_clu);
figure
% subplot(1,8,1:7)
%     imagesc(delta_res(s_ind,:))
imagesc(delta_clu)
set(gca,'XTick',[1 2 3 4],'XTicklabels',{'Red','Green','Blue','UV'},...
    'FontSize',25)
ylabel('Trace #','FontSize',25)
xlabel('Cone','FontSize',25)
colormap(parula)
%average within each cluster
%allocate memory for the averages
code_ave = zeros(clu_num,size(delta_clu,2));
%and for the number of items in each cluster
clu_count = zeros(clu_num,1);
%for all the clusters
for clu = 1:clu_num
    %average
    code_ave(clu,:) = mean(delta_clu(idx_clu==clu,:),1);
    %and count
    clu_count(clu) = sum(idx_clu==clu);
end


% figure
% imagesc(clu_ave)
%plot them
%     [sort_m,sort_i] = sortrows(code_ave,4);
[sort_m,sort_i] = sort(clu_count,'descend');

figure
bar(sort_m)
set(gca,'FontSize',25)
xlabel('Cluster','FontSize',25)
ylabel('Number of seeds','FontSize',25)

figure
%     imagesc(sort_m)
imagesc(code_ave(sort_i,:))
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
%     set(gca,'XTick',[2.5,6.5,10.5],'XTickLabels',{'Gain at cone','Quality','Phase shift'}...
%         ,'YTick',1:size(code_ave,1),'YTickLabel',mat2cell(sort_i,ones(length(sort_i),1)))
set(gca,'XTick',[1 2 3 4],'XTicklabels',{'Red','Green','Blue','UV'},...
    'YTick',1:size(code_ave,1),'YTickLabel',mat2cell(sort_i,ones(length(sort_i),1)),...
    'FontSize',25);
ylabel('Clusters','FontSize',25)
colormap(parula)

%also plot the clustergram of the averages
clustergram(normr_1(code_ave,1),'Colormap',parula,'Cluster','column','Symmetric',false)
%% Use regressors to find cell types

%first z-score the matrix to only highlight the values above a
%certain threhold
zscored_delta = zscore(delta_norm,0,1);
%     zscored_delta = delta_norm;
%     %define the threshold (50% percentile)
%     z_thres = prctile(zscored_delta,50,1);
%     %now threshold
%     zscored_delta = bsxfun(@gt,zscored_delta,z_thres);

%     %define the 1 2,3 and 4 hit weights
%     hit1 = 1;
%     hit2 = 1;
%     hit3 = 1;
%     hit4 = 1;
%manually define the set of regressors
%     reg_set = [hit1 0 0 0;0 hit1 0 0;0 0 hit1 0;0 0 0 hit1;hit2 hit2 0 0;hit2 0 hit2 0;...
%         hit2 0 0 hit2;0 hit2 hit2 0;0 hit2 0 hit2;0 0 hit2 hit2;hit3 hit3 hit3 0;...
%         hit3 hit3 0 hit3;hit3 0 hit3 hit3;0 hit3 hit3 hit3;hit4 hit4 hit4 hit4];

%interval of regressors to use
reg_interval = [0 1];

reg_set = zeros(length(reg_interval).^4-2,4);



bit_c = 1;
for bit1 = reg_interval
    for bit2 = reg_interval
        for bit3 = reg_interval
            for bit4 = reg_interval
                reg_set(bit_c,:) = [bit1 bit2 bit3 bit4];
                bit_c = bit_c + 1;
            end
        end
    end
end

reg_set = reg_set(2:end-1,:);
%get the number of regressors
reg_num = size(reg_set,1);
%get the number of profiles
prof_num = size(zscored_delta,1);
%     %allocate memory to store the correlation values
%     reg_mat = zeros(reg_num,prof_num);
%     %for all of the regressors
%     for regs = 1:reg_num
%         %show the current regressor
%         fprintf(strcat('Current regressor:',num2str(regs),'\r\n'))
%         %for all the profiles
%         for prof = 1:prof_num
%             %calculate the correlation between the profiles and regressor
%             reg_mat(regs,prof) = sum(reg_set(regs,:).*zscored_delta(prof,1:stim_num2));
%         end
%     end

reg_mat = corr(zscored_delta',reg_set');

close all
figure
imagesc(sortrows(reg_mat))
%% Test
test_vec = [10 2 5 2];
test_val = zeros(reg_num,1);
for regs = 1:reg_num
    test_val(regs) = sum(reg_set(regs,:).*test_vec);
end
%% Select the top correlation values for each regressor and plot
close all
%     figure
%     imagesc(reg_mat')
%
%     %define the number to pick
%     reg_subnum = prof_num;
%     %allocate memory to store them
%     reg_picks = zeros(reg_subnum,reg_num);
%     %for all the regressors
%     for regs = 1:reg_num
%         %sort and pick the top 20 indexes
%         [~,temp_ind] = sort(reg_mat(regs,:),'descend');
%         reg_picks(:,regs) = temp_ind(1:reg_subnum);
%     end
%
%     %plot the results
%     figure
%     imagesc(reg_picks)

%     %plot the correlation between the sets of profiles
%     figure
%     rho = corr(reg_picks);
%     imagesc(rho)

%assign each trace to a particular type based on the best fitting
%regressor
%allocate memory for the indexes
best_reg = zeros(prof_num,1);
%for all of the profiles
for prof = 1:prof_num
    %get the index of the max reg value (i.e. the ID of the regressor)
    %         [~,best_reg(prof)] = max(reg_picks(prof,:));
    [~,best_reg(prof)] = max(reg_mat(prof,:));
end

%also get a vector with the trace counts for each one
reg_count = zeros(reg_num,1);
%for all of the regressors
for regs = 1:reg_num
    %count the number of traces for this regressor
    reg_count(regs) = sum(best_reg==regs);
end

c_map = [0 0 1;0 0 0;1 0 0];
figure
h1 = subplot(9,1,1:6);
imagesc(reg_set')
h2 = subplot(9,1,7:9)
imagesc(reg_count')
colormap(h1,c_map)
colormap(h2,parula)
%% Cluster the correlation patterns

%get rid of the NaNs
delta_clu = reg_mat;
delta_clu(isnan(delta_clu)) = 0;
%define the vector of cluster numbers to try
%     clu_vec = [5 10 20 30];
clu_vec = [50 100 150];
%     clu_vec = [5 6 7 8 9 10];
%     clu_vec = [100 300 500 1000 1100]
%     clu_vec = [80 90 100 150 200 400 600];
%     clu_vec = [32 33 34 35 36 37 38 40];
replicates = 1;
[idx_clu,GMModel,clu_num] = sPCA_GMM_cluster_4(normr_1(delta_clu,0),[],[],[],[],bic_name,clu_vec,replicates);
%% Cluster the correlations using k-means
%get rid of the NaNs
delta_clu = reg_mat;
delta_clu(isnan(delta_clu)) = 0;
clu_num = 16;
idx_clu = kmeans(normr_1(delta_clu,0),clu_num);
%% Find and store the indexes of the top target percentile of a set of regressors
%define the subset of regressors to take
reg_subset = [9 10 11 12];

%get their number
reg_subnum = length(reg_subset);

%allocate memory for the stored subsets
reg_subind = zeros(prof_num,reg_subnum);

%define the target percentile to look at
perc_tar = 85;

%for all of the elements in the subset
for subc = 1:reg_subnum
    
    %define the target regressor
    reg_tar = reg_subset(subc);
    
    %get the percentile
    reg_perc = prctile(reg_mat(reg_tar,:),perc_tar);
    
    %find the indexes of the profiles above the percentile
    reg_subind(:,subc) = reg_mat(reg_tar,:)>reg_perc;
    reg_subind = logical(reg_subind);
end
%% Plot the anatomy of a given regressor
close all

plot_save = 0;
%define the target percentile to look at
perc_tar = 90;

%define the target regressor
reg_tar = 1;

%get the percentile
reg_perc = prctile(reg_mat(reg_tar,:),perc_tar);

%find the indexes of the profiles above the percentile
tar_profiles = reg_mat(reg_tar,:)>reg_perc;

%also get the indexes of the fish and original seed
ori_seed = fish_ori(tar_profiles,:);

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
        
        %             %if the plotting variable is on
        %             if plot_save == 1
        %                 C = fused_im;
        %                 if fish == 1 && z == 1
        %                     [fig_name,fig_path] = uiputfile('E:\Behavioral data\Matlab\AF_proc\Figures\*.*');
        %                 end
        %                 fig_full = fullfile(fig_path,strcat(fig_name,'_',num2str(fish),'.tif'));
        %                 if z ==1
        %                     imwrite(C,fig_full,'tif','WriteMode','overwrite')
        %                 else
        %                     imwrite(C,fig_full,'tif','Resolution',size(C),'WriteMode','append')
        %                 end
        %             end
    end
    %     %calculate the max projection and plot
    %     figure
    %     imagesc(max(max_pro,[],3))
end
%% Use the gains and the linear model to re-obtain the traces

%allocate memory to store the traces
new_trace = zeros(size(all_trace));
%and to store the pre-calc equations
stim_set = zeros(stim_num2,cone_num,time_num);
%for all the stimuli
for stim = 1:stim_num2
    %load the color info for this stimulus
    col_info = squeeze(col_out(stim,:,:));
    %reshape to put the two sides in different dimensions
    col_info = reshape(col_info,size(col_info,1),4,2);
    %determine the LED that was on
    [~,led_on] = max(col_info(1,:,1));
    %load the values from that led for all the cones following the
    %stimulus
    stim_set(stim,:,:) = squeeze(cone_exc(:,led_on,col_info(:,led_on,1)+1));
    %make the rests 0
    stim_set(:,:,1:20) = 0;
    stim_set(:,:,61:80) = 0;
    %         %for all the leds
    %         for leds = 1:led_num
    %             %load the set of excitations
    %             curr_exc = squeeze(cone_exc(stim,leds,col_info(:,stim,1)+1));
    %             %calculate the response using the calculated model
    %         end
end

%for all the traces
for traces = 1:size(delta_res,1)
    %solve the equations for the expected responses at every time point
    %for all the time points
    for timep = 1:time_num
        new_trace(traces,timep,:) = delta_res(traces,:)*stim_set(:,:,timep)';
    end
    
end

%reshape the matrix
new_trace = reshape(new_trace,size(new_trace,1),size(new_trace,2)*size(new_trace,3));
%% Calculate the divergence between model and data
close all

figure
imagesc(new_trace)

%calculate squared differences
%     square_d = 1 - (sum((conc_trace-new_trace).^2,2)./sum((bsxfun(@minus,conc_trace,mean(conc_trace,2))).^2,2));
square_d = sum((conc_trace-new_trace).^2,2);
figure
histogram(square_d,1000)
%     set(gca,'XScale','log')

figure
plot(new_trace(1,:))
hold('on')
plot(conc_trace(1,:))
%% Save Gains analysis output

%define the save path
save_path = 'E:\Behavioral data\Matlab\AF_proc\Analysis\Stage3_cluster\';

save_var = 0;
if save_var == 1
    
    %get the root of the save name
    [ori_name,~] = uiputfile(strcat(save_path,'*.*'));
    %save the clustering output
    save_clu = strcat(ori_name,'_gains.mat');
    save(fullfile(save_path,save_clu),'conc_trace','stim_num2',...
        'time_num','col_out','delta_norm','fish_num')
end

load_var = 0;

if load_var == 1
    [f_name,f_path] = uigetfile(save_path);
    load(fullfile(f_path,f_name))
end
