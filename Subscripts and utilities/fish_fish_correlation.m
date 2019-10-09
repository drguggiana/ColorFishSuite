% fish to fish correlation plots

%% Given the clusters, add the fish identity information
close all
%allocate memory for the cluster averages per fish
fish_clu = zeros(size(conc_trace,2),fish_num,clu_num);
%also allocate memory to store how many traces go into each average
fish_count = zeros(fish_num,clu_num);
%for all the fish
for fish = 1:fish_num
    %for all the clusters
    for clu = 1:clu_num
        %get the trace indexes from the corresponding fish and cluster
        tar_traces = fish_ori(:,1)==fish&idx_clu==clu;
        %if there are none
        if isempty(tar_traces)
            continue
        end
        %otherwise fill in the data
        %fill in the matrix
        fish_clu(:,fish,clu) = mean(conc_trace(tar_traces,:),1);
        
        %and the counting one too
        fish_count(fish,clu) = sum(tar_traces);
    end
end

%get rid of the NaNs from the empty clusters
fish_clu(isnan(fish_clu)) = 0;

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

%plot the counts for each fish and cluster
figure
imagesc(fish_count(:,ind_count))
xlabel('Cluster','Fontsize',25)
ylabel('Fish','Fontsize',25)
colormap(parula)
set(gca,'FontSize',25)
colorbar
%plot the profiles for each fish

figure

%for all the fish
for fish = 1:fish_num
    subplot(round(sqrt(fish_num)),ceil(sqrt(fish_num)),fish);
    imagesc(squeeze(fish_clu(:,fish,ind_count))')
    set(gca,'XTick',[],'YTick',[])
%     colormap(parula)
end

%get all the pairwise combinations ofthe fish
fish_comb = combnk(1:fish_num,2);

%get the number of combs
num_comb = size(fish_comb,1);

%allocate memory for the correlation coeffs
fish_corr = zeros(fish_num,fish_num,clu_num);
%for all the combinations
for combs = 1:num_comb
    %select the 2 targets
    tar1 = squeeze(fish_clu(:,fish_comb(combs,1),:));
    tar2 = squeeze(fish_clu(:,fish_comb(combs,2),:));
    %for each cluster
    for clu = 1:clu_num
        %calculate the pairwise correlations between each cluster average
        temp_corr = corrcoef(tar1(:,clu),tar2(:,clu));
        fish_corr(fish_comb(combs,1),fish_comb(combs,2),clu) = temp_corr(1,2);
    end
end

figure
imagesc(nanmean(fish_corr,3))
% colormap(parula)
colorbar
set(gca,'FontSize',25)
xlabel('Fish','FontSize',25)
ylabel('Fish','FontSize',25)

%allocate memory for the lumped info
stim_corr = zeros(clu_num,2);
%for all the clusters
for clu = 1:clu_num
    %temp variable
    temp_var = fish_corr(:,:,clu);
    %calculate the average across fish comparisons
    stim_corr(clu,1) = nanmean(temp_var(:));
    %and the std
    stim_corr(clu,2) = nanstd(temp_var(:));
end
figure
errorbar(1:clu_num,stim_corr(:,1),stim_corr(:,2)./sqrt(fish_num))