function [idx_clu,clu_num] = cluster_snr(snr_mat,clu_num,idx_clu,num_thres,stim_thres)

% %load the snr for the dataset
% snr_mat = load(name_cell{files},'snr_mat');
% snr_mat = snr_mat.snr_mat;

%calculate the average snr vector for each cluster

%allocate memory for the clusters
clu_snr = zeros(clu_num,size(snr_mat,2));

%for all the clusters
for clu = 1:clu_num
    
    %load the average matrix
    clu_snr(clu,:) = mean(snr_mat(idx_clu==clu,:),1);
end

% %plot a histogram of the averages
% figure
% for stim = 1:stim_num2
%     subplot(round(sqrt(stim_num2)),ceil(sqrt(stim_num2)),stim)
%     histogram(clu_snr(:,stim),100)
% end

%find a threshold by defining the 50th percentile on each stimulus
snr_thres = prctile(clu_snr,50,2);

%keep traces that pass the threshold in at least 1 stim
snr_thres_all = bsxfun(@gt,clu_snr,snr_thres);
snr_thres_vec = sum(snr_thres_all,2)>stim_thres;
%also modify the vector based on the number of traces per cluster
% num_thres = 10;
%for all the clusters
for clu = 1:clu_num
    %kill the ones with less than num_thres traces
    if sum(idx_clu==clu)<num_thres
        snr_thres(idx_clu==clu) = 0;
        idx_clu(idx_clu==clu) = 0;
        snr_thres_vec(clu) = 0;
    end
end

%modify the idx_clu vector and clu_num to reflect this fact
%find all the clusters that stay
new_clu = find(snr_thres_vec);
%create a new idx_clu to overwrite the old one
new_idx = zeros(size(idx_clu));
%initialize a new cluster counter
c_count = 1;
%for all the clusters
for clu = 1:clu_num
    %check whether the cluster number stays
    if sum(new_clu==clu)>0
        %renumber all instances in the new idx vector
        new_idx(idx_clu==clu) = c_count;
        %update the counter
        c_count = c_count + 1;
        %clusters not found will be left as zeros
    end
end

%overwrite the old variables
idx_clu = new_idx;
clu_num = numel(unique(idx_clu(idx_clu>0)));