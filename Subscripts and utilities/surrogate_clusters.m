% surrogate clusters
%% Calculate surrogate clusterings and assess cluster quality

%define the proportion of the data set to use
prop_data = 0.9;

%define the number of surrogates
num_surr = 20;

%fourier cluster?
four_clu = 0;
%allocate memory to store the indexes of the surrogates selected, the GMMs
%,the idx_clu of the surrogate data sets and the cluster averages
surr_cell = cell(num_surr,4);

%for all the surrs
for surr = 1:num_surr
    %if it's clustering the fourier components (i.e. no spca)
    if four_clu == 1
        clu_mat = delta_clu;
    else 
        clu_mat = conc_trace;
    end
    
    %generate a random subset of indexes of the main data set
    rand_ind = randperm(size(clu_mat,1));
    rand_ind = rand_ind(1:round(0.9*size(clu_mat,1)));
    %store the random subset of indexes
    surr_cell{surr,1} = rand_ind;
    %cluster the subset of the data set
%     %define the vector of cluster numbers to try
%     clu_vec = [5 10 20 30 40 50 70 100];
%     clu_vec = [80 90 100 150 200];
    replicates = 1;
    %if it's clustering the fourier components (i.e. no spca)
    if four_clu == 1
        [idx_surr,GMModel_surr,clu_num] = sPCA_GMM_cluster_3(clu_mat,ori_clu...
            ,[],[],[],[],bic_name,clu_vec,replicates);
    else
       
        %define the sPCA parameters to use
        bounds_top = 1:time_num:size(clu_mat,2);
        bounds_bottom = [bounds_top(2:end)-1,size(clu_mat,2)];
        bounds = [bounds_top;bounds_bottom];
        K = ones(1,stim_num2).*8;
        t_bins = ones(stim_num2,1).*10;
        pca_vec = ones(stim_num2,1).*1;
        [idx_surr,GMModel_surr,clu_num,pcs] = sPCA_GMM_cluster_3(clu_mat,clu_num,bounds...
            ,K,t_bins,pca_vec,bic_name,clu_vec,replicates);
    end
    %regardless, save the clustering output
    surr_cell{surr,2} = GMModel_surr;
    surr_cell{surr,3} = idx_surr;
    %calculate the cluster averages
    %allocate memory for the cluster average
    clu_ave_surr = zeros(clu_num,size(clu_mat,2));
    %for all the clusters
    for clu = 1:clu_num
        clu_ave_surr(clu,:) = mean(clu_mat(idx_surr==clu,:),1);
    end
    %store the average in the cell
    surr_cell{surr,4} = clu_ave_surr;
    
end
%% Save the surrogate cell (or load)
%define the save path
save_path = 'E:\Behavioral data\Matlab\AF_proc\Analysis\Stage3_cluster\';

save_var = 1;
if save_var == 1

    %get the root of the save name
    [ori_name,~] = uiputfile(strcat(save_path,'*.*'));
    %save the clustering output
    save_clu = strcat(ori_name,'_CONTROLclusters.mat');
    save(fullfile(save_path,save_clu),'surr_cell')
end
%% use the surrogate cell to assess the quality of the clusters

%allocate memory for the resulting correlation per cluster
clu_qual = zeros(num_surr,clu_num);

%if it's the fourier matrix
if four_var == 1
    real_clu = code_ave;
else
    real_clu = clu_ave;
end

%for all of the surrogates
for surr = 1:num_surr
    %load the current surrogate average cluster matrix
    temp_surr = surr_cell{surr,4};
    %calculate the correlation matrix between the current surrogate and the
    %original matrix
    rho = corr(real_clu',temp_surr');
    %generate a list of clusters of the first input to pair to the second
    clu_list = zeros(clu_num,1);
    
    %for all the clusters
    for clu = 1:clu_num
        %find the maximum in the correlation matrix
        [max_val,max_ind] = max(rho(:));
        %turn the linear index into x,y
        [max_row,~] = ind2sub([clu_num clu_num],max_ind);
        
        %store it in the output vector
        clu_list(max_row) = max_val;
        
        %NaN the entire row so that this cluster doesn't show up anymore
        rho(max_row,:) = NaN;
    end
    %store the final vector into the output matrix
    clu_qual(surr,:) = clu_list;
end

%calculate the median for each cluster and plot
med_surr = median(clu_qual,1);

close all

figure
imagesc(clu_qual)
ylabel('Replicate','FontSize',20)
xlabel('Cluster','FontSize',20)
figure
plot(med_surr)
ylabel('Correlation median across replicates','FontSize',20)
xlabel('Cluster','FontSize',20)