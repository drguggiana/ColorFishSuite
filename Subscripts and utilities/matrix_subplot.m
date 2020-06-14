function [rho_cell,combo_list,combo_number] = matrix_subplot(data,index_cell)


% % define the number of shuffles
% num_shuffles = 500;
% 
% % allocate memory for the legend handles
% legend_cell = cell(num_data,1);

% get the number of stimuli
stim_num = data(1).stim_num;
% get the combinations of stimuli and their number
combo_list = nchoosek(1:stim_num,2);
combo_number = size(combo_list,1);

% allocate memory for the correlations
rho_cell = cell(stim_num,stim_num,2);

% load the clusters and get the cluster averages
for datas = 1:2
    % load the data and clusters
    conc_trace = data(datas).conc_trace(index_cell{datas}==1,:);
    conc_trace = reshape(conc_trace,[],data(datas).time_num,stim_num);
    idx = data(datas).idx_clu(index_cell{datas}==1);
    clu_num = data(datas).clu_num;
    % for all the combinations
    for combo = 1:combo_number
        % load the particular color combination
        target_trace = reshape(conc_trace(:,21:60,combo_list(combo,:)),size(conc_trace,1),[]);

        % allocate memory for the averages
        temp_ave = zeros(clu_num,size(target_trace,2));
        for clu = 1:clu_num
            temp_ave(clu,:) = mean(target_trace(idx==clu,:),1);
        end
        % store the correlations
        [rho_cell{combo_list(combo,1),combo_list(combo,2),datas},pval1] = ...
            corr(temp_ave(:,1:40)',temp_ave(:,41:end)');
        
    end

end