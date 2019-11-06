% Model the interaction between regions

clearvars
close all
addpath(genpath('E:\Behavioral data\Matlab\AF_proc\ColorFishSuite'))

cluster_path = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\Stage3_cluster\';
fig_path = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Figures\Model\';

data = load_clusters(cluster_path);
%% Run the modelling 
% TODO: use the clusters per region
close all

% for all the data sets
for datas = 1:size(data,2)
    % get the region info
    region_info = data(datas).anatomy_info(:,1);
    % get the number of regions in the set
    unique_regions = unique(region_info);
    num_data = size(unique_regions(~isnan(unique_regions)),1);
    %get all the pairwise combinations of the fish
    region_comb = combnk(1:num_data,2);

    %get the number of combs
    num_comb = size(region_comb,1);
    % %get a vector to signal the rest periods
    % rest_vec = ones(20,1);
    % stim_vec = zeros(40,1);
    % 
    % rest_all = logical(repmat([rest_vec;stim_vec;rest_vec],stim_num2,1));

    % define the period of interest (0 pre, 1 stim, 2 post, 3 pre-post)
    period = 1;
    % get the target period labeled with ones
    rest_all = period_of_interest(period,data(datas).stim_num,1);
    %for all the combinations
    for combs = 1:num_comb
        %concatenate the clusters from both animals involved
%         tar1 = clu_all{fish_comb(combs,1),1};
%         tar2 = clu_all{fish_comb(combs,2),1}; 

        tar1 = data(datas).conc_trace(region_info==region_comb(combs,1),:);
        raw_trace2 = data(datas).conc_trace(region_info==region_comb(combs,2),:);
        tar2 = raw_trace2;
    %     %also load the raw traces
    %     raw_trace1 = load(name_cell{fish_comb(combs,1)},'conc_trace');
    %     raw_trace1 = raw_trace1.conc_trace;

%         raw_trace2 = load(name_cell{fish_comb(combs,2)},'conc_trace');
%         raw_trace2 = raw_trace2.conc_trace; 
% 
%         p_ind = randperm(size(raw_trace2,1));
%         train_trace = raw_trace2(p_ind(1:ceil(length(p_ind)/2)),:);
%         test_trace = raw_trace2(p_ind(ceil(length(p_ind)/2)+1:end),:);

        %for all the averages in 1, calculate models from the raw traces in 2
        %allocate memory to store the model results
        model_para = cell(size(tar1,1),1);
        %for all the averages
        for clu = 1:size(tar1,1)

            for stim = 1:data(datas).stim_num
    %         model_para{clu} = fitlm(raw_trace2',tar1(clu,:)','exclude',rest_all);
    %         model_para{clu} = fitglm(raw_trace2',tar1(clu,:)','exclude',rest_all);
    %         model_para{clu} = fitlm(tar2',tar1(clu,:)','linear','exclude',rest_all);
    %         figure
    %         plot(tar2(clu,:))
    %         hold('on')
    %         plot(model_para{clu}.Fitted)
    %        errorbar(1:length(tar2(clu,:)),model_para{clu}.Fitted,model_para{clu}.Residuals.Raw)

    %         model_para{clu} = fitglm(tar2',tar1(clu,:)','linear','exclude',rest_all,'Intercept',false);
    %         figure
    %         plot(tar2(clu,:))
    %         hold('on')
    %         plot(model_para{clu}.Fitted.Response)
    %         plot(model_para{clu}.Fitted.LinearPredictor)

                model_para{clu} = fitrlinear(raw_trace2',tar1(clu,:)','CrossVal','on');
                figure
                plot(tar2(clu,:))
                hold('on')
                plot(model_para{clu}.Y)
            end
    %         figure
    %         imagesc(model_para{clu}.CoefficientCovariance)
        end



    %         figure
    %         imagesc(model_para{1}.CoefficientCovariance)

    end
end