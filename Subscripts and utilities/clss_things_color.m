function [class_cell] = clss_things_color(conc_trace,loo,trace_frac,set_part,redec_num,shuff_label,classpcolor,bin_width,stim_num2,rep_num,period)

% close all

% %define whether to use leave one out cross validation (0 is no)
% loo = 1;
% %define the fraction of the traces to keep
% trace_frac = 1;
% %define the partition to use
% set_part = 1;
% %set the number of reps
% redec_num = 1;
% %define whether to shuffle labels (for neutral classification)
% shuff_label = 0;
% %define wether to use clusters or the whole data
% all_data = 0;
% %define the number of classes per color (1,3,5,or 8) (or 10,11 and 12 for the
% %p6p8 data)
% classpcolor = 5;
% %define the binning factor
% bin_width = 1;
%define the vector of lambdas to try
% lambda_vec = [0.1 0.3 1 3 10 30 100];
% lambda_vec = [300 1000 3000 10000 30000 100000];
% lambda_vec = 5;
lambda_vec = 'auto';
%select the learning method
template_log = templateLinear('Learner','svm','Lambda',lambda_vec);
learner_vec = {'svm','discriminant','knn','linear','naivebayes','tree',template_log};
learn_var = 4;

%allocate memory to store the results from the analysis
%1 conc conf matrix,2 conc diagonal frac,3 shuff ave and shuff std
%4 mut_inf and shuff_mutinf
class_cell = cell(3,1);


%remove the empty entries in the averages
conc_trace = conc_trace(sum(conc_trace,2)>0,:);

%normalize the entries (feature scaling for the svm)
% conc_trace = normr_1(conc_trace,[]);

%allocate memory for the results
redec_cell = cell(redec_num,2);

%for all the reps
for redec = 1:redec_num
    
    %define the target matrix
    tar_mat = conc_trace;
    %extract a random subset of the traces out
    rand_trace = randperm(size(tar_mat,1));
    
    %keep only that fraction of traces
    tar_mat = tar_mat(rand_trace(1:round(size(tar_mat,1)*trace_frac)),:);
    
    % tar_mat = raw_trace + rand(size(raw_trace)).*0.01;
    % tar_mat = clu_ave;
    %get the number of traces
    trace_num = size(tar_mat,1);
    % get the period of interest labeled with ones
    rest_all = period_of_interest(period,stim_num2,rep_num);
    %             rest_all = logical(ones(80*stim_num2,1));
    obs_num = sum(rest_all);
%     %calculate the time per stimulus
%     t_perstim = obs_num/stim_num2;
    %allocate memory to store the restless data
    stim_only = zeros(trace_num,obs_num);
    %extract the non-rest periods
    %for all the traces
    for tra = 1:trace_num
        %extract the stim only info
        stim_only(tra,:) = tar_mat(tra,rest_all);
    end
    
    %time-bin the input matrix
    
    %get the number of bins
    bin_num = ceil(obs_num/bin_width);
    % %get the indexes for each bin
    % [~,~,bin_ind] = histcounts(1:obs_num,'BinWidth',bin_width);
    %bin the array
    %allocate memory for the binned array
    stim_bin = zeros(trace_num,bin_num);
    
    %initialize a counter for the time
    time_c = 1;
    %for all the time points
    for bin_var = 1:bin_num
        %actually perform the binning
        stim_bin(:,bin_var) = mean(stim_only(:,time_c:time_c+bin_width-1),2);
        %update the time counter
        time_c = time_c + bin_width;
    end
    %replace the original matrix
    stim_only = stim_bin;
    %also update the observation counter time per stimulus
    obs_num = bin_num;
    t_perstim = obs_num/stim_num2;
    
    %if 5 color categories are desired
    if (classpcolor == 5)||(classpcolor == 6)
        %bin the matrix matching the redundant intensity points
        
        %allocate memory to store the averaged data
        ave_stim = zeros(trace_num,5,obs_num/8);
        %reshape the data to be able to average
        stim_resh = reshape(stim_only,trace_num,8,obs_num/8);
        %for each one of the averaging intervals
        for interv = 1:5
            switch interv
                case {1,5}
                    %select the corresponding points and average
                    ave_stim(:,interv,:) = stim_resh(:,interv+2,:);
                case 2
                    %select the corresponding points and average
                    ave_stim(:,interv,:) = stim_resh(:,4,:);
                case 3
                    %select the corresponding points and average
                    ave_stim(:,interv,:) = stim_resh(:,5,:);
                case 4
                    %select the corresponding points and average
                    ave_stim(:,interv,:) = stim_resh(:,6,:);
            end
        end
        %reshape the matrix back to 2D shape
        stim_only = reshape(ave_stim,trace_num,5*obs_num/8);
        %recalculate the obs_num
        obs_num = size(stim_only,2);
        %also redefine the time per stimulus
        t_perstim = obs_num/stim_num2;
    end
    
    %allocate memory for the stimulus labels
    stim_label = zeros(obs_num,1);
    %create the weights vector
    weight_vec = ones(size(stim_label));
    %select the appropriate labeling scheme
    switch classpcolor
        case 1
            %create a counter for the labels
            label_c = 1;
            %get the width of each stimulus
            stim_width = obs_num/stim_num2;
            %for all the stimuli
            for stim = 1:stim_num2
                %create a vector with the stimulus labels
                stim_label(label_c:label_c + stim_width-1) = stim-1;
                %update the counter
                label_c = label_c + stim_width;
            end
        case 5
            %define a custom stim label (based on the LED intensity levels,
            %look at col_out)
            %                 %5 classes per color
            %                 stim_label(1:4:end) = 3;%mid level
            %                 stim_label(7:8:end) = 5;%top
            %                 stim_label(3:8:end) = 1;%bottom
            %                 stim_label(2:8:end) = 2;%mid bottom 1
            %                 stim_label(4:8:end) = 2;%mid_bottom 2
            %                 stim_label(6:8:end) = 4;%mid top 1
            %                 stim_label(8:8:end) = 4;%mid top 2
            %                 %fill up the weights vector
            %                 weight_vec(7:8:end) = 0.8;
            %                 weight_vec(3:8:end) = 0.8;
            %5 color classes, disregarding the repeat levels (p17b)
            stim_label = mod(0:obs_num-1,5)+1;
        case 3
            %3 classes per color (p17b)
            stim_label(1:4:end) = 2;%mid level
            stim_label(7:8:end) = 3;%top
            stim_label(3:8:end) = 1;%bottom
            stim_label(2:8:end) = 1;%mid bottom 1
            stim_label(4:8:end) = 1;%mid_bottom 2
            stim_label(6:8:end) = 3;%mid top 1
            stim_label(8:8:end) = 3;%mid top 2
            %fill up the weights vector
            weight_vec(1:4:end) = 1.5;
        case 8
            %8 classes per color (p17b)
            stim_label(1:8:end) = 1;%mid level
            stim_label(5:8:end) = 5;%mid level2
            stim_label(7:8:end) = 7;%top
            stim_label(3:8:end) = 3;%bottom
            stim_label(2:8:end) = 2;%mid bottom 1
            stim_label(4:8:end) = 4;%mid_bottom 2
            stim_label(6:8:end) = 6;%mid top 1
            stim_label(8:8:end) = 8;%mid top 2
        case 10
            %this is actually stim type labelling for the p6p8 data
            %for both colors
            for c_type = 1:2
                stim_label(1+(t_perstim*(c_type*3-3)):(t_perstim*(c_type*3-2))) = 1;
                stim_label(1+(t_perstim*(c_type*3-2)):(t_perstim*(c_type*3-1))) = 2;
                stim_label(1+(t_perstim*(c_type*3-1)):(t_perstim*(c_type*3))) = 3;
            end
        case 11
            %this is color type labelling for the p6p8 data
            %for all 3 types of stimuli
            for s_type = 1:3
                stim_label(1+(t_perstim*((s_type*2-1)-1)):(t_perstim*(s_type*2-1))) = 1;
                stim_label(1+(t_perstim*(s_type*2-1)):(t_perstim*(s_type*2))) = 2;
            end
        case 12
            %this is both stim and color labelling for the p6p8 data
            %for all 3 types of stimuli
            for s_type = 1:3
                stim_label(1+(t_perstim*((s_type*2-1)-1)):(t_perstim*(s_type*2-1))) = 1+(s_type*2-2);
                stim_label(1+(t_perstim*(s_type*2-1)):(t_perstim*(s_type*2))) = 1+(s_type*2-1);
            end
        case  6
            %label only intensity (iffy cause of the actual luminance
            %levels)
            %for all 5 levels of intensity
            for i_type = 1:5
                stim_label(1+(i_type-1):5:end) = i_type;
            end
%             assignin('base','stim_label',stim_label)
            
    end
    %if not using 1 class per color
    if classpcolor > 1 && classpcolor < 10 && classpcolor ~= 6
        %get the highest number on the first category
        max_cat = max(stim_label);
        %for each stimulus increase the number of stimuli
        %for all the reps
%         for reps = 1:rep_num
            for stim = 2:4
                stim_label((stim-1)*t_perstim+1:stim*t_perstim) = ...
                    stim_label((stim-1)*t_perstim+1:stim*t_perstim)+(stim-1)*max_cat;
            end
%         end
    end
    %if shuffling of labels is on
    if shuff_label == 1
        stim_label = stim_label(randperm(length(stim_label)));
    end
%     assignin('base','stim_label',stim_label)
    %get the number of categories
    categ_num = length(unique(stim_label));
    
    %initialize the train and test label variables
    train_label = zeros(round(length(stim_label)*set_part),1);
    test_label = zeros(length(stim_label)-length(train_label),1);
    
    %shuffle only if not loo
    if ~loo
        %keep randomizing until the sets used contain members from all stimuli
        while length(unique(train_label))<categ_num || length(unique(test_label))<categ_num
            %define the training and test sets at random
            %get the vector of randomized indices
            % rng(0)
            rand_ind = randperm(obs_num);
            % rng('default')
            
            %get the sets
            train_set = stim_only(:,rand_ind(1:ceil(obs_num*set_part)))';
            train_label = stim_label(rand_ind(1:ceil(obs_num*set_part)));
            
            %if there is a partition
            if set_part ~= 1
                %use the rest of the data set for calculations
                test_set = stim_only(:,rand_ind(1+ceil(obs_num*set_part):end))';
                test_label = stim_label(rand_ind(1+ceil(obs_num*set_part):end));
            else
                %if not, use again the whole data set for testing
                test_set = train_set;
                test_label = train_label;
            end
        end
    else
        %get the sets
        train_set = stim_only';
        train_label = stim_label;
        
        test_set = train_set;
        test_label = train_label;
    end
    
    %if leave one out
    if loo == 1
        %set parallel computing parameters
        options = statset('UseParallel',0);
        %calculate the decoder
        af_deco = fitcecoc(train_set,train_label,'LeaveOut','on','Verbose',1,...
            'Learners',learner_vec{learn_var},'Coding','onevsall','Prior','Empirical',...
            'Weights',weight_vec,'Options',options);
        
        %use the decoder to predict the data
        %             af_pred = predict(af_deco,test_set);
        af_pred = kfoldPredict(af_deco);
    else
        %calculate the decoder
        af_deco = fitcecoc(train_set,train_label,'LeaveOut','off','Verbose',1,...
            'Learners',learner_vec{learn_var},'Coding','onevsall');
        
        %use the decoder to predict the data
        af_pred = predict(af_deco,test_set);
    end
    
    %calculate a confusion matrix
    redec_cell{redec,1} = confusionmat(test_label,af_pred);
    %and the fraction of values in the diagonal
    redec_cell{redec,2} = sum(diag(redec_cell{redec,1}))/sum(redec_cell{redec,1}(:));
    
end

%concatenate the conf matrices
redec_conf = cat(3,redec_cell{:,1});
%and the diagonal fractions
redec_frac = horzcat(redec_cell{:,2});
%calculate the average conf mat
conf_mat = mean(redec_conf,3);
%     %plot and report the averages (normalizing to sum per column)
%     figure
%     imagesc(normr_1(conf_mat,4))
% %     title('Confusion Matrix Mean','FontSize',20)
%     xlabel('Predicted Stimulus','FontSize',35)
%     ylabel('Observed Stimulus','FontSize',35)
% %     set(gca,'XTick',1:2:categ_num,'FontSize',20)
%     set(gca,'XTick',[],'YTick',[])
%     axis('image')

%     figure
%     imagesc(std(redec_conf,0,3))
%     title('Confusion Matrix STD','FontSize',20)
%     xlabel('Predicted Stimulus','FontSize',20)
%     ylabel('Observed Stimulus','FontSize',20)
%     set(gca,'XTick',1:2:categ_num,'FontSize',20)

%     ave_frac = mean(redec_frac)

%calculate average diag fraction of shuffled data
%define the number of reps
mut_rep = 1000;
%initialize the variable
diag_shuff = zeros(mut_rep,1);
%for all the reps
for reps = 1:mut_rep
    rand_ind = randperm(numel(conf_mat));
    shuff_mat = reshape(conf_mat(rand_ind),size(conf_mat));
    diag_shuff(reps) = sum(diag(shuff_mat))/sum(conf_mat(:));
end

shuff_mean = mean(diag_shuff);
shuff_std = std(diag_shuff);

%store the classifier output
class_cell{1} = redec_conf;
class_cell{2} = redec_frac;
class_cell{3} = [shuff_mean,shuff_std];
    
%     %calculate the mutual information between the observed and predicted
%     mut_inf = mutualinfo(test_label,af_pred);
%     %also, shuffle the prediction and calculate the mutual information. Repeat
%     %100 times
%     %define the number of reps
%     mut_rep = 100;
%     %allocate memory for the calculation
%     mut_vec = zeros(mut_rep,1);
%     %for all the times
%     for reps = 1:mut_rep
%         mut_vec(reps) = mutualinfo(test_label,af_pred(randperm(length(af_pred))));
%     end
%     
%     shuff_mutinf = mean(mut_vec);
%     %store the mutual info data
%         class_cell{datas,4} = [mut_inf,shuff_mutinf];



% %report the mean average fraction
% ave_diag = mean(vertcat(class_cell{:,2}),2);