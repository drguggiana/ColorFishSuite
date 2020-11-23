function [class_cell] = clss_things_color(conc_trace,loo,trace_frac,set_part,redec_num,shuff_label,classpcolor,bin_width,stim_num2,rep_num,period)

% close all

lambda_vec = 'auto';
%select the learning method
template_log = templateLinear('Learner','svm','Lambda',lambda_vec);
learner_vec = {'svm','discriminant','knn','linear','naivebayes','tree',template_log};
learn_var = 4;

%allocate memory to store the results from the analysis
%1 conc conf matrix,2 conc diagonal frac,3 shuff ave and shuff std
%4 predicted vectors, 5 label vector
class_cell = cell(6,1);

%allocate memory for the results
redec_cell = cell(redec_num,5);

%for all the reps
for redec = 1:redec_num
    
    %define the target matrix
    tar_mat = conc_trace;
    %extract a random subset of the traces out
    rand_trace = randperm(size(tar_mat,1));
    
    %keep only that fraction of traces
    tar_mat = tar_mat(rand_trace(1:round(size(tar_mat,1)*trace_frac)),:);

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
    t_perstim = obs_num/(stim_num2*rep_num);
    
    %if 5 color categories are desired
    switch classpcolor
        case {5,6}
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
            stim_width = obs_num/(stim_num2*rep_num);
            % for all the reps
            for reps = 1:rep_num
                %for all the stimuli
                for stim = 1:stim_num2
                    %create a vector with the stimulus labels
                    stim_label(label_c:label_c + stim_width-1) = stim-1;
                    %update the counter
                    label_c = label_c + stim_width;
                end
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
        case {10,22}
            %this is actually stim type labelling for the p6p8 data
            % for all the reps
            for reps = 1:rep_num
%                 %for both colors
%                 for c_type = 1:2
                stim_label(1+t_perstim*6*(reps-1):...
                    t_perstim*2+t_perstim*6*(reps-1)) = 1;
                stim_label(1+t_perstim*2+t_perstim*6*(reps-1):...
                    t_perstim*4+t_perstim*6*(reps-1)) = 2;
                stim_label(1+t_perstim*4+t_perstim*6*(reps-1):...
                    t_perstim*6+t_perstim*6*(reps-1)) = 3;
%                 end
            end
        case {11,23}
            %this is color type labelling for the p6p8 data
            % for all the reps
            for reps = 1:rep_num
                %for all 3 types of stimuli
                for s_type = 1:3
                    stim_label(1+(t_perstim*((s_type*2-1)-1))+t_perstim*6*(reps-1):...
                        (t_perstim*(s_type*2-1))+t_perstim*6*(reps-1)) = 1;
                    stim_label(1+(t_perstim*(s_type*2-1))+t_perstim*6*(reps-1):...
                        (t_perstim*(s_type*2))+t_perstim*6*(reps-1)) = 2;
                end
            end
        case {12,24}
            %this is both stim and color labelling for the p6p8 data

            % for all the reps
            for reps = 1:rep_num
                %for all 3 types of stimuli
                for s_type = 1:3
                    stim_label(1+(t_perstim*((s_type*2-1)-1))+t_perstim*6*(reps-1):...
                        (t_perstim*(s_type*2-1))+t_perstim*6*(reps-1)) = 1+(s_type*2-2);
                    stim_label(1+(t_perstim*(s_type*2-1))+t_perstim*6*(reps-1):...
                        (t_perstim*(s_type*2))+t_perstim*6*(reps-1)) = 1+(s_type*2-1);
                end
            end
        case  6
            %label only intensity (iffy cause of the actual luminance
            %levels)
            %for all 5 levels of intensity
            for i_type = 1:5
                stim_label(1+(i_type-1):5:end) = i_type;
            end
        case 13
            % label only UV and red, exclude the other colors
            % get rid of the other colors
            stim_only = stim_only(:,[1:t_perstim,3*t_perstim+1:4*t_perstim]);
            % set the labels
            %             stim_label = stim_label(1:2*t_perstim);
            stim_label = cat(2,ones(1,t_perstim),2.*ones(1,t_perstim));
            % also change the weight
            weight_vec = ones(2*t_perstim,1);
        case 14
            % label only CK colors
            % get rid of the other colors
            stim_only = stim_only(:,1:2*t_perstim);
            % set the labels
            %             stim_label = stim_label(1:2*t_perstim);
            stim_label = cat(2,ones(1,t_perstim),2.*ones(1,t_perstim));
            % also change the weight
            weight_vec = ones(2*t_perstim,1);
        case 15
            % label only GR colors
            % get rid of the other colors
            stim_only = stim_only(:,2*t_perstim+1:4*t_perstim);
            % set the labels
            %             stim_label = stim_label(1:2*t_perstim);
            stim_label = cat(2,ones(1,t_perstim),2.*ones(1,t_perstim));
            % also change the weight
            weight_vec = ones(2*t_perstim,1);
        case 16
            % label only flash colors
            % get rid of the other colors
            stim_only = stim_only(:,4*t_perstim+1:6*t_perstim);
            % set the labels
            %             stim_label = stim_label(1:2*t_perstim);
            stim_label = cat(2,ones(1,t_perstim),2.*ones(1,t_perstim));
            % also change the weight
            weight_vec = ones(2*t_perstim,1);
        case 17
            if bin_width > 1
                error('bin width must be 1')
            end
            % leave only max UV and 127 RGB to have more matching colors
            stim_only = reshape(stim_only,trace_num,[],stim_num2,3);
            stim_only = cat(3,stim_only(:,[1,9,17,25,33],1:3,:),stim_only(:,[7,15,23,31,39],4,:));
            stim_only = reshape(stim_only,trace_num,[]);
            stim_label = repmat(cat(2,ones(1,5),2.*ones(1,5),3.*ones(1,5),4.*ones(1,5)),1,3);
            weight_vec = ones(size(stim_label,2),1); 
        case 18
            % label only pre stim and stim CK
            % get the length of a single cycle
            cycle_length = size(stim_only,2)/(stim_num2*rep_num);
            
            % get rid of the unused data
            stim_only = reshape(stim_only,trace_num,cycle_length,stim_num2,rep_num);
            stim_only = reshape(stim_only(:,:,[1 2],:),trace_num,[]);
            stim_label = repmat(cat(2,ones(1,cycle_length/2),2.*ones(1,cycle_length/2)),1,rep_num*2);
            weight_vec = ones(size(stim_label,2),1);
        case 19
            % label only pre stim and stim GR
            % get the length of a single cycle
            cycle_length = size(stim_only,2)/(stim_num2*rep_num);
            
            % get rid of the unused data
            stim_only = reshape(stim_only,trace_num,cycle_length,stim_num2,rep_num);
            stim_only = reshape(stim_only(:,:,[3 4],:),trace_num,[]);
            stim_label = repmat(cat(2,ones(1,cycle_length/2),2.*ones(1,cycle_length/2)),1,rep_num*2);
            weight_vec = ones(size(stim_label,2),1);
        case 20
            % label only pre stim and stim FL
            % get the length of a single cycle
            cycle_length = size(stim_only,2)/(stim_num2*rep_num);
            
            % get rid of the unused data
            stim_only = reshape(stim_only,trace_num,cycle_length,stim_num2,rep_num);
            stim_only = reshape(stim_only(:,:,[5 6],:),trace_num,[]);
            stim_label = repmat(cat(2,ones(1,cycle_length/2),2.*ones(1,cycle_length/2)),1,rep_num*2);
            weight_vec = ones(size(stim_label,2),1);
        case 21
            % label only red and UV in p17b dataset
            % get the length of a single cycle
            cycle_length = size(stim_only,2)/(stim_num2*rep_num);
            
            % get rid of the unused data
            stim_only = reshape(stim_only,trace_num,cycle_length,stim_num2,rep_num);
            stim_only = reshape(stim_only(:,:,[1 4],:),trace_num,[]);
            stim_label = repmat(cat(2,ones(1,cycle_length/2),2.*ones(1,cycle_length/2)),1,rep_num*2);
            weight_vec = ones(size(stim_label,2),1);
    end
    
    % if it's the stimuli without flash, exclude the flash
    if any([22,23,24]==classpcolor)
        % reshape
        stim_only = reshape(stim_only,[],t_perstim,stim_num2,rep_num);
        stim_label = reshape(stim_label,1,t_perstim,stim_num2,rep_num);
        % redefine the number of stimuli
        stim_num2 = 4;
        % remove the flash stimuli
        stim_only = reshape(stim_only(:,:,1:stim_num2,:),[],t_perstim*stim_num2*rep_num);
        stim_label = reshape(stim_label(1,:,1:stim_num2,:),1,t_perstim*stim_num2*rep_num);
    end
    %if not using 1 class per color
    if classpcolor > 1 && classpcolor < 10 && classpcolor ~= 6
        %get the highest number on the first category
        max_cat = max(stim_label);
        %for each stimulus increase the number of stimuli
        for stim = 2:4
            stim_label((stim-1)*t_perstim+1:stim*t_perstim) = ...
                stim_label((stim-1)*t_perstim+1:stim*t_perstim)+(stim-1)*max_cat;
        end
    end
    %if shuffling of labels is on
    if shuff_label == 1
        stim_label = stim_label(randperm(length(stim_label)));
    end
%     %get the number of categories
%     categ_num = length(unique(stim_label));
    
    %initialize the train and test label variables
%     train_label = zeros(round(length(stim_label)*set_part),1);
%     test_label = zeros(length(stim_label)-length(train_label),1);

    %get the sets
    train_set = stim_only';
    train_label = stim_label;

%     test_set = train_set;
    test_label = train_label;
    
    %set parallel computing parameters
    options = statset('UseParallel',0);
    %if leave one out
    if loo == 1

        %calculate the decoder
        af_deco = fitcecoc(train_set,train_label,'LeaveOut','on','Verbose',0,...
            'Learners',learner_vec{learn_var},'Coding','onevsall','Prior','Empirical',...
            'Weights',weight_vec,'Options',options);
        
        %use the decoder to predict the data
        %             af_pred = predict(af_deco,test_set);
        af_pred = kfoldPredict(af_deco);
    else
        %calculate the decoder
        af_deco = fitcecoc(train_set,train_label,'KFold',5,'Verbose',0,...
            'Learners',learner_vec{learn_var},'Coding','onevsall','Prior','Empirical',...
            'Options',options);
        
        %use the decoder to predict the data
        af_pred = kfoldPredict(af_deco);
    end
    
    %calculate a confusion matrix
    redec_cell{redec,1} = confusionmat(test_label,af_pred);
    %and the fraction of values in the diagonal
    redec_cell{redec,2} = sum(diag(redec_cell{redec,1}))/sum(redec_cell{redec,1}(:));
    % also save the predicted and label vectors
    redec_cell{redec,3} = af_pred;
    redec_cell{redec,4} = test_label;
    % save the weights (importance index according to Stefanini et al.
    % 2020)
    % get the number of cells
    cell_number = size(train_set,2);
    % get the number of folds
    fold_number = size(af_deco.Trained,1);
    % allocate memory for the weights across folds
    fold_weights = zeros(fold_number,cell_number);
    % for all the folds
    for folds = 1:fold_number
        % get the current fold
        current_fold = af_deco.Trained{folds};
        % allocate memory for the individual learner weights
        learner_number = size(current_fold.BinaryLearners,1);
        learner_weights = zeros(learner_number,cell_number);
        % for all the learners
        for learners = 1:learner_number
            % store the current learner weights
            learner_weights(learners,:) = current_fold.BinaryLearners{learners}.Beta;
        end
        % average the learner weights
        fold_weights(folds,:) = mean(abs(learner_weights),1);
        
    end   
    
    % average across folds and save
    redec_cell{redec,5} = mean(fold_weights,1)';
end

%concatenate the conf matrices
redec_conf = cat(3,redec_cell{:,1});
%and the diagonal fractions
redec_frac = horzcat(redec_cell{:,2});
%calculate the average conf mat
conf_mat = mean(redec_conf,3);

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
class_cell{4} = cat(2,redec_cell{:,3});
class_cell{5} = cat(2,redec_cell{:,4});
class_cell{6} = cat(2,redec_cell{:,5});