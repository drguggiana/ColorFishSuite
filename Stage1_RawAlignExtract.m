%% clean up
clearvars
close all
% addpath(genpath('E:\Behavioral data\Matlab\AF_proc\ColorFishSuite'))
addpath(genpath('E:\Behavioral data\Matlab\AF_proc\ColorFishSuite'))
%% Define miscellaneous constants

%maximum allowed area for a seed
thres_area = 150;
%minimum seed area
thres_minarea = 10;

skip_stim = [];
comb_vec = [];

% % maximum allowed area for a seed
% thres_area = 10;
% %minimum seed area
% thres_minarea = 6;

%define whether to save files
save_var = 1;
%define whether to activate the plots
plot_var = 0;
% %define whether to cluster or not
% clu_var = 0;
% %define whether or not to use the old error filtering
% err_var = 0;
%% Load the files and define paths

%get the folder where the image files are
% tar_path_all = uipickfiles('FilterSpec','F:\20160603_2Pbackup');
tar_path_all = uipickfiles('FilterSpec','J:\2P');


%get the number of experiments selected
num_exp = length(tar_path_all);

%for all the experiments
for exp_ind = 1:num_exp
    
    %get the path of the current experiment
    tar_path = tar_path_all{exp_ind};

    %define the path for saving the file with the BIC info
    [~,tar_name] = fileparts(tar_path);
%     bic_name = strcat('F:\20160603_2Pbackup\Proc_data\',tar_name);
%     bic_name2 = strcat('F:\20160603_2Pbackup\Proc_data\',tar_name,'_wBEHAV');
    %determine whether it's an RF file
    str_parts = strsplit(tar_name,'_');
    if strcmp(str_parts{3},'p18')||strcmp(str_parts{3},'p18proto')
        rf_var = 1;
        %define the amount of mode taking when aligning
        mode_wind = 4;
    else
        rf_var = 0;
        %define the amount of mode taking when aligning
        mode_wind = 5;
    end
    %% And parse them based on their filenames

    % [file_info,stim_num,rep_num,z_num,tar_files] = parser_1(tar_path);
    % [file_info,stim_num,rep_num,z_num,tar_files] = parser_2(tar_path);
    [file_info,stim_num,rep_num,z_num,tar_files] = parser(tar_path);

    %using the first file, find out the number of time points in the trace
    first_info = imfinfo(fullfile(tar_path,tar_files{1}));
    time_num = length(first_info);
    %% Coordinate the RF stimulus processing

    %load the stimulus colors
    [~,~,stim_list,~,~,col_ave] = Tail_proc(tar_path,rf_var);
    %also replace the file_info stim column with the actual stim one coming
    %from the tail tracking data
    file_info(:,3) = stim_list;
    % %if there are any RF stimuli in the program
    % if ~isempty(skip_stim)
    %     %get the number of frames excluded
    %     rf_num = length(skip_stim);
    %     %and calculate how many frames will there be for correlation
    %     corr_frames = stim_num*time_num;
    %     %and the total number of frames in the singlez stack
    %     total_frames = corr_frames + rf_num*rep_num*time_num;
    %     %for the seeding, determine the frames to skip due to the exclusion
    %     frame_vec = ones(total_frames,1);
    %     %for all the skipped stimuli
    %     for sk_stim = 1:rf_num
    %         %mark the positions to skip with 0
    %         frame_vec((skip_stim(sk_stim)-1)*time_num+1:skip_stim(sk_stim)*time_num) = 0;
    %     end
    %     %finally define the number of frames of the seeded matrix
    %     seed_frames = total_frames - rf_num*time_num;
    %     %redefine stim_num
    %     nstim_num = round(seed_frames./time_num);
    % else
        %define all the frame constants as the same number
        corr_frames = stim_num*time_num;
        total_frames = corr_frames;
        seed_frames = corr_frames;
        %the allowed frame vector
        frame_vec = ones(corr_frames,1);
        %and the stim number stays constant
        nstim_num = stim_num;
    % end
    %% Run the z section main for loop
    %allocate memory to store the seeds per z-section and their image
    seed_cell = cell(z_num,1);
    im_cell = cell(z_num,1);
    %also for the traces
    trace_cell = cell(z_num,1);
    % for the single reps
    trace_cell_reps = cell(z_num,1);
    %for the shifts
    shift_cell = cell(z_num,1);
    %for the behavioral interpolations
    behav_cell = cell(z_num,1);
    %for the averaged frames
    ave_frame = cell(z_num,1);
    %and for the bout counts per stim and z
    bout_count = cell(z_num,1);
    %and for the SNR
    snr_cell = cell(z_num,1);

    %define the pre and post periods
    pre_time = false(time_num,1);
    stim_time = false(time_num,1);
    post_time = false(time_num,1);

    %if its not at RF paradigm
    if rf_var == 0
        pre_time(6:0.25*time_num) = 1;
        stim_time(0.25*time_num+1:0.75*time_num) = 1;
        post_time(0.75*time_num+1:end) = 1;
    else %if it is
        %total number of frames per bar
        bar_frames = 6;
        %define the frame with the bar
        trigger_frame = 2;
        %for all the frames per bar
        for frames = 1:bar_frames
            if any(frames == trigger_frame)
                stim_time(trigger_frame:bar_frames:time_num) = 1;
            else
                pre_time(frames:bar_frames:time_num) = 1;
            end

        end

    end
    %allocate memory to store the correlation stack and the shuff one
    corr_stack = zeros(first_info(1).Height,first_info(1).Width,z_num);
    % corr_shuff = zeros(first_info(1).Height,first_info(1).Width,z_num);
    %allocate memory for a cutoff vector
    seed_cutoff = ones(z_num,1);

    %if it's the first iteration
    if exp_ind == 1
        %activate the parallel pool
        gcp = parpool;
    end
    %for all the z sections
    parfor z = 1:z_num
        %% Calculate dfof, compress replicates, accumulate and align
    %     profile on
        [singlez_mat,shift_cell{z},ave_frame{z},snr_mat,allreps_mat] = ...
            aligner(tar_path,tar_files,stim_num,rep_num,file_info,pre_time,z...
            ,skip_stim,mode_wind);
    %     profile off
    %     profile viewer
        %NaN the edges of the image to avoid spurious correlations
        singlez_mat(1:5,:,:) = NaN;
        singlez_mat(:,1:5,:) = NaN;
        singlez_mat(end-4:end,:,:) = NaN;
        singlez_mat(:,end-4:end,:) = NaN;
        % do the same with the single rep traces
        allreps_mat(1:5,:,:,:) = NaN;
        allreps_mat(:,1:5,:,:) = NaN;
        allreps_mat(end-4:end,:,:,:) = NaN;
        allreps_mat(:,end-4:end,:,:) = NaN;
        %% Run the correlation software
        %only in the frames that were averaged or moded. Don't do it on the
        %ones added at the end from the exclusion
        corr_stack(:,:,z) = Stack_process(singlez_mat(:,:,1:corr_frames),0);
        %% OFF Shuffle each pixel along the time dimension
    %     shuff_mat = Shuffle(singlez_mat,3);
    %     %run the correlation on the shuffled data
    %     corr_shuff(:,:,z) = Stack_process_1(shuff_mat,0);
        %% Determine the seed threshold based on the shuff stack

    %     for c = 0:0.001:1
    %         if (sum(sum(corr_stack(:,:,z)>c)) / sum(sum(corr_shuff(:,:,z)>c))) >= 70000
    %             seed_cutoff(z) = c;
    %             break
    %         end
    %     end
        %linearize the correlation frame for percentile 
        lin_corr = corr_stack(:,:,z);
    %     seed_cutoff(z) = prctile(reshape(corr_stack(:,:,z),...
    %         [size(corr_stack,1)*size(corr_stack,2),1]),90);
        seed_cutoff(z) = prctile(lin_corr(:),95);
        %print the resulting value
        fprintf(strcat('The cutoff for z =',num2str(z),' is:',num2str(seed_cutoff(z)),'\n'))
        %% Seeding function

        %define the seeding parameters

        %minimum correlation to start a seed
        thres_seed = seed_cutoff(z);
        %minimum correlation to expand a seed
        thres_nb = seed_cutoff(z);
    %     %maximum allowed area for a seed
    %     thres_area = 5;
    %     %minimum seed area
    %     thres_minarea = 3;


        [seed_cell{z},im_cell{z}] = Seeder(corr_stack(:,:,z),...
            thres_seed,thres_nb,thres_area,thres_minarea);
        %% Use the seeds to process all of the slices in this z-section for signals
        %get the seed number in this z section
        seed_num = size(seed_cell{z},2);

        %allocate memory for storing the traces
        seed_currz = zeros(seed_num,seed_frames);
        % and for the rep traces
        seed_currz_reps = zeros(seed_num,seed_frames,rep_num);
        %and for the snr output
        snr_currz = zeros(seed_num,stim_num);
        %initialize a frame counter
        frame_counter = 1;
        %for all the frames
        for frames = 1:total_frames
    %         fprintf(strcat('Current frame: ',num2str(frames),'\n'));
            %if the frame is on the skip stimuli
            if frame_vec(frames) == 0
                %skip the stimulus
                continue
            end
            %load the current frame
            curr_frame = singlez_mat(:,:,frames);
            
            %for all the seeds in this z section
            for seed = 1:seed_num
                %add the intensities for each frame and store in the seed
                seed_currz(seed,frame_counter) = mean(curr_frame(seed_cell{z}(seed).pxlist));
                % also the single rep intensities
                % for all the reps
                for reps = 1:rep_num
                    % load the reps frame
                    curr_frame_reps = squeeze(allreps_mat(:,:,reps,frames));
                    seed_currz_reps(seed,frame_counter,reps) = mean(curr_frame_reps(seed_cell{z}(seed).pxlist));
                end
            end
            %update the frame counter
            frame_counter = frame_counter + 1;
        end

        %also calculate the snr for each seed by averaging across the snr of
        %its voxels

        %for all the seeds in this z section
        for seed = 1:seed_num
            %for all the stimuli
            for stim = 1:stim_num
                %load the frame for the current stim
                curr_stim = snr_mat(:,:,stim);
                %calculate the average of the seed voxels
                snr_currz(seed,stim) = mean(curr_stim(seed_cell{z}(seed).pxlist));
            end
        end
        %store the reshaped version in the output cell
        trace_cell{z} = reshape(seed_currz,[seed_num,time_num,nstim_num]); 
        % store the single reps too
        trace_cell_reps{z} = reshape(seed_currz_reps,[seed_num,rep_num,time_num,nstim_num]);
        %and the snr of the seeds in this z
        snr_cell{z} = snr_currz;
    end

    %turn the average cell into a stack
    ave_stack = cat(3,ave_frame{:});
    %NaN the edges of the image since they are not used
    ave_stack(1:5,:,:) = NaN;
    ave_stack(:,1:5,:) = NaN;
    ave_stack(end-4:end,:,:) = NaN;
    ave_stack(:,end-4:end,:) = NaN;
    %if it's the last iteration
    if exp_ind == num_exp
        %delete the parallel pool (closing it)
        delete(gcp)
    end
    % error('Stop here')
    %% Concatenate across sections
    close all

    %Method 2: shuffle and 2 stds

    %concatenate the trace cell info
    all_trace = vertcat(trace_cell{:});
    all_trace_reps = vertcat(trace_cell_reps{:});

    %create a vector with the z of each seed
    z_seed = zeros(size(all_trace,1),1);
    %initialize a counter
    z_count = 1;
    %for all the zs
    for z = 1:z_num
        %get the number of seeds in this z
        z_seednum = size(trace_cell{z},1);
        %turn the corresponding positions to the z in the map vector
        z_seed(z_count:z_seednum+z_count-1) = z;
        %update the counter
        z_count = z_count + z_seednum;
    end

    
    %concatenate the behavior correlated to fluo info
    all_behav = horzcat(behav_cell{:});
    %also for the seed_cell
    seed_concat = horzcat(seed_cell{:})';

    %% Assemble the trace matrix (excluding RF stimuli and combining checkers and OMR versions)
    close all

    %define a vector with the stimuli to exclude
    exclude_vec = {21:nstim_num,skip_stim};

    %then exclude the RF traces from clustering and combine the target
    %stimuli
    [conc_trace,~,stim_num2,~,col_out,~,~,all_trace_reps] = ...
        exc_comb(all_trace,exclude_vec,comb_vec,[],col_ave,[],[],all_trace_reps);

    %get the number of seeds
    seed_num = size(conc_trace,1);

    figure
    imagesc(normr_1(conc_trace,0))
    % imagesc(log(abs(conc_trace)))
    % figure
    % imagesc(bout_signal)
    %% Get the anatomy information
    
    [anatomy_info,area_all] = anatomy_extractor(tar_path,[round(vertcat(seed_concat(:).centroid)),z_seed]);
    %% Save analysis output
    if save_var == 1

        %define the save path
        save_path = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\Stage1_extraction';
        %define the root of the save name
        [~,ori_name,~] = fileparts(tar_path);
    %     %if the clustering is on
    %     if clu_var == 1
    %         %save the clustering output
    %         save_clu = strcat(ori_name,'_clusters.mat');
    %         save(fullfile(save_path,save_clu),'pcs','GMModel','idx_clu','clu_num')
    %     end
        % save the trace cell extracted from the fluo data
        save_trace = strcat(ori_name,'_traces.mat');
        save(fullfile(save_path,save_trace),'conc_trace','trace_cell','seed_cell',...
            'shift_cell','ave_stack','seed_concat','z_seed','snr_cell','all_trace_reps','anatomy_info','area_all')
    %     %also save the seed cell, behav cell, shift cell, bout_count and bout stim
    %     save_behav = strcat(ori_name,'_behav.mat');
    %     save(fullfile(save_path,save_behav),'behav_cell','bout_count','bout_stim'...
    %         ,'bout_signal','bta_out','corr_out')
        %and save variables for plotting and such
        save_plot = strcat(ori_name,'_plot.mat');
        save(fullfile(save_path,save_plot),'time_num','stim_num2','col_out')
        %finally, write a tiff file with the average anatomy stack

        %for all the zs
        for z = 1:z_num
            %             subplot(round(sqrt(z_num)),ceil(sqrt(z_num)),z)
            %allocate memory for the seed locations
            C = uint16(ave_stack(:,:,z));
            save_fig = strcat(ori_name,'_anato','.tif');
            fig_full = fullfile(save_path,save_fig);
            if z ==1
                imwrite(C,fig_full,'tif','WriteMode','overwrite')
            else
                imwrite(C,fig_full,'tif','WriteMode','append')
            end
        end
    end
end
%% OFF Load analysis files
% %define the load path
% load_path = 'E:\Behavioral data\Matlab\AF_proc\Clustering results';
% %provide the target loadding file name
% ori_name = '20160610_SynG6s_3';
% %save the clustering output
% load_clu = strcat(ori_name,'_clusters.mat');
% load(fullfile(load_path,load_clu),'pcs','GMModel','idx_clu','clu_num')
% 
% % save the trace cell extracted from the fluo data
% load_trace = strcat(ori_name,'_traces.mat');
% load(fullfile(load_path,load_trace),'trace_cell','seed_cell','shift_cell','ave_stack')
% %also save the seed cell, behav cell, shift cell, bout_count and bout stim
% load_trace = strcat(ori_name,'_behav.mat');
% load(fullfile(load_path,load_trace),'behav_cell','bout_count','bout_stim')