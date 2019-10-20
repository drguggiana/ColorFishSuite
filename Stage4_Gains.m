% gain analysis
clearvars
close all force
addpath(genpath('E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\'))
%% Load clusters
in_path = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\Stage3_cluster\';
% [clu_ave, idx_clu, clu_number, conc_trace, stim_num2, time_num, col_out, fish_ori, file_name] = load_clusters(in_path,{'_clusters.mat'});

main_str = load_clusters(in_path);

%% Get the useful variables from the structure

conc_trace = main_str(1).conc_trace;
time_num = main_str(1).time_num;
stim_num2 = main_str(1).stim_num;
col_out = main_str(1).col_out;
%% Extract the Fourier power at the stimulus frequency

[four_cat, qual_cat] = fourier_extraction(conc_trace,time_num,stim_num2);
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

% %append the fourier quality for clustering (both in rbgu)
% delta_clu = cat(2,delta_signals,four_qual(:,5:8));
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
%% Now determine the gain from each cone to each cell
% figure
% errorbar(1:cone_num,mean(delta_signals,1),std(delta_signals,0,1)./sqrt(11))
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
            delta_signals(:,stim) = four_cat(:,stim+4);
    
%     %use the full swing of the signal
%     delta_signals(:,stim) = signal_swing(:,stim);
    
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

% add the gain info to the structure
main_str(1).delta_norm = delta_norm;
%% Save the gains
save_path = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\Stage4_gain\';

%save the clustering output
save_clu = strcat(main_str(1).name,'_gains.mat');
% save(fullfile(save_path,save_clu),'idx_clu','clu_number','delta_norm',...
%     'conc_trace','stim_num2','time_num','col_out','fish_ori')

save(fullfile(save_path,save_clu),'main_str')