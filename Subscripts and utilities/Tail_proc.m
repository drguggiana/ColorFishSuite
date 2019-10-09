function [tail_cell,exp_map,stim_list,time_vec,tail_vec,col_ave] = Tail_proc(tar_path,rf_var)
%% Load the tail data from the corresponding tdms file

[color_mat,idpos_mat,stim_mat,~,defsum_vec,time_vec] = Tail_load_1(tar_path);
%% Get the time vector, zero it and scale it to seconds (approximate since
%I'm not sure what's going on with the time stamps)
time_vec = (time_vec-time_vec(1));
time_diff = diff(time_vec);
%if there is 1 glitch in the time trace
if sum(time_diff<0)==1
    %if it's only one, correct it and move on
    neg_ind = find(time_diff<0);
    time_vec(neg_ind+1:end) = time_vec(neg_ind) + cumsum([mode(time_diff);time_diff(neg_ind+1:end)]);
elseif sum(time_diff<0) > 1
    %if it's more than one, kill the program and throw an error message
    error('Too many errors in the time trace')
end
time_vec = (1/200.7).*time_vec./mean(time_diff); %using an estimated frame rate
%% ID bouts

% %center the trace based on the mode (in case the original angle was
% %skewed)
% ori_center = defsum_vec-mode(defsum_vec);

%use a rolling mode
mode_wind = 1000;
sub_sam = 10;
roll_mode = arrayfun(@(i) mode(defsum_vec(1+i-(mode_wind/2):i+mode_wind/2))...
    ,mode_wind/2:sub_sam:length(defsum_vec)-(mode_wind/2));
ori_center = arrayfun(@(i) defsum_vec((mode_wind/2)+sub_sam*(i-1)-sub_sam/2:...
    (mode_wind/2)+sub_sam*(i-1)+sub_sam/2-1)-roll_mode(i),1:length(roll_mode),'UniformOutput',0);
ori_center = vertcat(ori_center{:});
start_mode = defsum_vec(1:(mode_wind/2)-sub_sam/2-1);
% end_mode = defsum_vec(end-(mode_wind/2)-sub_sam/2-1:end);
end_mode = defsum_vec((mode_wind/2)+sub_sam*(length(roll_mode)-1)+sub_sam/2:end);
ori_center = [mode(start_mode).*ones(length(start_mode),1)...
    ;ori_center;mode(end_mode).*ones(length(end_mode),1)];

% %use a high pass filter to get rid of illumination changes
% high_pass = designfilt('highpassiir','StopbandFrequency',0.1...
%     ,'PassbandFrequency',10,'SampleRate',200);
% ori_center = filter(high_pass,defsum_vec); 

% % eliminate the changes in baseline by subtracting the mode of the trace
% % every n frames
% wind_size = 1000;
% defsum_std = defsum_vec;
% defsum_std(abs(defsum_vec)>1*std(defsum_vec)) = 0;
% ori_roll = filtfilt(ones(1,wind_size)/wind_size,1,defsum_std);
% ori_center = defsum_vec - ori_roll;

%now low pass with a rolling average
wind_size = 3;
f_trace = filtfilt(ones(1,wind_size)/wind_size,1,ori_center);

%send the bandpassed trace as the tail_vec variable
tail_vec = f_trace;

%and then run a sliding standard deviation window
std_size = 6;
%calculate the standard deviation of the tail trace with a sliding window
std_trace = movingstd(f_trace,std_size);

%finally threshold to bring things to 0 except the bouts
std_ctt = 5;
std_thres = std_ctt.*mean(std_trace);
fstd_trace = std_trace.*(std_trace>std_thres);

%and generate a binary version of the trace above for finding the bout
%starts and ends
bin_trace = diff(std_trace>std_thres);
%get the starts and ends
b_start = find(bin_trace>0);
b_end = find(bin_trace<0);

%if the first end is before the first start, eliminate it
if b_end(1)<b_start(1)
    b_end = b_end(2:end);
end
%also, if the last start is after the last end, also kill it
if b_start(end)>b_end(end)
    b_start = b_start(1:end-1);
end
%% Bout calculations
%calculate bout isi and merge bouts with isi's below a threshold
isi_thres = 0.05;
%create sub vectors of the starts and ends for merging later
start_sub = b_start(2:end);
end_sub = b_end(1:end-1);
b_isi = time_vec(start_sub) - time_vec(end_sub);

b_start = [b_start(1);start_sub(b_isi>isi_thres)];
b_end = [end_sub(b_isi>isi_thres);b_end(end)];

%eliminate bouts that are too short (while calculating the bout durations)
dur_thres = 0.1; %in seconds
b_duration = time_vec(b_end) - time_vec(b_start);

b_start = b_start(b_duration>dur_thres);
b_end = b_end(b_duration>dur_thres);
%and update both the isi and dur vectors for plotting
b_duration = b_duration(b_duration>dur_thres);
start_sub = b_start(2:end);
end_sub = b_end(1:end-1);
b_isi = time_vec(start_sub) - time_vec(end_sub);
b_freq = 1./b_isi;

%get the number of bouts
n_bouts = length(b_start);
%calculate the swim vigor per bout and the bout direction (as taken from
%the resulting angle after adding positive and negative values)
%allocate memory for the vigor
b_vigor = zeros(n_bouts,1);
b_dir = zeros(n_bouts,2);

%also calculate a low passed version of the original trace for
%determination of bout parameters
wind_size = 5;
l_trace = filtfilt(ones(1,wind_size)/wind_size,1,defsum_vec-mode(defsum_vec));

%for all the bouts
for bout = 1:n_bouts
    %integrate the absolute angle for each bout for the vigor
    b_vigor(bout) = sum(abs(f_trace(b_start(bout):b_end(bout))));

    %and quantify the deflection of the first tail movement to detemine bout
    %direction (from Kuo-Hua's paper)
    
    %acquire the bout trace
    temp_trace = f_trace(b_start(bout):b_end(bout));    
%     temp_trace = l_trace(b_start(bout):b_end(bout));    
    
    %save the max amplitude
    [~,ind] = max(abs(temp_trace));
    b_dir(bout,1) = abs(temp_trace(ind));
    
%     %and also the summed tail curvature (signed)
%     b_dir(bout,2) = sum(temp_trace);
%     b_dir(bout,2) = sum(temp_trace(1:ind));
    
    %find the peaks
    [~,peak_trace] = findpeaks(abs(temp_trace));
    
%     %OPTION 1 save the first peak
%     b_dir(bout,2) = temp_trace(peak_trace(1));
    %OPTION 2 pick largest peak
    %if there are peaks detected
    if ~isempty(peak_trace)
        %if the bout is less than 3 peaks
        if length(peak_trace)<3
            %select however many peaks there are
            tar_peak = length(peak_trace);
        else %if more peaks, grab the first 3
            [~,tar_peak] = max(abs(temp_trace(peak_trace(1:3))));
        end
        %if the highest is not the first peak (most cases)
        if tar_peak > 1
            %if the difference between the highest peak and the previous one is
            %less than 50% of the highest peak
            if (abs(temp_trace(peak_trace(tar_peak))) - abs(temp_trace(peak_trace(tar_peak-1))))<...
                    0.5*abs(temp_trace(peak_trace(tar_peak)))
                %select the previous peak instead
                tar_peak = tar_peak - 1;
            end
        end
        %then store the peak deflection in the output vector
%         b_dir(bout,2) = temp_trace(peak_trace(tar_peak));
        b_dir(bout,2) = sum(temp_trace(1:peak_trace(tar_peak)));
    else %if empty, store a NaN and continue
        b_dir(bout,2) = NaN;
        continue
    end
end

%concatenate the starts and ends for export (in actual time)
tail_temp = [time_vec(b_start),time_vec(b_end),b_vigor,b_dir];
%% Plot traces and stages of analysis

close all
% plot_offset = 10000;
% plot_range = 1+plot_offset:16000+plot_offset;
plot_range = 1:length(std_trace);
figure
subplot(4,1,1)
plot(time_vec(plot_range),defsum_vec(plot_range))
hold('on')
% plot(time_vec(plot_range),ori_center(plot_range))
% plot(get(gca,'XLim'),[0 0],'.-')
% plot(get(gca,'XLim'),2.*[std(defsum_vec),std(defsum_vec)],'.-k')
% plot(get(gca,'XLim'),-2.*[std(defsum_vec),std(defsum_vec)],'.-k')
% plot(time_vec(plot_range),l_trace(plot_range))

% wind_size = 10;
% l_trace = filtfilt(ones(1,wind_size)/wind_size,1,defsum_vec);
% plot(time_vec(plot_range),defsum_vec(plot_range))
title('High-passed - Mode')

subplot(4,1,2)
plot(time_vec(plot_range),f_trace(plot_range))
title(strcat('Filtered trace - wind size: ',num2str(wind_size)))

subplot(4,1,3)
plot(time_vec(plot_range),std_trace(plot_range))
title(strcat('Std trace - std size: ',num2str(std_size)))

subplot(4,1,4)
plot(time_vec(plot_range),fstd_trace(plot_range))
title(strcat('Thres trace - std ctt: ',num2str(std_ctt)))

figure
plot(time_vec(plot_range),defsum_vec(plot_range))
hold('on')
plot(time_vec(plot_range),fstd_trace(plot_range))
plot(time_vec(b_start),fstd_trace(b_start),'g*')
plot(time_vec(b_end),fstd_trace(b_end),'r*')

figure
subplot(2,2,1)
histogram(b_duration,100)
title('Bout duration (s)')
subplot(2,2,2)
histogram(b_isi,100)
title('Bout isi (s)')
set(gca,'XScale','log','YScale','log')

subplot(2,2,3)
histogram(b_freq,1000)
set(gca,'XScale','log','YScale','log')
title('Bout frequency (Hz)')

figure
subplot(1,2,1)
histogram(b_dir(:,1),100)
xlabel('Bout Peak Deflection','Fontsize',10)
subplot(1,2,2)
histogram(b_dir(:,2),100)
xlabel('Bout Cumm Angle','Fontsize',10)
%% OFF Filter trace from the times where the scope saves frames

%the scope computer sends a constant 0 when it's initializing, therefore
%generating aberrant behavior in the stim computer at these transition
%periods. Hence I use the mode of the stim vector and the length of the
%intervals to filter out the aberrant traces

% close all
%correct the stimulus trace

% idpos_vec = double(idpos_mat(2,:));
% idpos_map = true(1,size(idpos_vec,2));
% idpos_bin = diff(idpos_vec)==0;
% idpos_cc = regionprops(idpos_bin,{'Area','PixelIdxList'});
% 
% %go through the components and NaN the short ones
% for comps = 1:length(idpos_cc)
%     %get the pixel list
%     pix_list = idpos_cc(comps).PixelIdxList;
%     %if the component is shorter than 1000
%     if idpos_cc(comps).Area <1000
%         %destroy the component
%         idpos_map([pix_list;pix_list(end)+1]) = 0;
%     else
%         %otherwise mark it as legit
%         idpos_map(pix_list(end)+1) = 1;
%     end
%         
% end
% 
% 
% % %filter all the file outputs with this vector
% % plot(idpos_bin)
% % hold('on')
% % plot(idpos_mat(2,:))
% % plot(uint16(idpos_map).*idpos_mat(2,:))
% 
% stim_mat = stim_mat(:,idpos_map)+1;
% time_vec = time_vec(idpos_map);
% idpos_mat = idpos_mat(:,idpos_map);
% idpos_mat(2,:) = idpos_mat(2,:)+1;
%% Generate a list of the stimuli at the fish in terms of stacks (not camera frames)
% idpos_map = find(diff(idpos_mat(2,:)));
% idpos_temp = find(diff(idpos_map)<1000);
% idpos_new = zeros(size(idpos_mat,2),1);
% for ids = 1:length(idpos_temp)
%     idpos_new( = 
% end

%adjust the ranges for the stim and frame indicators to start on 1
stim_mat = stim_mat + 1;
idpos_mat(2,:) = idpos_mat(2,:) + 1;

%get the set of available z positions, reps and stim
uni_z = double(unique(stim_mat(1,:)));
uni_rep = double(unique(stim_mat(2,:)));
uni_time = double(unique(stim_mat(3,:)));
uni_scstim = double(unique(stim_mat(4,:)));

%get the frame number (scope frames)
frame_num = length(uni_z)*length(uni_rep)*length(uni_time)*length(uni_scstim);

%allocate memory for the stim list
stim_list = zeros(length(uni_z)*length(uni_rep)*length(uni_scstim),1);

%initialize a counter for the presentations
pres_count = 1;
%for every z
for z = uni_z
    %get the current z coordinates
    z_vec = stim_mat(1,:)==z;
    %for every rep
    for rep = uni_rep
        %get the rep coordinates on this z
        rep_vec = stim_mat(2,:)==rep&z_vec;
        %correct the idpos vector from the randomization occuring on the
        %first frame of the presentation
        %for all the stimuli
        for stim = uni_scstim
            %get the frames corresponding to this stimulus
            stim_vec = stim_mat(4,:)==stim&rep_vec;
            %for the first and last scope frame in this rep and z
            for frames = min(uni_time)
                %get the coordinates of the frame
                tar_vec = stim_mat(3,:)==frames&stim_vec;
                %make the stimulus the mode of this frame (which is what it
                %should be)
                idpos_mat(2,tar_vec) = mode(idpos_mat(2,stim_vec));
            end
        end
        %find the unique values for the stim in this presentation
        [uni_val,~,~] = unique(idpos_mat(2,rep_vec),'stable');
        %find the mode of the stim codes for this presentation
        stim_list(pres_count:pres_count+max(uni_scstim)-1) = uni_val;
            
        %update the counter
        pres_count = pres_count + max(uni_scstim);
    end
end
%% Generate a vector that marks each individual frame

% close all
%calculate the boundaries of each frame based on the frames vector
frame_bound = [1,find(diff(double(stim_mat(3,:)+1))~=0),length(stim_mat)] ;

%create a vector that tags each frame with an individual number
frame_vec = zeros(size(stim_mat,1),1);

%for all of the frames
for frame = 1:frame_num
    frame_vec(frame_bound(frame)+1:frame_bound(frame+1)) = frame;
end
%allocate memory for every single frame in the scope and indicators of the
%z, rep, stim and frame
exp_map = zeros(frame_num,5);
%also allocate memory for storing the colors
col_map = zeros(frame_num,8);

%for all the frames
for frame = 1:frame_num
    %get the index to use
    c1 = frame_bound(frame)+5;
    %load the info from each frame into the map
    exp_map(frame,:) = [stim_mat(1:3,c1)',idpos_mat(:,c1)'];
    %also load the color info
    col_map(frame,:) = color_mat(:,c1)';
end

%generate a representative stimulus map for the color channels (for future
%plotting) 
%get the number of stimuli
stim_num = length(uni_scstim);
%get the number of time slices
time_num = length(uni_time);
%allocate memory for the representative color map
col_ave = zeros(stim_num,time_num,8);
%create a vector for selection of the first rep and z
first_repz = exp_map(:,1)==1&exp_map(:,2)==1;
%for all the stimuli
for stim = 1:stim_num
    %load the stim color info into the map
    col_ave(stim,:,:) = col_map(first_repz&exp_map(:,5)==stim,:);
end

% %plot the results
% figure
% for stim = 1:stim_num
%     subplot(round(sqrt(stim_num)),ceil(sqrt(stim_num)),stim)
%     imagesc(squeeze(col_ave(stim,:,:)))
% end


%filter the RF accounting for the rest periods

%if it's an RF paradigm
if rf_var == 0
    %generate a vector with the rest periods marked (frames 1-19 and 61-80 for
    %80 frame setups)
    one_pres = [zeros(20,1);ones(40,1);zeros(20,1)];
    %repeat the vector for the number of presentations
    rest_vec = repmat(one_pres,length(stim_list),1);
    %apply the vector to the position column of the map only
    exp_map(:,4) = exp_map(:,4).*rest_vec;
end
% %plots for debugging
% figure
% imagesc(exp_map)
% figure
% imagesc(col_map)
% figure
% plot(idpos_mat(1,:))
% figure
% plot(exp_map(:,4));
%% Generate a cell containing all of the bouts of the animal per frame (camera)

%separate the bouts into their respective frames (since the starts in
%tail_temp are noted in real time, gotta use the time_vec to bin them
%instead of just the frame_bounds in camera frames)
[~,~,bin] = histcounts(tail_temp(:,1),time_vec(frame_bound));

%allocate memory for a cell with the bouts
tail_cell = cell(frame_num,2);

%for all the scope frames
for frames = 1:frame_num
    %if there are bouts in this scope frame
    if ~isempty(bin==frames)
        tail_cell{frames,1} = tail_temp(bin==frames,:);
    end
    %record the start and end of the scope frame in camera frames
    tail_cell{frames,2} = [frame_bound(frames),...
        frame_bound(frames+1)-1];
end
%% OFF Plot all the bouts for each stimulus (during stimulation period)

% figure
% %get a vector marking the stimuli during each frame but excluding the rest
% %times
% stim_only = exp_map(:,5).*rest_vec;
% %initialize a counter for the cummulative bout angle
% cumm_ang = zeros(stim_num,1);
% %for all the stim
% for stim = 1:3%stim_num
% %     subplot(round(sqrt(stim_num)),ceil(sqrt(stim_num)),stim)
%     subplot(3,1,stim)
%     %get the bouts during the interval of interest
%     stim_frames = tail_cell(stim_only==stim,1);
%     %isolate the bouts only
%     stim_bouts = vertcat(stim_frames{:});
%     %get the number of bouts
%     frames_num = length(stim_bouts);
%     %for all the bouts
%     for frames = 1:frames_num
% %         %get the current stretch of angles
% %         curr_frame = defsum_vec(stim_bouts{frames}(1):stim_bouts{frames}(2));
%         %filter out the forward swims
%         if stim_bouts(frames,4)<60
%             continue
%         end
%         %find the indices of the bout in frames
%         b_on = find(stim_bouts(frames,1)==time_vec);
%         b_off = find(stim_bouts(frames,2)==time_vec);
%         %get the angle trace for the bout
%         curr_frame = defsum_vec(b_on:b_off);
%         %plot the bout
%         plot(curr_frame,'*-')
%         hold('on')
%         
%         %calculate the cumm turn angle
%         cumm_ang(stim) = cumm_ang(stim) + sum(curr_frame);%/frames_num;
%     end
% end
% 
% figure
% plot(cumm_ang,'*')
%% OFF Calculate pref indexes
% close all
% figure
% %get a vector marking the stimuli during each frame but excluding the rest
% %times
% stim_only = exp_map(:,5).*rest_vec;
% %allocate memory for the PI
% PI_vec = zeros(stim_num,1);
% 
% %for all the stim
% for stim = 1:stim_num
% %     subplot(round(sqrt(stim_num)),ceil(sqrt(stim_num)),stim)
%     subplot(stim_num,1,stim)
%     %get the bouts during the interval of interest
%     temp_bouts = tail_cell(stim_only==stim,1);
%     temp_bouts = vertcat(temp_bouts{:});
%     %filter the forward bouts out
%     temp_bouts(abs(temp_bouts(:,4))<60,:) = 0;
%     
%     plot(temp_bouts(:,4),'*')
%     hold('on')
%     plot(temp_bouts(:,5),'*')
%     plot(get(gca,'XLim'),[0 0],'k--')
%     
%     tar_col = 5;
%     %calculate the PI
% %     PI_vec(stim) = (sum(temp_bouts(:,tar_col)>0)-sum(temp_bouts(:,tar_col)<0))/...
% %         (sum(temp_bouts(:,tar_col)>0)+sum(temp_bouts(:,tar_col)<0));
%     turn_vec = sign(temp_bouts(:,5)).*temp_bouts(:,4);
%     PI_vec(stim) = (sum(temp_bouts(turn_vec>0,4))-sum(temp_bouts(turn_vec<0,4)))/...
%         (sum(temp_bouts(turn_vec>0,4))+sum(temp_bouts(turn_vec<0,4)));
%     
% end
% figure
% plot(PI_vec,'*')
%% OFF Debugging plots
% close all
% figure
% plot(idpos_mat(1,:))
% hold('on')
% plot(idpos_mat(2,:))
% plot(color_mat(1,:),'r')
% plot(color_mat(4,:),'m')
% 
% plot(-color_mat(5,:),'r')
% plot(-color_mat(8,:),'m')