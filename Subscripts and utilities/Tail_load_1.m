function [color_mat,idpos_mat,stim_mat,curv_vec,defsum_vec,time_vec] = Tail_load_1(tar_path)
%% Load the target file
%format the path from the parent
tar_file = strcat(tar_path,'.tdms');
%read the tdms file
[ConvertedData,~,ChanNames,~,~] = convertTDMS(0,tar_file);
%% Put the data on a cell array for easier access

%get the number of channels
chan_num = length(ChanNames{1});

%allocate memory for the data
data_cell = cell(chan_num,1);

%also get the names of the data in the original structure in a cell array
ori_names = {ConvertedData.Data.MeasuredData.Name};

%for all the channels
for chans = 1:chan_num
    %find which field corresponds to the data on the list of names
    field_vec = strcmp(ori_names,ChanNames{1}(chans));
    %store the data in the main cell
    data_cell{chans} = ConvertedData.Data.MeasuredData(field_vec).Data;
end
%% Decode the stimulus info

%consider that the setup starts recording as soon as the software is
%started, but the stimulus information comes in only after the scope starts
%saving. Hence I'll trim the tail trajectory files to match the color files
%in size since those were recorded frame by frame

% Order of the cells is:
% 1 Color 1
% 2 Color 2
% 3 ID Pos
% 4 Stim_info
% 5 Curvature
% 6 DeflectionSum
% 7 TailX
% 8 TailY
% 9 Time

% close all

%define the frame number (for the camera, not the scope)
frame_num = size(data_cell{1},1);
%% Decode the color info

%output after the flipud
% 1 Red right
% 2 Green right
% 3 Blue right
% 4 UV right
% 5 Red left
% 6 Green left
% 7 Blue left
% 8 UV left

%allocate memory for the color information
color_mat = zeros(8,frame_num);

%for both color channels
for chans = 1:2
    switch chans
        case 1
            col_range = 1:4;
        case 2
            col_range = 5:8;
    end
    %decode the color channel
    color_temp = reshape(typecast(uint32(data_cell{chans}),'uint8'),4,frame_num);
%     color_temp = typecast(uint32(data_cell{chans}),'uint8');
    color_mat(col_range,:) = flipud(double(color_temp(1:4,:)));
end
%% Decode the ID Pos and stim info

%ADDENDUM, CHECK
idpos_vec = data_cell{3}(end-frame_num+1:end);

%decode the idpos and stim channel
% idpos_mat = reshape(typecast(uint32(data_cell{3}),'uint16'),2,frame_num);
idpos_mat = reshape(typecast(uint32(idpos_vec),'uint16'),2,frame_num);

%row 1: position of random dot
%row 2: stimulus (actual)
stim_mat = reshape(typecast(uint32(data_cell{4}),'uint8'),4,frame_num);
%row 1: z position
%row 2: rep
%row 3: frame
%row 4: stim from scope
%% Trim the tail tracking info and time vector

curv_vec = data_cell{5}(end-frame_num+1:end);
defsum_vec = data_cell{6}(end-frame_num+1:end);
time_vec = data_cell{9}(end-frame_num+1:end);
%% OFF Add and plot the tail coordinates and a movie
% close all
% 
% coord_mat = cat(3,reshape(data_cell{7},6,length(data_cell{7})/6)...
%     ,reshape(data_cell{8},6,length(data_cell{8})/6));
% 
% coord_mat = coord_mat(:,end-frame_num+1:end,:);
% 
% xlim_vec = [min(min(coord_mat(:,:,1))),max(max(coord_mat(:,:,1)))];
% ylim_vec = [min(min(coord_mat(:,:,2))),max(max(coord_mat(:,:,2)))];
% 
% % figure
% % for f = 1:100:frame_num
% %     plot(coord_mat(:,f,1),coord_mat(:,f,2),'*')
% %     if f == 1
% %         set(gca,'NextPlot','replacechildren',...
% %             'XLimMode','Manual','YLimMode','Manual',...
% %             'XLim',xlim_vec,'YLim',ylim_vec)
% %     end
% %     
% %     pause(0.01)
% %     
% % end
% 
% plot_mat = sqrt(coord_mat(:,:,1).^2 + coord_mat(:,:,2).^2);
% figure
% for t = 1:6
%     subplot(6,1,t)
% %     plot(plot_mat(t,:))
% % figure
%     plot(squeeze(coord_mat(t,:,1)))
%     hold('on')
%     plot(squeeze(coord_mat(t,:,2)))
%     set(gca,'XTick',[])
% end
% 
% %calculate the correlation between the second segment and all the other
% %ones
% %allocate memory to store the correlation
% corr_vec = zeros(4,1);
% %for the 3rd to the last segment
% for seg = 3:6
%     temp = corrcoef(coord_mat(seg,:,1),coord_mat(2,:,1));
%     corr_vec(seg-2) = temp(2,1);
% end
% figure
% plot(corr_vec)
% figure
% plot(std(squeeze(coord_mat(:,:,1)),0,2))
%% Exclude defined tale segments and re-calculate the sum angle

close all
%define the segments to exclude
exc_vec = logical([1 1 1 1 1 0]);
%load the tail points
coord_tempx = data_cell{7}(1:frame_num*6);
coord_tempy = data_cell{8}(1:frame_num*6);
coord_mat = cat(3,reshape(coord_tempx,6,length(coord_tempx)/6)...
    ,reshape(coord_tempy,6,length(coord_tempy)/6));
%trim the trace as with the other ones (because the setup starts recording
%before the onset of the experiments and it accumulates extra points)
coord_mat = coord_mat(:,end-frame_num+1:end,:);
%exclude the target points
excmat = coord_mat(exc_vec,:,:);
%extract the baseline
base_point = squeeze(coord_mat(1,1,:));
tip_point = [mode(squeeze(coord_mat(6,:,1)));mode(squeeze(coord_mat(6,:,2)))];
base_delta = diff([base_point,tip_point],1,2);
base_angle = (atan2(base_delta(1),base_delta(2)));
%extract the angle between each segment at each timepoint

%calculate the deltas between consecutive segments in x and y
delta_mat = diff(excmat,1,1); 

%now calculate the angles between the segments
angle_mat = atan2(delta_mat(:,:,1),delta_mat(:,:,2))-base_angle;

%and finally add the resulting angles after conversion to degrees
defsum_vec2 = sum(rad2deg(angle_mat),1)';
defsum_vec = defsum_vec2;

% figure
% plot(defsum_vec)
% hold('on')
% plot(defsum_vec2)
% histogram(delta_mat(:),100)
% set(gca,'YScale','log')
%% OFF Plot the output (for quick evaluation of the file)
% figure
% 
% subplot(2,1,1)
% plot(time_vec,curv_vec)
% title('Tail Curvature')
% subplot(2,1,2)
% plot(time_vec,defsum_vec)
% title('Tail Deflection')
%% OFF Plot the color information

% %define the color labels
% col_labels = 'rbgmrbgm';
% 
% figure
% for chans = 1:2
%     switch chans
%         case 1
%             sign_col = 1;
%             col_range = 1:4;
%         case 2
%             sign_col = -1;
%             col_range = 5:8;
%     end
%     for cols = col_range
%         plot(sign_col.*color_mat(cols,:),col_labels(cols))
%         hold('on')
%     end
% end