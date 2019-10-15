% Fourier analysis
function [four_cat, qual_cat] = fourier_extraction(conc_trace,time_num,stim_num2)
%% Fourier assignment

%define the parameters for the assignment
startFrame = [1,23,61];
endFrame = [20,60,80];
des_freq = 0.15625;
frame_rate = 1/0.800;
aggregate = 1;

all_trace = reshape(conc_trace,size(conc_trace,1),time_num,stim_num2);

% tar_trace = all_trace(10,:,1);

%allocate memory for the fourier peaks
four_peaks = zeros(size(all_trace,1),stim_num2,3);
%and for the quality
four_qual = zeros(size(all_trace,1),stim_num2,3);

%for all the stimuli
for stim = 1:stim_num2
    fprintf(strcat('Stim:',num2str(stim),'\n'))

    %for all the traces
    for tracevar = 1:size(all_trace,1)
        %for all the intervals
        for intval = 1:3
            [four_peaks(tracevar,stim,intval),four_qual(tracevar,stim,intval)] = ...
                AssignFourier(all_trace(tracevar,:,stim),startFrame(intval)...
                ,endFrame(intval),des_freq,frame_rate,aggregate,0,intval);
        end
    end
end

%reshape the output for easy handling
four_cat = reshape(four_peaks,size(all_trace,1),stim_num2*3);
qual_cat = reshape(four_qual,size(all_trace,1),stim_num2*3);
% figure
% imagesc(four_cat(s_ind,:))
%% OFF Instead of Fourier assignment, calculate a different metric

% %define the parameters for the assignment
% startFrame = 21;
% endFrame = 60;
% % des_freq = 0.15625;
% % frame_rate = 1/0.800;
% % aggregate = 1;
% 
% all_trace = reshape(conc_trace,size(conc_trace,1),time_num,stim_num2);
% 
% % tar_trace = all_trace(10,:,1);
% 
% %allocate memory for the full swing
% signal_swing = zeros(size(all_trace,1),stim_num2);
% 
% %for all the stimuli
% for stim = 1:stim_num2
%     fprintf(strcat('Stim:',num2str(stim),'\n'))
% 
%     %for all the traces
%     for tracevar = 1:size(all_trace,1)
%         %calculate the max swing as the difference between the 5th and 95th
%         %percentiles in the signal
%         signal_swing(tracevar,stim) = prctile(all_trace(tracevar,startFrame:endFrame,stim),88) - ...
%             prctile(all_trace(tracevar,startFrame:endFrame,stim),12);
%      
%     end
% end
% 
% close all
% figure
% subplot(1,8,1:7)
% imagesc(all_trace(:,:,1))
% subplot(1,8,8)
% imagesc(signal_swing(:,1))
% %% Cluster and plot the fourier results
% if four_var == 1
%     close all
% 
% %     %% test of the assignment with the regressor
% %     test_val = AssignFourier_1(wave_sine,21,60,des_freq,frame_rate,aggregate,1,2);
%     %% Select the top Fourier peaks for display on each cone
%     close all
%     %define the number to select
%     num_peaks = 20;
% 
%     %for all the stimuli
%     for stim = 1:stim_num2
%         figure
%         subplot(8,1,1:7)
%         %get the indexes of the top num_peaks peaks during stim period
%         [~,peak_ind] = sort(four_peaks(:,stim,2),'descend');
%         %set a counter for the height of the plot
%         height_c = 0;
%         %for all the desired peaks
%         for peak_c = 1:num_peaks
%     %         subplot(round(sqrt(num_peaks)),ceil(sqrt(num_peaks)),peak_c)
%             plot(conc_trace(peak_ind(peak_c),:)+height_c)
%             %update the height counter
%             height_c = height_c + 0.2;
%             hold('on')
% 
%             for stim = 1:stim_num2
%                 plot([stim*time_num,stim*time_num],get(gca,'YLim'),'k')
%             end
%             set(gca,'XLim',[0 time_num*stim_num2])
%     %         image_stim_1(conc_trace(peak_ind(1:num_peaks),:),time_num,stim_num2,col_out,2)
%         end
%         top_lim = get(gca,'YLim');
%         set(gca,'YLim',[-0.1 top_lim(2)],'YTick',[],'XTick',[])
%         subplot(8,1,8)
%         hold('on')
%         %plot the stimulus
%         col_cat = reshape(permute(col_out,[2 1 3]),size(col_out,1)*size(col_out,2),size(col_out,3));
%         col_cat = 0.1.*(col_cat-127)./max(col_cat(:)-127);
%         col_label = [1 0 0;0 1 0;0 0 1;1 0 1];
%         for side = 1:1
%             for chan = 1:4
%                 plot(col_cat(:,chan+4*(side-1))+height_c,'Color',col_label(chan,:))
%             end
%             height_c = height_c + 0.6;
%         end
%         set(gca,'YTick',[],'XTick',[],'XLim',[0 size(col_cat,1)])
%         set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
% 
%     end
%     %% Cluster the Fourier peaks
% 
%     clu_num = 0;
%     %define the vector of cluster numbers to try
%     clu_vec = [2 3 4 5 10 20 30 40 50 70];
%     replicates = 1;
%     [idx_clu,GMModelF,clu_num] = sPCA_GMM_cluster_3(four_cat(:,5:8),clu_num,[],[],[],[],bic_name,clu_vec,replicates);
%     %% Plot the Fourier clustering results
% 
%     close all
%     %sort the indexes and the fish IDs
%     [~,s_ind] = sort(idx_clu);
%     figure
%     % subplot(1,8,1:7)
%     imagesc(four_cat(s_ind,:))
%     ylabel('Trace #')
%     set(gca,'XTick',[],'YTick',[])
%     % subplot(1,8,8)
%     % imagesc(fish_ori(s_ind))
% 
%     %average within each cluster
%     %allocate memory for the averages
%     clu_aveF = zeros(clu_num,size(four_cat,2));
%     %and for the number of items in each cluster
%     clu_count = zeros(clu_num,1);
%     %for all the clusters
%     for clu = 1:clu_num
%         %average
%         clu_aveF(clu,:) = mean(four_cat(idx_clu==clu,:),1);
%         %and count
%         clu_count(clu) = sum(idx_clu==clu);
%     end
% 
% 
%     figure
%     imagesc(clu_aveF)
% 
% 
%     %average within each cluster
%     %allocate memory for the averages
%     clu_ave = zeros(clu_num,size(conc_trace,2));
%     %and for the number of items in each cluster
%     clu_count = zeros(clu_num,1);
%     %for all the clusters
%     for clu = 1:clu_num
%         %average
%         clu_ave(clu,:) = mean(conc_trace(idx_clu==clu,:),1);
%         %and count
%         clu_count(clu) = sum(idx_clu==clu);
%     end
% 
% 
%     figure
%     imagesc(clu_ave)
% 
%     figure
% 
%     bar(clu_count)
% end