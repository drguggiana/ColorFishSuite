% Fourier analysis
function [four_cat, qual_cat] = fourier_extraction(conc_trace,time_num,stim_num2,des_freq,frame_rate)
%% Fourier assignment

%define the parameters for the assignment
startFrame = [1,23,61];
endFrame = [20,60,80];
% des_freq = 0.15625;
% frame_rate = 1/0.800;
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