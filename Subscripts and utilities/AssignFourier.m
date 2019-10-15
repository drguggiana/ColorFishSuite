function [four_out,qual_out] = AssignFourier(tar_trace,startFrame,endFrame,des_freq,frame_rate,aggregate,plot_flag,intval)
%     """
%         graph: The unit to which the timeseries fourier transform should be assigned
%         startFrame: The start-frame in the timeseries of the considered fragment
%         endFrame: The end-frame in the timeseries of the considered fragment
%         suffix: Gets added to the end of the attribute name to which transforms will be assigned
%         aggregate: If True AND the timeseries is a multiple of 2 period lengths
%         transforms will be computed on an average of beginning and end half of the timeseries
%     """
%anti-aliasing
filtered = gauss1d(tar_trace',frame_rate/4);%,frame_rate/2)
filtered = filtered(startFrame:endFrame);

% figure
% plot(filtered)
% filtered = tar_trace';
% filtered = filtered(startFrame:endFrame);
%TODO: Somehow make the following noise reduction more generally applicable...
%if the length of filtered is divisble by 2, break into two blocks and average for noise reduction
if aggregate == 1
    %Test if we can aggregate: Find the period length pl in frames. If the length of filtered
    %is a multiple of 2 period lengths (size = 2* N * pl), reshape and average across first
    %and second half to reduce noise in transform (while at the same time reducing resolution)
    pl = round(1 / des_freq * frame_rate);
    if mod(length(filtered)/pl,2)==0
        filtered = mean(reshape(filtered,2,floor(length(filtered)/2)),1);
    end
end
fft_trace = real(fft(filtered));
freqs = linspace(0,frame_rate,size(fft_trace,1));
mag = abs(fft_trace);
[~,ix] = min(abs(des_freq-freqs));%index of bin which contains our desired frequency
four_out = mag(ix);

qual_out = four_out./sum(mag);

%     setattr(graph,"fft_"+suffix,fft)
%     setattr(graph,"freqs_"+suffix,freqs)
%     setattr(graph,"fourier_ratio_"+suffix,mag[ix]/mag.sum())

% assignin('base','fft_trace',fft_trace)
% assignin('base','freqs',freqs)
% assignin('base','filtered',filtered)


if plot_flag == 1
    if intval == 1
        figure
    end
    plot(freqs,fft_trace)
    hold('on')
%     plot(ix,mag(ix)./sum(mag),'o')
end

function filt_out = gauss1d(tar_trace,sigma)

filt_size = sigma*4;
x = linspace(-filt_size / 2, filt_size / 2, filt_size);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize
% assignin('base','gaussFilter',gaussFilter)
filt_out = filter(gaussFilter,1,tar_trace);