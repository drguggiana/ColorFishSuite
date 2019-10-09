% calcium kernel evaluation

%% Calculate a calcium kernel based on the onset of each stimulus

close all

%define whether to calculate from data or input a custom kernel
data_kern = 0;
%if data should be used to calculate the kernel
if data_kern == 1
    %run the function. it will output the average "kernel" of the fastest
    %decaying stimulus
    kernel_out = kernel_calc_1(conc_trace,stim_num2,time_num,0,'mean');
else
    %provide code to calculate a custom kernel
    kernel_out = exp(-(1:20).*0.9);
end

% %use when approximating the kernel
% conv_kernel = kernel_out;

figure
plot(kernel_out)
%and save the kernel
%define the save path
save_path = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\ConvolutionKernels';

save_var = 1;
if save_var == 1

    %get the root of the save name
    [ori_name,~] = uiputfile(strcat(save_path,'*.*'));
    %save the clustering output
    save_clu = strcat(ori_name,'_kernel.mat');
    save(fullfile(save_path,save_clu),'kernel_out',...
        'stim_num2','time_num','name_cell')
end
%% Load kernels for data processing (axonal vs soma)

%load the kernels
%define the load path
load_path = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\ConvolutionKernels\';

%load the axonal kernel
deconv_name = uigetfile(strcat(load_path,'*.mat'),'Select the axonal kernel');
deconv_kernel = load(strcat(load_path,deconv_name),'kernel_out');
deconv_kernel = deconv_kernel.kernel_out;

%load the nuclear kernel
reconv_name = uigetfile(strcat(load_path,'*.mat'),'Select the axonal kernel');
reconv_kernel = load(strcat(load_path,reconv_name),'kernel_out');
reconv_kernel = reconv_kernel.kernel_out;

%load a convolution only kernel
conv_name = uigetfile(strcat(load_path,'*.mat'),'Select the axonal kernel');
conv_kernel = load(strcat(load_path,conv_name),'kernel_out');
conv_kernel = conv_kernel.kernel_out;
%% Plot the kernels
close all
%define the time vector
time_vec = 0.5:0.5:10;%for imaging at 2 Hz with a 20 frame kernel
figure
plot(time_vec,deconv_kernel)
hold('on')
plot(time_vec,reconv_kernel)
plot(time_vec,conv_kernel)
title('Approximate calcium kernels','FontSize',20)
xlabel('Time','FontSize',20)
ylabel('A.U.','FontSize',20)
legend({'Syn-GC6s','H2B-GC6s','Delay-Syn-GC6s'})
set(gca,'FontSize',20)
%% Use the kernels for modifying the data before clustering
close all

%define whether to reconv or not
reconv_var = 2;
%get the number of total time frames (per trace)
trace_time = size(conc_trace,2);
%deconvolve the data with the axonal kernel
%allocate memory for the deconvolved traces
deconv_trace = zeros(size(conc_trace));
%get the number of traces
trace_num = size(deconv_trace,1);
%for all the traces
for traces = 1:trace_num
    [deconv_trace(traces,:),~] = deconv([conc_trace(traces,:),zeros(1,length(deconv_kernel)-1)],deconv_kernel);
end

%if reconv is desired
switch reconv_var
    case 1
    %allocate memory for the reconvolved traces
    reconv_trace = zeros(size(deconv_trace));
    %now convolve each trace with the nuclear kernel
    %for all the traces
    for traces = 1:trace_num
        temp = conv(deconv_trace(traces,:),reconv_kernel);
        reconv_trace(traces,:) = temp(1:trace_time);
    end

    %plot both data sets
    figure
    imagesc(conc_trace)
    figure
    imagesc(deconv_trace)
    figure
    imagesc(reconv_trace)

    figure
    tar_trace = 1;
    plot(conc_trace(tar_trace,:),'r')
    hold('on')
    plot(deconv_trace(tar_trace,:),'b')
    plot(reconv_trace(tar_trace,:),'m')

    %save the reconvolved trace for clustering
%     conc_trace = reconv_trace;

    case 0
    %if no reconv, just load the deconv trace for clustering
%     conc_trace = deconv_trace;
    case 2
        %convolution only
        %allocate memory for the reconvolved traces
        conv_trace = zeros(size(conc_trace));
        %now convolve each trace with the nuclear kernel
        %for all the traces
        for traces = 1:trace_num
            temp = conv(conc_trace(traces,:),conv_kernel);
            conv_trace(traces,:) = temp(1:trace_time);
        end
        
        figure
        imagesc(conv_trace)
        
        figure
        tar_trace = 1;
        plot(conc_trace(tar_trace,:),'r')
        hold('on')
        plot(conv_trace(tar_trace,:),'b')
        %replace the old conc_trace matrix
        conc_trace = conv_trace;

end