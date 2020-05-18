% calcium kernel evaluation
%% clean up
clearvars
close all force
load('paths.mat')
addpath(genpath(paths(1).main_path))
%% Load the files and define paths

% define whether to calculate kernels or use them
calculate_kernels = 0;
% define whether to use data or save an artificial kernel
data_kernel = 0;

%get the folder where the image files are
tar_path_all = uipickfiles('FilterSpec',...
    paths(1).stage2_path);

%get the number of experiments selected
num_exp = length(tar_path_all);

%define the list of labels to sort the files
% label_list = {'_thres.mat','_thresBEH.mat'};
label_list = {'_thres.mat'};

%get the program used
[~,name_whole,~] = fileparts(tar_path_all{1});
name_parts = strsplit(name_whole,'_');

%get the number of each type of data (the round is to avoid the nonscalar
%warning for the for loop)
num_data = round(num_exp./length(label_list));

%allocate memory for the different types of files
name_cell = cell(num_data,length(label_list));

%for the types of files
for f_type = 1:length(label_list)
    %get the coordinates of the file names
    name_map = ~cellfun(@isempty,strfind(tar_path_all,label_list{f_type}));
    %store them in the corresponding layer of the name cell
    name_cell(:,f_type) = tar_path_all(name_map);
end
%% Load traces and constants

% % load the main file
% main_file = load(name_cell{1});
% 
% % get the components
% conc_trace = main_file.conc_trace;
% stim_num = main_file.stim_num2;
% time_num = main_file.time_num;

load(name_cell{1});
%% Calculate a calcium kernel based on the onset of each stimulus

if calculate_kernels == 1
    close all
    
    % allocate a structure to save the kernel info
    kernel_str = struct([]);
    %if data should be used to calculate the kernel
    if data_kernel == 1
        %run the function. it will output the average "kernel" of the fastest
        %decaying stimulus
        [kernel_str(1).kernel,fit_val]= kernel_calc_1(conc_trace,stim_num2,time_num,0,'mean');
        % get the name for saving
        [~,ori_name] = fileparts(name_cell{1});
        kernel_str(1).ori_name = ori_name(1:end-6);
        % get the parameters
        kernel_str(1).kernel_length = length(kernel_str(1).kernel);
        kernel_str(1).kernel_scale = fit_val(1,1);
        kernel_str(1).kernel_tau = fit_val(2,1);

    else
        % define the length
        kernel_str(1).kernel_length = 20;
        % define the time constant
        kernel_str(1).kernel_tau = 1.2;
        % define the scale
        kernel_str(1).kernel_scale = 1;
        %provide code to calculate a custom kernel
        kernel_str(1).kernel = kernel_str(1).kernel_scale.*exp(...
            -(1:kernel_str(1).kernel_length).*kernel_str(1).kernel_tau);
        kernel_str(1).ori_name = 'Artificial';

    end

    % %use when approximating the kernel
    % conv_kernel = kernel_out;

    figure
    plot(kernel_str(1).kernel)
    
    %and save the kernel
    %define the save path
    save_path = paths.kernel_path;

    % assemble the file name
    save_clu = strjoin({kernel_str(1).ori_name,'len',num2str(kernel_str(1).kernel_length),...
        'scale',num2str(kernel_str(1).kernel_scale),...
        'tau',num2str(kernel_str(1).kernel_tau),'kernel.mat'},'_');

        save(fullfile(save_path,save_clu),'kernel_str')

end
%% Load kernels for data processing (axonal vs soma)

if calculate_kernels == 0
    %load the kernels
    %define the load path
    load_path = paths.kernel_path;
    
    % select the kernels to use
    kernel_names = uipickfiles('FilterSpec',load_path);
    % get the number of kernels
    num_kernels = size(kernel_names,2);
    % allocate memory for the kernels
    kernel_cell = cell(num_kernels,1);
    % load them
    for kernels = 1:num_kernels
        % load the kernel
        temp_kernel = load(strcat(kernel_names{kernels}));
        kernel_cell{kernels} = temp_kernel.kernel_str;
    end
    % concatenate the structures
    kernel_str = vertcat(kernel_cell{:});
    %% Plot the kernels
    close all
    % define the fontsize
    fontsize = 10;
    %define the time vector
    time_vec = 0.5:0.5:10;%for imaging at 2 Hz with a 20 frame kernel
    figure
    % for all the kernels
    for kernels = 1:num_kernels
        plot(time_vec,kernel_str(kernels).kernel);
        hold on
    end
        
    % plot(time_vec,deconv_kernel)
    % hold('on')
    % plot(time_vec,reconv_kernel)
    % plot(time_vec,conv_kernel)
    title('Approximate calcium kernels','FontSize',fontsize)
    xlabel('Time','FontSize',fontsize)
    ylabel('A.U.','FontSize',fontsize)
    legend({kernel_str.ori_name},'Interpreter','None')
    % legend({'Syn-GC6s','H2B-GC6s','Delay-Syn-GC6s'})
    set(gca,'FontSize',fontsize)
    %% Use the kernels for modifying the data before clustering
    close all
    
    % define the kernels to use for deconvolution and convolution
    deconv_kernel = kernel_str(contains({kernel_str.ori_name},'p17b_syngc6s')).kernel;
    conv_kernel = kernel_str(contains({kernel_str.ori_name},'Artificial')).kernel;
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
            
            [conv_kernel]= kernel_calc_1(conv_trace,stim_num2,time_num,0,'mean');


            figure
            imagesc(conv_trace)

            figure
            tar_trace = 1;
            plot(conc_trace(tar_trace,:),'r')
            hold('on')
            plot(conv_trace(tar_trace,:),'b')
            %replace the old conc_trace matrix
            conc_trace = conv_trace;
            
            figure
            plot(time_vec,conv_kernel);
            hold on
            plot(time_vec,kernel_str(contains({kernel_str.ori_name},'p17b_syngc6s')).kernel)
            plot(time_vec,kernel_str(contains({kernel_str.ori_name},'p17b_gc6s')).kernel)
            legend({'Convolved','Original','Target'})

    end
    %% save the modified dataset
    
    %define the path for saving the concatenated files
    thres_path = paths(1).stage2_path;
    %define the saving path
    %     [thres_name,thres_fpath] = uiputfile(thres_path);
    
%     thres_name = strcat(fname{3},'_',fname{2});
    % get the tau of the artificial kernel
    tau = kernel_str(contains({kernel_str.ori_name},'Artificial')).kernel_tau;
    [~,thres_name] = fileparts(name_cell{1});
    thres_name = thres_name(1:end-6);
    %save the corrected variables into a new file, already grouped
    save_name = strcat(thres_name,'_conv_tau_',num2str(tau),'_thres.mat');
    save(fullfile(thres_path,save_name),'conc_trace','fish_ori','cat_stack_all',...
        'cat_seed_all','cat_z_all','time_num','stim_num2','col_out','snr_mat','cat_anatomy_all','cat_reps')
end