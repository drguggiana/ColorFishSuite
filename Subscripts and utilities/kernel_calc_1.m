function kernel_out = kernel_calc_1(conc_trace,stim_num2,time_num,plot_flag,type_var)

% get the number of clusters
trace_num = size(conc_trace,1);
%allocate memory for the cluster average reshaped by trials
trace_form = zeros(trace_num,time_num,stim_num2);

%define the color labels
color_label = [1 0 0;0 1 0;0 0 1;1 0 1];
%define the frame_rate in Hz
fr = 0.5; 

%for all the clusters
for traces = 1:trace_num
    trace_form(traces,:,:) = reshape(conc_trace(traces,:),time_num,stim_num2);
end

%define the extent of the initial pulse
pulse_ext = 20;

%grab the beginning of every trial
trial_start = trace_form(:,1:pulse_ext,:);


%average across stimuli and traces
overall_ave = mean(mean(trial_start,3),1);
%averages across stimuli
stim_ave = squeeze(mean(trial_start,1));

%calculate decay constant
%grab only the decay from the max
decay_mat = stim_ave(5:end,:);
%remove the curve offsets
decay_mat = bsxfun(@minus,decay_mat,min(decay_mat,[],1));
%create the exponential fit type
exp_fit = fittype(@(a,b,x) a*exp(-b*x));

%allocate memory to store the fit parameters
fit_val = zeros(2,stim_num2);
%and the model itself
fit_model = cell(stim_num2,1);
%create the x vector
x_vec = 0:fr:size(decay_mat,1)*fr - fr;
%for all the stimuli
for stim = 1:stim_num2
    %fit an exponential to the decay
    temp_fit = fit(x_vec',decay_mat(:,stim),exp_fit,'StartPoint',[0.01 0.1]);
    %get the coefficient values
    fit_val(:,stim) = coeffvalues(temp_fit);
    %and store the model
    fit_model{stim} = temp_fit;
    
end

if plot_flag
    %plot the average
    %define the time axis based on the frame rate
    time_x = 0:fr:pulse_ext*fr-fr;
    figure
    plot(time_x,overall_ave)
    figure
    %for all the stimuli
    for stim = 1:stim_num2
    %     subplot(2,2,stim)
        plot(time_x,stim_ave(:,stim),'Color',color_label(stim,:))
        hold('on')

    end

    figure
    %for all the stimuli
    for stim = 1:stim_num2
        %also plot the trace and the fit
        subplot(2,2,stim)
        plot(x_vec,decay_mat(:,stim),'o','Color',color_label(stim,:))
        hold('on')
        plot(fit_model{stim})
    end
    %plot the parameters
    figure
    %for all the parameters
    for paras = 1:size(fit_val,1)
        subplot(2,2,paras)
        plot(1:stim_num2,fit_val(paras,:),'o')
    end
end

%pick the type of kernel to output
switch type_var
    case 'min'
        %grab the slowest decay and output that as a kernel
        [~,max_decay] = min(fit_val(2,:));
        kernel_out = normr_1(stim_ave(:,max_decay),1);
    case 'max'
        %grab the fastest decay and output that as a kernel
        [~,max_decay] = max(fit_val(2,:));
        kernel_out = normr_1(stim_ave(:,max_decay),1);
    case 'mean'
        %grab the mean decay and output that as a kernel
        kernel_out = normr_1(overall_ave,1);
    case 'int1'
        %grab the mean decay, normalize between 0 and integrated area 1
        kernel_out = normr_1(overall_ave,1)./sum(normr_1(overall_ave,1));
    
end
