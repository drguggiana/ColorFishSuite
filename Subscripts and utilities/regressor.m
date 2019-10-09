% regressor analysis

%% Regressor correlation

if reg_var == 1
    %define the regressor

    %define the number of excitation levels to be calculated
    num_lev = 40;
    %define the number of cycles
    wave_cyc = 5;

    %get the vector to calculate the sine with
    wave_vec = 2*pi()*wave_cyc*linspace(0,1,num_lev);
    %finally, calculate the sine and add the edges
    edge_vec = zeros(1,20);
    % wave_sine = sin(wave_vec);
    % wave_sine = [edge_vec,sin(wave_vec),edge_vec];
    wave_sine = [edge_vec,ones(1,num_lev),edge_vec];
    % wave_sine = [edge_vec,1,zeros(1,39),edge_vec];
    % wave_sine = ([edge_vec,ones(1,num_lev),edge_vec].*-1)+1;

    all_trace = reshape(conc_trace,size(conc_trace,1),time_num,stim_num2);

    %allocate memory for the fourier peaks
    corr_peaks = zeros(size(all_trace,1),stim_num2);

    %for all the stimuli
    for stim = 1:stim_num2
        fprintf(strcat('Stim:',num2str(stim),'\n'))

        %for all the traces
        for tracevar = 1:size(all_trace,1)
            corr_peaks(tracevar,stim) = ...
                max(abs(xcorr(wave_sine,all_trace(tracevar,:,stim))));
        end
    end

    figure
    imagesc(corr_peaks)
    % corr_cat = reshape(corr_peaks,size(all_trace,1),stim_num2);
    %% Select the top Corr peaks for display on each cone
    close all
    %define the number to select
    num_peaks = 20;

    %for all the stimuli
    for stim = 1:stim_num2
        figure
        subplot(8,1,1:7)
        %get the indexes of the top num_peaks peaks during stim period
        [~,peak_ind] = sort(corr_peaks(:,stim),'descend');
        %set a counter for the height of the plot
        height_c = 0;
        %for all the desired peaks
        for peak_c = 1:num_peaks
    %         subplot(round(sqrt(num_peaks)),ceil(sqrt(num_peaks)),peak_c)
            plot(conc_trace(peak_ind(peak_c),:)+height_c)
            %update the height counter
            height_c = height_c + 0.2;
            hold('on')

            for stim = 1:stim_num2
                plot([stim*time_num,stim*time_num],get(gca,'YLim'),'k')
            end
            set(gca,'XLim',[0 time_num*stim_num2])
    %         image_stim_1(conc_trace(peak_ind(1:num_peaks),:),time_num,stim_num2,col_out,2)
        end
        top_lim = get(gca,'YLim');
        set(gca,'YLim',[-0.1 top_lim(2)],'YTick',[],'XTick',[])
        subplot(8,1,8)
        hold('on')
        %plot the stimulus
        col_cat = reshape(permute(col_out,[2 1 3]),size(col_out,1)*size(col_out,2),size(col_out,3));
        col_cat = 0.1.*(col_cat-127)./max(col_cat(:)-127);
        col_label = [1 0 0;0 1 0;0 0 1;1 0 1];
        for side = 1:1
            for chan = 1:4
                plot(col_cat(:,chan+4*(side-1))+height_c,'Color',col_label(chan,:))
            end
            height_c = height_c + 0.2;
        end
        set(gca,'YTick',[],'XTick',[],'XLim',[0 size(col_cat,1)])
        set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

    end
end
