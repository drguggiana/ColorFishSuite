function [conc_trace, stim_num2, col_out, fname] = dataset_fixes(fname,conc_trace,stim_num2,col_out,time_num)

% get the number of traces
trace_num = size(conc_trace,1);
% reshape the incoming data to isolate stimuli in a single dimension
conc_trace = reshape(conc_trace,trace_num,time_num,stim_num2,[]);
% get the number of reps
rep_num = size(conc_trace,4);
%if it's a gain program
if any(strcmp(fname{3},{'p17b','p21'}))==1
    %flip the middle stimuli in conc_trace to match wavelength
    conc_trace = conc_trace(:,:,[1,3,2,4],:);
elseif strcmp(fname{3},'p20')==1
    %flip the middle stimuli in conc_trace to match wavelength
    conc_trace = conc_trace(:,:,[1,3,2,4,5],:);
elseif strcmp(fname{3},'p22')==1
    %flip the middle stimuli in conc_trace to match wavelength
    conc_trace = conc_trace(:,:,[1,3,2,4,5,7,6,8],:);
elseif strcmp(fname{3}, 'p6')==1
    % if it's p6, remove the stimuli not to be used in the rest of the
    % analysis when comparing to p8
 
    %stimuli to exclude
    keep_stim = [5 7 9 11 22 24]; % excludes all but 6 stim out of 24
    
    % take only the stimuli in the list
    conc_trace = conc_trace(:,:,keep_stim,:);

    % figure
    % imagesc(conc_trace)
    % figure
    % imagesc(exc_trace)

    %rewrite the stim_num2 variable
    stim_num2 = size(conc_trace,3);
    %also modify the color code variable
    if ~isempty(col_out)
        col_out = col_out(keep_stim,:,:);
    end
    % also change the program variable to p8 for saving
    fname{3} = 'p8';
        
end

% reshape the input data back to 2D
conc_trace = reshape(conc_trace,trace_num,time_num*stim_num2*rep_num);