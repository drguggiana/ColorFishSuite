function [conc_trace,exc_trace,stim_num2,bout_signal,col_out,bta_out,corr_out] = ...
    exc_comb1(all_trace,exclude_vec,comb_vec,bout_stim,col_ave,bta_behav,corr_behav)

% %concatenate the cell entries to get all the seeds into a 3D matrix
% all_trace = vertcat(trace_cell{:});
%get the stim_num
stim_num = size(all_trace,3);
%if there are RF stimuli
if ~isempty(comb_vec)
    %this function excludes the traces in exclude_vec from the matrix to be
    %clustered and averages the ones in comb_vec together
    
   
    %create a matrix with the exclusions only and one without them
    exc_trace = all_trace(:,:,exclude_vec{1});
    %generate a vector without the exclusions
    exc_temp = ones(stim_num,1);
    exc_temp(exclude_vec{1}) = 0;
    nonexc_trace = all_trace(:,:,logical(exc_temp));
    
    %get the stim num from the non-trace data
    stim_num_nt = size(bout_stim,1);
    %create a sperate exclusion vector for the data coming from the
    %behavior side
    exc_temp2 = ones(stim_num_nt,1);
    exc_temp2(exclude_vec{2}) = 0;
    %also exclude the stimuli from the behavior traces
    nonexc_beh = bout_stim(logical(exc_temp2),:,:);
    %the color matrix
    nonexc_col = col_ave(logical(exc_temp2),:,:);
    %the btas
    point_num = size(bta_behav,2)./stim_num_nt;
    nonexc_bta = reshape(bta_behav',point_num,stim_num_nt,size(bta_behav,1));
    nonexc_bta = nonexc_bta(:,logical(exc_temp2),:);
    %and the corr coeffs
    nonexc_corr = corr_behav(logical(exc_temp2),:);
    %define the new stim_num
    stim_num_ex = size(nonexc_trace,3);

    %get the number of combine events
    comb_num = size(comb_vec,1);
    %get the number of noncombined stimuli
    noncomb_num = stim_num_ex - numel(comb_vec);
    %allocate memory to store the new matrix with the combined traces
    comb_trace = zeros(size(nonexc_trace,1),size(nonexc_trace,2),size(nonexc_trace,3)-comb_num);
    %also allocate memory for the combined behavior data
    comb_beh = zeros(size(nonexc_beh,1)-comb_num,size(nonexc_beh,2),size(nonexc_beh,3));
    %the color data
    comb_col = zeros(size(nonexc_beh,1)-comb_num,size(col_ave,2),size(col_ave,3));
    %the bta data
    comb_bta = zeros(size(nonexc_bta,1),size(nonexc_bta,2)-comb_num,size(nonexc_bta,3));
    %and the corr coeffs
    comb_corr = zeros(size(nonexc_corr,1)-comb_num,size(nonexc_corr,2));
    
    %get a matrix with all the stim that are not being combined
    noncomb_vec = ones(stim_num_ex,1);
    noncomb_vec(comb_vec(:)) = 0;
    comb_trace(:,:,1:noncomb_num) = nonexc_trace(:,:,logical(noncomb_vec));
    %do the same with the behavior traces
    comb_beh(1:noncomb_num,:,:) = nonexc_beh(logical(noncomb_vec),:,:);
    %with the color matrix
    comb_col(1:noncomb_num,:,:) = nonexc_col(logical(noncomb_vec),:,:);
    %with the btas
    comb_bta(:,1:noncomb_num,:) = nonexc_bta(:,logical(noncomb_vec),:);
    %and the corr coeffs
    comb_corr(1:noncomb_num,:) = nonexc_corr(logical(noncomb_vec),:);
    
    %this will have then PT CQ, then DF, OMR0 and OMR90

    %for all the combination cases
    for comb_case = 1:comb_num
        %load the combination in the matrix
        comb_trace(:,:,noncomb_num+comb_case) = mean(nonexc_trace(:,:,comb_vec(comb_case,:)),3);
        %in the behavior matrix
        comb_beh(noncomb_num+comb_case,:,:) = mean(nonexc_beh(comb_vec(comb_case,:),:,:),1);
        %in the color matrix
        comb_col(noncomb_num+comb_case,:,:) = mean(nonexc_col(comb_vec(comb_case,:),:,:),1);
        %in the btas
        comb_bta(:,noncomb_num+comb_case,:) = mean(nonexc_bta(:,comb_vec(comb_case,:),:),2);
        %and in the corr coeffs
        comb_corr(noncomb_num+comb_case,:) = mean(nonexc_corr(comb_vec(comb_case,:),:),1);
    end
    %define the new stim_num
    stim_num2 = size(comb_trace,3);
    %finally concatenate across time to get single traces
    conc_trace = reshape(comb_trace,[size(comb_trace,1),size(comb_trace,2)*size(comb_trace,3)]);
    %and send out the behavior trace
    bout_signal = comb_beh;
    %the color matrix
    col_out = comb_col;
    %the btas
    bta_out = reshape(comb_bta,stim_num2*point_num,size(comb_bta,3))';
    %and the corr coeffs
    corr_out = comb_corr;
else
    conc_trace = reshape(all_trace,size(all_trace,1),size(all_trace,2)*size(all_trace,3));
    bout_signal = bout_stim;
    stim_num2 = stim_num;
    exc_trace = [];
    col_out = col_ave;
    bta_out = bta_behav;
    corr_out = corr_behav;
end