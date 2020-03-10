function sort_index = sort_by_average(data, target_field)

% sort traces by max average stim response
% get the number of traces
trace_num = size(data.(target_field),1);
% calculate the average responses
% reshape the matrix
conc_reshape = reshape(data.(target_field),trace_num,data.time_num,data.stim_num);
ave_resp = squeeze(mean(conc_reshape,2));
% reshape the matrix back
[~,sort_index] = sortrows(ave_resp);