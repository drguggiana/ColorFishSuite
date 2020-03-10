function [region_data,num_regions] = region_split(conc_trace,regions,name,varargin)

% check if the flag is there
if size(varargin,1) == 0 || varargin{1} == 0
    combined = 0;
else
    combined = varargin{1};
end

% check combined
if combined == 0
    if contains(name,'syn')
        %define the region labels
        reg_label = {'N/A','AF4','AF5','AF6','AF7','AF8','AF9','AF10','All'};
        reg_map = [0 4 5 6 7 8 9 10];
    else
        %define the region labels
        reg_label = {'N/A','L-TcN','R-TcN','L-TcP','R-TcP','L-Cb','R-Cb','L-Hb','R-Hb','L-Pt','R-Pt','All'};
        reg_map = [0:10];
    end

    % turn NaNs into 0
    regions(isnan(regions)) = 0;
    % get a list of the regions present
    list_regions = unique(regions);
    % get the number of regions
    num_regions = length(list_regions)-1;
    % allocate memory for the regions
    region_data = cell(num_regions,3);

    % for all the regions (excluding the sero)
    for region = 1:num_regions
        % get the traces from this region and store
        region_data{region,1} = conc_trace(regions==list_regions(region+1),:,:);
        % store also the label
        region_data{region,2} = reg_label{reg_map==list_regions(region+1)};
        % and the index of the region traces
        region_data{region,3} = regions==list_regions(region+1);
    end

else
    if contains(name,{'Syn', 'syn'})
        %define the region labels
        reg_label = 'RGCs';
    else
        %define the region labels
        reg_label = 'Tectum';
    end
    
    region_data = {conc_trace,reg_label,ones(size(conc_trace,1),1)};
    num_regions = 1;
    
end