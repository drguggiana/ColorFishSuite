function [region_data,num_regions] = region_split(conc_trace,regions,name,varargin)

% check if the combine flag is there
if size(varargin,2) == 0 
    combined = 0;
else
    combined = varargin{1};
end
% check if the region selection flag is there
if size(varargin,2) > 1
    region_selection = varargin{2};
else
    region_selection = [];
end


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
% check if a subselection was provided
if ~isempty(region_selection)
    % modify the set of regions present based on the selection
    % for all the regions in the selection
    for region = 1:size(list_regions,1)
        % get the region name
        region_name = reg_label{list_regions(region)==reg_map};
        % if the region is not on the list, make it zero
        if all(~contains(region_selection,region_name))
            regions(regions==list_regions(region)) = 0;
        end
    end
    % get list_regions again
    list_regions = unique(regions);
end

% check combined
if combined == 1
    % edit the list of regions and the index vector to be the same region
    % for all detected
    list_regions = [0, 1];
    regions(regions>0) = 1;
    reg_map = [0, 1];
    if contains(name,{'Syn', 'syn'})
        %define the region labels
        reg_label = {'N/A','RGCs'};
    else
        %define the region labels
        reg_label = {'N/A','Tectum'};
    end

end

% get the number of regions
num_regions = length(list_regions)-1;
% allocate memory for the regions
region_data = cell(num_regions,3);

% for all the regions (excluding the zero)
for region = 1:num_regions
    % get the traces from this region and store
    region_data{region,1} = conc_trace(regions==list_regions(region+1),:,:);
    % store also the label
    region_data{region,2} = reg_label{reg_map==list_regions(region+1)};
    % and the index of the region traces
    region_data{region,3} = regions==list_regions(region+1);
end
 

 