function [loaded_structure] = load_clusters(in_path)
% load structures of data from this project

% if it's a cell, assume it's full paths
if iscell(in_path)
    tar_path_all = in_path;
else
    % detect whether the path has an extension or not
    [~,file_name] = fileparts(in_path);
    
    % if no extension, use the path to open uipickfiles
    if isempty(file_name)
        %get the folder where the image files are
        tar_path_all = uipickfiles('FilterSpec',in_path);
    else
        % encapsulate in a cell and use as full path
        tar_path_all = {in_path};
    end
end

%get the number of each type of data (the round is to avoid the nonscalar
%warning for the for loop)
num_data = length(tar_path_all);

% initialize a cell to hold the data
data_cell = cell(num_data,1);

% allocate memory for the field names
fields_cell = cell(num_data,1);

% for all the files
for files = 1:num_data
    
   % fill up the structure
   data_cell{files} = load(tar_path_all{files},'main_str');
   data_cell{files} = data_cell{files}.main_str;
   % store the field names
   fields_cell{files} = fieldnames(data_cell{files});
end

% remove fields that don't overlap

% get the min number of fields
[~,min_field] = min(cellfun(@length,fields_cell));
% run through the fields and erase the extra
for files = 1:num_data
    % if it's the min_field, skip
    if files == min_field
        continue
    end
    % get the difference between fields
    delete_fields = setdiff(fields_cell{files},fields_cell{min_field});
    % for all the target fields
    for fields = 1:length(delete_fields)
        % delete those fields
        data_cell{files} = rmfield(data_cell{files},delete_fields{fields});
    end
    
end

% concatenate the data for output
loaded_structure = cat(2,data_cell{:});