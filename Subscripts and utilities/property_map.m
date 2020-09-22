function maps_cell = property_map(coord,indexes,registered_anatomy,property,stim_num,full_ref_dim,varargin)

% take the given seed map, or use a default
if length(varargin) >= 1
    seed_map = varargin{1};
else
    seed_map = 1:size(property,1);
end

% if given, use a threshold
if length(varargin) >= 2
    threshold = varargin{2};
else
    threshold = NaN;
end

% allocate memory to store each map
maps_cell = cell(stim_num,1);
for color = 1:stim_num
    maps_cell{color} = zeros(full_ref_dim);
end

% initialize a counter
counter = 1;
% run through all the seeds
for seeds = seed_map

    % get the indexes
    index_vector = coord(indexes==seeds&registered_anatomy>0);

    % for each color
    for color = 1:stim_num
        % if there's a threshold, check it
        if ~isnan(threshold)
            if threshold(seeds,color) == 0
                continue
            end
        end
        % add the value at that position of the given property
        maps_cell{color}(index_vector) = ...
            maps_cell{color}(index_vector) + abs(property(counter,color));

    end
    % update the seed counter
    counter = counter + 1;
    
end