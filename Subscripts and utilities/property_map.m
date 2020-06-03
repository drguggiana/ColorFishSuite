function maps_cell = property_map(coord,indexes,registered_anatomy,property,stim_num,full_ref_dim,varargin)

% % get the reshaped activity
% property = reshape(data.conc_trace,[],data.time_num,data.stim_num);
% % take only the stimulation time
% property = property(:,21:60,:);
% % take the absolute average
% property = squeeze(mean(abs(property),2));

if length(varargin) >= 1
    seed_map = varargin{1};
else
    seed_map = 1:size(property,1);
end

% % get the number of stimuli
% stim_num = data.stim_num;
% allocate memory to store each map
maps_cell = cell(stim_num,1);
for color = 1:stim_num
    maps_cell{color} = zeros(full_ref_dim);
end

% initialize a counter
counter = 1;
% run through all the seeds
for seeds = seed_map
    index_vector = coord(indexes==seeds&registered_anatomy>0);
    %         % accumulate the gain for the pixels of each seed
    %         gain_maps{max_idx(seeds)}(index_vector) = ...
    %             gain_maps{max_idx(seeds)}(index_vector) + max_gain(seeds);
    % for each color
    for color = 1:stim_num
        maps_cell{color}(index_vector) = ...
            maps_cell{color}(index_vector) + abs(property(counter,color));
        %             prc = prctile(gain_maps{color}(:),90);
        %             gain_maps{color}(gain_maps{color}>prc) = prc;
    end
    % update the seed counter
    counter = counter + 1;
    
end