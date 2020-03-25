function [distances] = distance_calculation(data, idx_group, varargin)
% Calculate the distance between the indexed group separately for fish in
% the data

% check for the randomization flag
if size(varargin,2) == 1
    randomize = varargin{1};
else
    randomize = 0;
end

%get the number of groups
group_num = length(unique(idx_group));
%the fish origin vector
xy_mat = data.xy_seed;
z_mat = data.z_seed;
fish_ori = data.fish_ori;
%and get the number of fish
fish_num = size(unique(fish_ori(:,1)),1);
%allocate memory to store the distributions per cluster
group_dist = cell(group_num,fish_num);


%for all the clusters
for group = 1:group_num
    % if randomization is zero
    if randomize == 0
        % get the logical vector for this group
        idx_logic = idx_group==group;
    else
        % randomize the indexes picked
        num_seeds = sum(idx_group==group);
        rand_seed = randperm(length(idx_group),num_seeds);
        idx_zeros = zeros(length(idx_group),1);
        idx_zeros(rand_seed) = 1;
        idx_logic = idx_zeros>0;
    end
    %for all the fish
    for fish = 1:fish_num
        
      
        %get the coordinates corresponding to this cluster and fish
        fish_coord = cat(1,xy_mat(idx_logic&fish_ori(:,1)==fish).centroid);

        fish_z = z_mat(idx_logic&fish_ori(:,1)==fish);
        % get the z present in this fish
        z_list = unique(fish_z);
        % allocate memory for the distances
        z_cell = cell(length(z_list),1);
        % for all the z
        for z = 1:length(z_list)
            %calculate the distance distribution
            z_cell{z} = pdist(fish_coord(fish_z==z_list(z)))';
        end
        %concatenate the distances for this fish
        group_dist{group,fish} = vertcat(z_cell{:});
    end
    
end

distances = group_dist;