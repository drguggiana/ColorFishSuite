function fishave_cell = reformat_fish(data,fish_ori,num_files,im_info_cell,aff_cell,ref_info,target)

% get the seed coordinates
cat_seed_all = data.xy_seed;
cat_z_all = data.z_seed;

% define the frame size based on the dataset
if contains(data.name,'p8_SynG6s')
    frame_size = 256;
else
    frame_size = 320;
end

%allocate memory for the registered fish
fishave_cell = cell(num_files,2);

%get the number of fish in the coordinate file
num_fish_coord = length(unique(fish_ori(:,1)));
 
%allocate memory for the final structure
seeds_px = vertcat(cat_seed_all(:).pxlist);
fish_xyz = zeros(size(seeds_px,1),3);

%initialize a seed counter
seed_c = 1;
%and also allocate memory for a vector to replace fish_ori
new_ori = zeros(size(seeds_px,1),1);
%and the stimulus vector
new_stim = zeros(size(seeds_px,1),1);
%for all the seeds
for seeds = 1:size(cat_seed_all,1)
    switch target
        case 'centroid'
            % use the centroids
            num_pix = 1;
            list_pix = cat_seed_all(seeds).centroid;
            sub_x = list_pix(1);
            sub_y = list_pix(2);
        case 'seeds'
            %get the number of pixels in this seeds
            num_pix = cat_seed_all(seeds).area;
            %and the actual pixels
            list_pix = cat_seed_all(seeds).pxlist;
            %transform the indexes to subindexes
            [sub_y,sub_x] = ind2sub([frame_size frame_size],list_pix);
    end
    %flip the x indexes
    sub_y = frame_size - sub_y + 1;
    %load the transformed indexes and the z into the fish_xyz
    %matrix
    fish_xyz(seed_c:seed_c+num_pix-1,:) = ...
        [sub_x,sub_y,ones(num_pix,1).*cat_z_all(seeds)];
    %load the fish of origin and stimulus into the new vector
    new_ori(seed_c:seed_c+num_pix-1) = fish_ori(seeds,1);
    %do the same with the seed ID
    new_stim(seed_c:seed_c+num_pix-1) = repmat(seeds,num_pix,1);
    
    %update the index
    seed_c = seed_c + num_pix;
    
end
% transform dimensions to match the orientation of the registered data
% allocate memory for the data for each stack
stack_cell = cell(num_fish_coord,1);
% get the z size of every stack for later subtraction
z_sizes = [0;cumsum(cellfun(@(x) size(x,1),im_info_cell))];
% for each stack
for stack = 1:num_fish_coord
    % get the im_info
    im_info = im_info_cell{stack};
    % get the traces for this stack
    stack_traces = fish_xyz(new_ori==stack,:);
    % get the z dimensions of the original stack
    original_z = size(im_info,1);
    % create a mapping between then non-inverted and the inverted z
    z_map = original_z:-1:1;
    % correct the trace z to within the stack
    stack_traces(:,3) = stack_traces(:,3)-z_sizes(stack);
    % replace the z in the traces
    stack_traces(:,3) = z_map(stack_traces(:,3));
%     %invert z
%     fish_xyz(:,3) = max(fish_xyz(:,3)) - fish_xyz(:,3)+1;
% 
%     %renumber z
%     fish_xyz(:,3) = mod(fish_xyz(:,3)-1,max(cat_z_all)/num_fish_coord)+1;

    %invert x and y
    stack_traces(:,[2 1]) = stack_traces(:,[1 2]);

    %rescale dimensions based on the micron/pixel divergence

    % apply it to the stack and correct for the removed border (always 5
    % pixels)
    stack_traces(:,[1 2]) = (stack_traces(:,[1 2])-5)/im_info(1).XResolution;
    % get the slice spacing
    spacing_string = im_info(1).ImageDescription;
    spacing_coord = strfind(spacing_string,'spacing=');
    dot_coord = strfind(spacing_string(spacing_coord:end),'.')+spacing_coord;
    z_spacing = str2double(spacing_string(spacing_coord+8:dot_coord(1)));
    % convert with the slice spacing
    stack_traces(:,3) = stack_traces(:,3)*z_spacing;
    % store the traces in the cell
    stack_cell{stack} = stack_traces;
end

% concatenate the stacks
fish_xyz = cat(1,stack_cell{:});

%initialize a counter for the fish
fish_c = 1;
%for all the files
for files = 1:num_files
    %get the coordinates for the current fish
    curr_fish = fish_xyz(new_ori==files,:);

    %transform the volume
    new_fishave_xyz = transformPointsInverse(aff_cell{fish_c},curr_fish);
    
    %bring the xy coordinates back to pixels
    new_fishave_xyz(:,[2 1]) = (new_fishave_xyz(:,[1 2]))*ref_info(1).XResolution;
    % get the z spacing
    spacing_string = ref_info(1).ImageDescription;
    spacing_coord = strfind(spacing_string,'spacing=');
    dot_coord = strfind(spacing_string(spacing_coord:end),'.')+spacing_coord;
    ref_z_spacing = str2double(spacing_string(spacing_coord+8:dot_coord(1)));
    % rescale also z
    new_fishave_xyz(:,3) = new_fishave_xyz(:,3)/ref_z_spacing;

    % apply the ref conversion and store the data in the fish cell
    fishave_cell{fish_c,1} = new_fishave_xyz; 
    fishave_cell{fish_c,2} = new_stim(new_ori==files,:);
    
    %update the counter
    fish_c = fish_c + 1;
    
end