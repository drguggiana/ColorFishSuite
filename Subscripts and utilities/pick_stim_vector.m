function stim_vector = pick_stim_vector(stim_set,stim_num,stim_name)
%get the stimulus vector to use based on variable and data set selection
% get the stimulus set
switch stim_set
    case 1
        % standard full stim set
        stim_vector = 1:stim_num;
    case 2
        % alternate red uv only
        if contains(stim_name,'p17b')
            stim_vector = [1,4];
        else
            error('Must use p17b as a dataset')
        end
    case 3
        % alternate red uv CK only
        if contains(stim_name,'p17b')
            error('Must use p8 as a dataset')
        else
            stim_vector = [1,2];
        end
    case 4
        % alternate red uv GR only
        if contains(stim_name,'p17b')
            error('Must use p8 as a dataset')
        else
            stim_vector = [3,4];
        end
    case 5
        % alternate red uv FL only
        if contains(stim_name,'p17b')
            error('Must use p8 as a dataset')
        else
            stim_vector = [5,6];
        end
end