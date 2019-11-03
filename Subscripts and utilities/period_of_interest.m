function rest_all = period_of_interest(period, stim_num, rep_num)

%get a vector to signal the rest periods
rest_vec = zeros(20,1);
stim_vec = ones(20,1);
% select the period
switch period
    case 0 % pre_stim
        rest_all = logical(repmat([stim_vec;rest_vec;rest_vec;rest_vec],stim_num*rep_num,1));
    case 1 % stim
        rest_all = logical(repmat([rest_vec;stim_vec;stim_vec;rest_vec],stim_num*rep_num,1));
    case 2 % post_stim
        rest_all = logical(repmat([rest_vec;rest_vec;rest_vec;stim_vec],stim_num*rep_num,1));
    case 3 % pre and post_stim
        rest_all = logical(repmat([stim_vec;rest_vec;rest_vec;stim_vec],stim_num*rep_num,1));
end