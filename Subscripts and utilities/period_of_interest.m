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
    case 4 % beginning of stim period
        rest_all = logical(repmat([rest_vec;stim_vec;rest_vec;rest_vec],stim_num*rep_num,1));
    case 5 % end of stim period
        rest_all = logical(repmat([rest_vec;rest_vec;stim_vec;rest_vec],stim_num*rep_num,1));
    case 6 % first 10 time frames
        rest_all = logical(repmat([rest_vec;stim_vec(1:10);rest_vec(1:10);rest_vec;rest_vec],stim_num*rep_num,1));
    case 7 % last 10 time frames
        rest_all = logical(repmat([rest_vec;rest_vec;rest_vec(1:10);stim_vec(11:end);rest_vec],stim_num*rep_num,1));
    case 8 %9 frames of pre stim and 9 frames in the middle of stim
        rest_all = logical(repmat([rest_vec(1:11);stim_vec(1:9);rest_vec;stim_vec(1:9);...
            rest_vec(1:11);rest_vec],stim_num*rep_num,1));
    case 9 %9 frames of pre stim and 9 frames at the start of stim
        rest_all = logical(repmat([rest_vec(1:11);stim_vec(1:9);stim_vec(1:9);rest_vec(1:11);...
            rest_vec;rest_vec],stim_num*rep_num,1));
end