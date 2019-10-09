function out_mat = normr_1(in_mat,row_flag)
%normalize the input matrix row by row (subtract the min and then divide by
%the new max)

if isempty(row_flag)
    row_flag = 0;
end
switch row_flag
    %normalize across rows
    case 0
        temp_var = bsxfun(@minus,in_mat,min(in_mat,[],2));
        out_mat = bsxfun(@rdivide,temp_var,max(temp_var,[],2));
        %make NaNs 0
        out_mat(isnan(out_mat)) = 0;
    %normalize across the entire matrix
    case 1
        out_mat = (in_mat-min(in_mat(:)))./(max(in_mat(:))-min(in_mat(:)));
    %normalize across columns
    case 2
        temp_var = bsxfun(@minus,in_mat,min(in_mat,[],1));
        out_mat = bsxfun(@rdivide,temp_var,max(temp_var,[],1));
        %make NaNs 0
        out_mat(isnan(out_mat)) = 0;
    %normalize across columns with sum of the column
    case 3
        temp_var = bsxfun(@minus,in_mat,min(in_mat,[],1));
        out_mat = bsxfun(@rdivide,temp_var,sum(temp_var,1));
        %make NaNs 0
        out_mat(isnan(out_mat)) = 0;
    %normalize across rows with sum of the row
    case 4
        temp_var = bsxfun(@minus,in_mat,min(in_mat,[],2));
        out_mat = bsxfun(@rdivide,temp_var,sum(temp_var,2));
        %make NaNs 0
        out_mat(isnan(out_mat)) = 0;
end