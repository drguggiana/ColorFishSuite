function [file_info,stim_num,rep_num,z_num,tar_files] = parser(tar_path)

%load the files in the path
tar_files = dir(tar_path);

%get rid of the dots and notes/elog files
tar_files = tar_files(3:end-2);
tar_files = {tar_files.name}';

%get the number of files
file_num = size(tar_files,1);

%allocate memory to store the info about each file
file_info = zeros(file_num,4);

%for all the files
for files = 1:file_num
    %get the current file path
    curr_file = tar_files{files};
    
    %find the dot identifiers
    dots = strfind(curr_file,'_');
    
    
    %extract the pres, rep and stim
    p_file = str2double(curr_file(1:dots(1)-1));
    rep_file = str2double(curr_file(dots(1)+1:dots(2)-1));
    stim_file = str2double(curr_file(dots(2)+1:dots(3)-1));
   
    z_file = str2double(curr_file(dots(3)+1:dots(4)-1));

    %store the info in a common cell
    file_info(files,:) = [p_file,rep_file,stim_file,z_file];
end

%get rid of the 0 indexes in file_info
file_info(:,1:3) = file_info(:,1:3) + 1;

% %sort the files by stimulus
% [~,stim_sort] = sort(file_info(:,3),'ascend');

%get the number of reps used in the study
rep_num = max(file_info(:,2));

%get the number of stimuli used
stim_num = max(file_info(:,3));
%get the number of z_sections in the data set
[z_uni,~,ic] = unique(file_info(:,4));
z_num = length(z_uni);
%also replace the actual z's with identifiers in numerical order
file_info(:,4) = ic;
%and sort the file_info matrix and the tar_files based on presentation
%order
[~,order_vec] = sort(file_info(:,1));
file_info = file_info(order_vec,:);
tar_files = tar_files(order_vec);
