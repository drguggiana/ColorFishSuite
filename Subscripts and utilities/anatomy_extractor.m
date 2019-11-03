function [anatomy_info,area_all] = anatomy_extractor(tar_path,c_cell)

%% Load the anatomy information
    
%     if isempty(strfind(name_cell{1},'p17b')) == 0
%         stim_num2 = 4;
%     end

if contains(tar_path,'syn')
    %define the region labels
    reg_label = {'N/A','AF4','AF5','AF6','AF7','AF8','AF9','AF10'};
    reg_map = [0 4 5 6 7 8 9 10];
else
    %define the region labels
    reg_label = {'N/A','L-TcN','R-TcN','L-TcP','R-TcP','L-Cb','R-Cb','L-Hb','R-Hb','L-Pt','R-Pt'};
    reg_map = [0:10];
end

%get the number of regions
reg_num = length(reg_label);

%calculate the number of time bins
%     temp_conc = load(name_cell{1},'conc_trace');
%     temp_conc = temp_conc.conc_trace;

%     time_num = size(temp_conc,2)/stim_num2;

%     clear('temp_conc')
%% Define the AFs in each file (ultra ghetto, please make better when possible)

%define the af_list depending on the experiment

%allocate memory for the list
af_list = cell(39,1);
%go file by file defining the AFs present
af_list{1} = cell(2,1);
af_list{1}{1} = '20161001_syngc6s_p17b_1_anato';
af_list{1}{2} = [10 9 8 7];

af_list{2} = cell(2,1);
af_list{2}{1} = '20161001_syngc6s_p17b_1b_anato';
af_list{2}{2} = [6 5 4];

af_list{3} = cell(2,1);
af_list{3}{1} = '20161001_syngc6s_p17b_2_anato';
af_list{3}{2} = [10 9 8 7];

af_list{4} = cell(2,1);
af_list{4}{1} = '20161001_syngc6s_p17b_2b_anato';
af_list{4}{2} = [10 9 8 7 6 5 4];

af_list{5} = cell(2,1);
af_list{5}{1} = '20161002_syngc6s_p17b_1_anato';
af_list{5}{2} = [10 9 8 7];

af_list{6} = cell(2,1);
af_list{6}{1} = '20161002_syngc6s_p17b_2_anato';
af_list{6}{2} = [10 9];

af_list{7} = cell(2,1);
af_list{7}{1} = '20161002_syngc6s_p17b_2b_anato';
af_list{7}{2} = [10 9 8 7 6 5 4];

af_list{8} = cell(2,1);
af_list{8}{1} = '20161003_syngc6s_p17b_1_anato';
af_list{8}{2} = [10 9 8 7];

af_list{9} = cell(2,1);
af_list{9}{1} = '20161003_syngc6s_p17b_1b_anato';
af_list{9}{2} = [6 5 4];

af_list{10} = cell(2,1);
af_list{10}{1} = '20161004_syngc6s_p17b_1_anato';
af_list{10}{2} = [10 9 8 7];
%allocate memory for the list (based on the number of fish)
%         af_list = cell(num_data,1);
%fill up the list
%for all the fish
%         for fish = 1:num_data
%             af_list{fish} = cell(2,1);
%             [~,base_name,~] = fileparts(name_cell{fish});
%             af_list{fish}{1} = strcat(base_name(1:end-9),'_anato');
%             af_list{fish}{2} = 1:10;
%         end

% initialize a file counter (have to initialize after the last manual entry
% from above)
file_counter = 11;

% get the files left in the directory
list_dir = dir('E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\Stage1_extraction');
% for all the files in the folder that are not syngc6s or p8
for files = 1:size(list_dir,1)
    % if the file name contains p8 or syngc6s, skip it
    if contains(list_dir(files).name,{'p8','syngc6s'})
        continue
    elseif contains(list_dir(files).name,'_anato')
    %     [~,base_name,~] = fileparts(tar_path);
    %     af_list{file_counter}{1} = strcat(base_name(1:end-9),'_anato');
        % introduce a cell in this position
        af_list{file_counter} = cell(2,1);
        af_list{file_counter}{1} = list_dir(files).name;
        af_list{file_counter}{2} = 1:10;
        % increment the counter
        file_counter = file_counter + 1;
    end
end    
    
%% Load the AF files

%the structure of the files is that per each fish there are several files
%corresponding to the different imaging sessions. Each one of these files
%contains an nxm cell where n is z and m is AF. Within this cell there are
%1x2 cells that contain a px2 vector with the coordinates of the
%corresponding ROI and a hxw logical matrix containing the ROI in mask form

%define the path where the files can be found
af_path = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\anatomy_files';

%allocate memory to store the matrices
% af_cell = cell(num_data,1);

%the fish names
% fish_name = cell(num_data,1);
%and to store the list cells
% af_sublist = cell(num_data,1);

% %for all the fish
% for fish = 1:num_data
    
%load the file name of the fish to find out which folders to select
[~,f_name,~] = fileparts(tar_path);
% f_name =f_name;
% f_name = strcat(f_name(1:end-7),'_ROI.mat');
%     %split the name to get the folder names (fish and pre/post)
%     f_parts = strsplit(f_name,'_');
%     %assemble the path to load files from
%     f_path = strcat(af_path,'\',f_parts{1},'\',f_parts{2});

%     f_path = af_path;
%get the files in the path
af_files = dir(af_path);
af_files = {af_files(3:end).name};
%     %get the number of files
%     file_num = size(af_files,2);
%     %allocate memory for the file contents
%     af_subcell = cell(file_num,1);
%     %for all the files
%     for files = 1:file_num
%         %load the file
%         af_subcell{files} = load(fullfile(f_path,af_files{files}),'a_store');
%         af_subcell{files} = af_subcell{files}.a_store;
%     end

%     %also search for the list cell corresponding to this fish
%     fish_name{fish} = strcat(f_parts{1},'_',f_parts{2});
%     %extract just the names from the list
%     %get the number of items in the list
%     i_number = size(af_list,1);
%     %allocate memory for just the names
%     i_names = cell(i_number,1);
%     %for all the fish
%     for list_i = 1:i_number
%         i_names{list_i} = af_list{list_i}{1};
%     end
%     list_c = strcmp(fish_name{fish},i_names);
%
%     %store the af info in the cell
%     af_sublist{fish} = af_list{list_c}{2};


%     %concatenate and add to the main cell
%     af_cell{fish} = af_subcell;
% get the name of the ROI file
roi_file = contains(af_files,strcat(f_name,'_ROI.mat'));

% if the file exists, extract the info, otherwise exit the function
if sum(roi_file) > 0
    %load the file of interest
    f_inter = load(fullfile(af_path,af_files{roi_file}),'a_store');
    af_cell = f_inter.a_store;
    % af_cell{fish} = f_inter;
    % end
    %% Combine the AF and seed coordinate info to yield a matrix with AF and z

    % %allocate memory for the matrices per fish
    % af_fish = cell(num_data,1);

    %get a cell with just the af_list names
    af_names = horzcat(af_list{:});
    af_names = af_names(1,:)';
    % %for all of the fish
    % for fish = 1:num_data

    %assemble a volume with the labeled regions for each fish

    %     %get the number of files
    %     file_num = size(af_cell{fish},1);
    %     %allocate memory to store the temporary matrix
    %     af_filescell = cell(file_num,1);
    %     %allocate memory to store the volumes per file
    %     af_vfiles = cell(file_num,1);
    %     %for each file
    %     for files = 1:file_num
    %determine the number of z sections
    z_num = size(af_cell,1);
    %also determine the number of AFs
    af_num = size(af_cell,2);
    %for all the z sections
    for z = 1:z_num
        %if the cell is empty, skip to the next iterations
        if isempty(af_cell{z})
            continue
        else
            %determine the dimensions of the image
            im_height = size(af_cell{z}{2},1);
            im_width = size(af_cell{z}{2},2);
            break
        end
    end
    %allocate memory for the volume
    af_vol = zeros(im_height,im_width,z_num);
    %for every z section
    for z = 1:z_num
        %before performing the assignment of the traces, correct for AF
        %overlap
        %allocate memory for a temp frame
        temp_frame = zeros(im_height,im_width);
        %for all the AFs
        for afc = 1:af_num
            %if there are no ROIs there
            if isempty(af_cell{z,afc})
                continue
            end
            %accumulate the frames, so I can find the places with
            %overlap
            temp_frame = temp_frame + af_cell{z,afc}{2};
        end

        %find the places with overlap
        overlap_vec = find(temp_frame>1);
        %if there is overlap
        if ~isempty(overlap_vec)
            %for all the AFs
            for afc = 1:af_num
                %if there are no ROIs there
                if isempty(af_cell{z,afc})
                    continue
                end
                %NEED TO ASSIGN TO ONE OF THEM RANDOMLY!!!
                %eliminate those voxels from every af
                af_cell{z,afc}{2}(overlap_vec) = 0;
            end
        end
        %now actually assign the voxels to the respective AFs
        %for all the AFs
        for afc = 1:af_num
            %if there are no ROIs there
            if isempty(af_cell{z,afc})
                continue
            end
            %load the volume
            af_vol(:,:,z) = af_vol(:,:,z) + afc.*af_cell{z,afc}{2};
        end
    end
    %         %store the volume
    %         af_vfiles{files} = af_vol;
    %get the number of seeds in this file
    seed_num = size(c_cell,1);

    %allocate memory for the output matrix
    af_filemat = zeros(seed_num,2);
    %load the z
    af_filemat(:,2) = c_cell(:,3);
    %and index the corresponding volume to get the AF
    %for all the seeds
    for seeds = 1:seed_num
        %load the points
        p = c_cell(seeds,:);
        %and load the AF for the point
        af_filemat(seeds,1) = af_vol(p(2),p(1),p(3));
        %if the AF is zero, scan the surroundings of the points
        if af_filemat(seeds,1) == 0
            %allocate a vector to contain the surroundings
            val_vec = zeros(9,1);
            vec_c = 1;
            for x = -1:1
                for y = -1:1
                    val_vec(vec_c) = af_vol(p(2)+x,p(1)+y,p(3));
                    vec_c = vec_c + 1;
                end
            end
            %define the target AF as the most common for the seed
            af_filemat(seeds,1) = mode(val_vec(val_vec~=0));
        end

    end
    %         %store the matrix
    %         af_filescell{files} = af_filemat;
    %     end

    %concatenate all the files and store in the final output cell
    %     af_fish{fish} = cat(1,af_filescell{:});

    af_fish = af_filemat;

    % %get the name of the fish
    % [~,fish_name,~] = fileparts(name_cell{fish});
    % %make the names compatible with the list
    % fish_name = strcat(fish_name(1:end-7),'_anato');
    % %check whether the file is in the list
    % list_vec = strcmp(fish_name,af_names)==1;

    list_vec = contains(af_names,strcat(f_name,'_anato'));
    %if it's a set of files with AFs
    if any(list_vec) == 1
        %find the AFs in the file
        af_shortvec = af_list{list_vec}{2};
        %for all of the AFs in this file
        for afs = 1:length(af_shortvec)
            %go through the AFs changing them to their actual label
            af_fish(af_fish(:,1)==afs) = af_shortvec(afs);
        end
    end

    % get the info out
    anatomy_info = af_fish;

    %% Calculate the area of each region

    % %allocate memory for the areas
    % area_cell = cell(num_data,1);

    %get the number of z sections
    z_num = size(af_cell,1);
    %get the number of regions
    reg_here = size(af_cell,2);
    %allocate memory for the area per region
    area_reg = zeros(reg_num,1);
    %for each z section
    for z = 1:z_num
        %for each region
        for regs = 1:reg_here
            %if it's empty, skip it
            if isempty(af_cell{z,regs})
                continue
            end
            %get the absolute number of the region (and its index)
            abs_reg = reg_map==af_list{list_vec}{2}(regs);
            %add the area of the region in this plane
            area_reg(abs_reg) = area_reg(abs_reg) + sum(sum(af_cell{z,regs}{2}));
        end
    end

    %calculate the total area per region
    area_all = area_reg;
else
    anatomy_info = [];
    area_all = [];
end
% end