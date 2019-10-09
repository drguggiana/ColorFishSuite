function [singlez_mat,shift_mat,ave_frame,err_mat] = ...
    aligner_7(tar_path,tar_files,stim_num,rep_num,...
    file_info,pre_time,tar_z,skip_stim,align_factor)

%using the first file, find out the number of time points in the trace
first_info = imfinfo(fullfile(tar_path,tar_files{1}));
time_num = length(first_info);
%also get the image dimensions
im_width = first_info(1).Width;
im_height = first_info(1).Height;

%allocate memory for an average frame
ave_frame = zeros(im_width,im_height);

%define the target z (replace with for loop once code runs)
z = tar_z;
%get the number of excluded stimuli
exc_num = length(skip_stim);

%allocate memory for the averaged and aligned stacks for this z section
singlez_mat = zeros(im_height,im_width,time_num*(stim_num+exc_num*rep_num),'single');
%allocate memory for the output SNR matrix
err_mat = zeros(im_height,im_width,stim_num);
%also initialize a matrix to store the shifts in x and y
shift_mat = zeros(rep_num,time_num,stim_num,2);
%initialize a stack counter
st_counter = 1;
%and an exclusion counter
ex_counter = time_num*stim_num+1;
%for all of the stimuli
for stim = 1:stim_num
    
    fprintf(strcat('Stim: ',num2str(stim),'\n'))
    %allocate memory to store the stacks
    singlestim_mat  = zeros(im_height,im_width,rep_num,time_num,'uint16');
    %find the file coordinates involved with this stimulus
    curr_coord = find(file_info(:,3)==stim&file_info(:,4)==z);
    %for all the reps
    for reps = 1:rep_num
       %for all the time points
        for timep = 1:time_num     
            %load each frame
            singlestim_mat(:,:,reps,timep) = imread(fullfile(tar_path,tar_files{curr_coord(reps)}),timep);
        end
    end
    %create a temporary copy of singlestim in double precision
    singlestim_copy = double(singlestim_mat);
    %calculate the error
    err_mat(:,:,stim) = var(squeeze(mean(singlestim_copy,3)),0,3)./mean(squeeze(var(singlestim_copy,0,4)),3);
    %delete the copy
    clear('singlestim_copy');
    
    %calculate the average frame for anatomy
    ave_frame = ave_frame + mean(mean(singlestim_mat,4),3)/stim_num;
    %create an averaged stack for alignment only
    %get the resulting number of slices
    ave_slices = time_num/align_factor;
    %allocate memory for the stack
    ave_stack = zeros(im_height,im_width,rep_num,ave_slices);
    %for all the new slices
    for slices = 1:ave_slices
        %calculate the stack
        ave_stack(:,:,:,slices) = ...
            mode(singlestim_mat(:,:,:,1+align_factor*(slices-1):align_factor*slices),4);
    end
    
    %define the max allowed shift
    max_shift = 3;
    %align the stack before averaging
    [astack,shift_mat(:,:,stim,1),shift_mat(:,:,stim,2)] = ...
        realign_1(singlestim_mat,1,max_shift,uint16(ave_stack));
 
    %turn astack into a single and subtract baseline from microscope
    %WARNING, SYSTEM SPECIFIC
    astack = single(astack) - 99;
    
    %use bsxfun to calculate dfof
    pre_ave = mean(astack(:,:,:,pre_time),4);
%     pre_ave = prctile(astack,50,4);
%     pre_all = astack(:,:,:,pre_time);
%     figure
%     histogram(pre_all(:),1000)
    dfof_mat = bsxfun(@minus,astack,pre_ave);
    dfof_mat = bsxfun(@rdivide,dfof_mat,pre_ave);
    dfof_mat(isnan(dfof_mat)) = 0;
    
    %if it's a skip (RF stim)
    if any(skip_stim==stim)
        %take the mode of the reps
        singlez_mat(:,:,st_counter:stim*time_num) = squeeze(mode(dfof_mat,3));

        %and also position the concatenated reps for response calculation
        %at the end of the stack
        %for all the reps
        for reps = 1:rep_num
            singlez_mat(:,:,ex_counter:ex_counter+time_num-1) = squeeze(dfof_mat(:,:,reps,:));
            %update the exclusion counter
            ex_counter = ex_counter + time_num;
        end
        %if not
    else
        %average the reps together
        singlez_mat(:,:,st_counter:stim*time_num) = squeeze(mean(dfof_mat,3));
    end
    %update the stim counter
    st_counter = st_counter + time_num;
    
end

function [astack,x_shifts,y_shifts] = realign_1(tstack,shape_flag,max_shift,ave_stack)

%if the shape flag is active
if shape_flag == 1
    %get the original number of time slices and reps
    ori_rep = size(tstack,3);
    ori_time = size(tstack,4);
    %and of the average stack
    ori_ave_time = size(ave_stack,4);
    
    %reshape the matrix into 3 dimensions only
    tstack = reshape(tstack,[size(tstack,1),size(tstack,2),ori_rep*ori_time]);
    %and reshape the ave matrix too
    ave_stack = reshape(ave_stack,...
        [size(ave_stack,1),size(ave_stack,2),ori_rep*ori_ave_time]);
end

t_num = size(tstack,3);
x_shifts = zeros(t_num,1);
y_shifts = zeros(t_num,1);
astack = tstack;
sstack = uint16(sum(ave_stack,3));
%get the number of time points after averaging
t_avenum = size(ave_stack,3);
%calculate the number of slices being averaged
ave_slices = t_num/t_avenum;
for t = 1:t_avenum
    [xshift,yshift] = ComputeAlignmentShift(ave_stack,t,sstack);

    
    if xshift == 0 && yshift == 0
        continue
    end
    if abs(xshift)>max_shift || abs(yshift)>max_shift
        fprintf(strcat('Warning. Slice_',num2str(t)...
            ,'_requires shift greater than_',num2str(max_shift),'_pixels. Not shifted\n'))
        continue
    end
    [xs,xt] = Shift2Index(xshift,size(astack,1));
    [ys,yt] = Shift2Index(yshift,size(astack,2));
    
    %get the range for this shift (the slices encompassed in the average)
    t_range = 1+(t-1)*ave_slices:t*ave_slices;
    newImage = zeros(size(astack,1),size(astack,2),ave_slices);
    newImage(xt(1):xt(2),yt(1):yt(2),:) = ...
        astack(xs(1):xs(2),ys(1):ys(2),t_range);
    %also update the vector with the shifts for each frame
    x_shifts(t_range) = xshift;
    y_shifts(t_range) = yshift;


    astack(:,:,t_range) = newImage;

end
%if the shape flag is on
if shape_flag == 1
    %reshape the matrix into 3 dimensions only
    astack = reshape(astack,[size(astack,1),size(astack,2),ori_rep,ori_time]);
    x_shifts = reshape(x_shifts,[ori_rep,ori_time]);
    y_shifts = reshape(y_shifts,[ori_rep,ori_time]);
end

function [xshift,yshift] = ComputeAlignmentShift(astack,index,sstack)
   
%For the slice in stack identified by index computes
%the x (row) and y (column) shift that corresponds to
%the best alignment of stack[index,:,:] to the re-
%mainder of the stack
    
%mod to speed up the code
sum_stack = sstack-astack(:,:,index);

exp_x = size(astack,1)/2 + 1;
exp_y = size(astack,2)/2 + 1;%these are the indices in the cross-correlation matrix that correspond to 0 shift

A = real(fft2(sum_stack));
B = real(fft2(astack(:,:,index)));

c = ifftshift(ifft2(A.*conj(B)));

% imagesc(c)
% assignin('base','sum_stack',sum_stack)
% assignin('base','astack',astack(:,:,1))
[~,xy_temp] = max(c(:));
[x,y] = ind2sub(size(c),xy_temp);
xshift = x-exp_x;
yshift = y-exp_y; 
% xshift = x-1;
% yshift = y-1;

function [source,target] = Shift2Index(s_shift,s_size)
       
%Translates a given shift into the appropriate
%source and target indices

if s_shift<0
    %coordinate n in source should be n-shift in target
    source = [1-s_shift,s_size];
    target = [1,s_size+s_shift];
elseif s_shift>0
    %coordinate n in source should be n+1 in target
    source = [1,s_size-s_shift];
    target = [1+s_shift,s_size];
else
    source = [1,s_size];
    target = [1,s_size];
end