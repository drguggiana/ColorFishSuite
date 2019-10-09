function [im_corr] = Stack_process_1(astack,time_flag)
%the idea of this code is to process a single z stack into a matrix with
%averaged seeds in time. The function has to align the time traces of each
%pixel, calculate their temporal correlation and then assemble seeds based
%on these correlations (averaging the pixels involved while keeping the
%locations of the seeds available).I might have to come up with a better
%indicator than the temporal correlation but I'll start here
%A ton of the functions came basically straight from Martin Haesemeyer

%If it's the first correlation stage (across time)
if time_flag == 1
    %get the dimensions of the input matrix (parameters of the stack
    im_height = size(astack,1);
    im_width = size(astack,2);
%     time_num = size(astack,3);
    stim_num = size(astack,4);


    %allocate memory to store the correlation matrices
    im_corr = zeros(im_height,im_width,stim_num);

    %for all the stimuli
    for stim = 1:stim_num
        %get the current time stack
        tstack = squeeze(astack(:,:,:,stim));
    %     %align the stack across time
    %     [seeds_out,xshift,yshift] = realign_1(tstack);

    %     %perform quality control on the alignment of the stack
    %     [corrs,slices] = CorrelationControl_1(tstack,20);

        %perform the correlation of the stack across time
    %     profile on
    %     im_corr = AvgNeighbhorCorrelations(tstack,1);

        im_corr(:,:,stim) = Flat_corr_1(tstack);
    %     size(a)
    %     profile off
    %     profile viewer
    end
else %or the second (across stimuli)
    
    %allocate memory to store the correlation matrices
    im_corr = Flat_corr_1(astack);
end


function [corrs,slices] = CorrelationControl_1(tstack, nFrames)
   
% Sub-divides stack into nFrames blocks and cross-correlates
% each summed block to first reporting the 0-shift correlation to id
% potential movement artefacts. All slices will be z-scored to prevent
% differences in (raw) correlation values based on intensity (bleaching etc)



[h,w,nSlices] = size(tstack);
if floor(nSlices/nFrames) < 2
    error('Need to identify at least two nFrames sized sub-stacks in the stack');
end
ix0_x = h-1;
ix0_y = w-1;%coordinates of 0-shift correlation
slices = zeros(h,w,floor(nSlices/nFrames));
corrs = zeros(floor(nSlices/nFrames)-1,1);
for i = 1:floor(nSlices/nFrames)
    slices(:,:,i) = sum(tstack(:,:,1+nFrames*(i-1):nFrames*i),3);
    if i > 1
        temp_corr = xcorr2(zslice_fun(slices(:,:,1)),zslice_fun(slices(:,:,i)));
        corrs(i-1) = temp_corr(ix0_x,ix0_y);
    end
end

function zslice = zslice_fun(tar_slice)
zslice = (tar_slice-mean(tar_slice(:)))/std(tar_slice(:));

function im_corr = AvgNeighbhorCorrelations(tstack,dist)

% Returns a 2D image which for each pixel in stack
% has the average correlation of that pixel's time-series
% with the timeseries of it's neighbors up to dist pixels
% away from the current pixel.
% Predicate is an optional function that takes an x/y coordinate
% pair as argument and returns whether to compute the correlation.

if dist<1
    error('Dist has to be at least 1')
end
im_corr = zeros(size(tstack,1),size(tstack,2));

%get the max number of possible pixel values
max_pix = max(size(tstack,1),size(tstack,2));
%create a string map of the possible pizels
im_map = cell(max_pix,1);
%for all of the values
for val = 1:max_pix
    im_map{val} = num2str(val);
end

corr_buff = containers.Map('KeyType','char','ValueType','double');%buffers computed correlations to avoid computing the same pairs multiple times!
for x = 1:size(tstack,1)
    for y = 1:size(tstack,2)
%         if (not predicate is None) and (not predicate(x,y))%pixel is excluded
%             continue
%         end
        c_sum = zeros(2,1);
        for dx = -dist:dist%dx=dist inclusive!
            for dy = -dist:dist                    
                if dx==0 && dy==0%original pixel
                    continue
                end
                if x+dx<1 || y+dy<1 || x+dx>=size(im_corr,1) || y+dy>=size(im_corr,2)%outside of image
                    continue
                end
                p_src = [im_map{[x,y]}];
                p_des = [im_map{[x+dx,y+dy]}];
                if isKey(corr_buff,[p_src,p_des])
                    c_sum(1,1) = c_sum(1,1) + corr_buff([p_src,p_des]);
                    c_sum(2,1) = c_sum(2,1) + 1;
                else
                    cval = prcorr2(tstack(x,y,:),tstack(x+dx,y+dy,:));
%                     cval = cval_temp(1,2);
                    corr_buff([p_des,p_src]) = cval;
                    c_sum(1,1) = c_sum(1,1) + cval;
                    c_sum(2,1) = c_sum(2,1) + 1;
                end
            end
        end
        if ~isempty(c_sum) && ~all(isnan(c_sum))
            im_corr(x,y) = c_sum(1,1)./c_sum(2,1);
        end
    end
end
im_corr(isnan(im_corr)) = 0;

function im_corr = Flat_corr_1(tstack)

translation = [-1,-1;-1,0;-1,1;0,-1;0,1;1,-1;1,0;1,1]; % 8 directions
% [im_height,im_width,time_num] = size(astack);

%allocate memory for the correlation directions
corr_store = zeros(size(tstack,1),size(tstack,2),8);

az = bsxfun(@minus, tstack, mean(tstack,3));
a2 = az .^ 2;
for i = 1:8
%     sstack = imtranslate(tstack,translation(i,:));

    sstack = zeros(size(tstack));
    sstack(2+translation(i,1):end-1+translation(i,1)...
        ,2+translation(i,2):end-1+translation(i,2),:) = ...
        tstack(2:end-1,2:end-1,:);
%     A = reshape(tstack,[im_width*im_height,time_num])';
%     B = reshape(sstack,[im_width*im_height,time_num])';

    % Remove means
    
    bz = bsxfun(@minus, sstack, mean(sstack,3));
    % Standard Pearson correlation coefficient formula
    
    b2 = bz .^ 2;
    ab = az .* bz;
    corr_store(:,:,i) = sum(ab,3) ./ sqrt(sum(a2,3) .* sum(b2,3));
    
    % below is the faster ~equivalent (super close) of diag(corr(A,B))
%     An=bsxfun(@minus,A,mean(tstack(:),1));
%     Bn=bsxfun(@minus,B,mean(sstack(:),1));
%     An=bsxfun(@times,An,1./sqrt(sum(An.^2,1)));
%     Bn=bsxfun(@times,Bn,1./sqrt(sum(Bn.^2,1)));
%     C=sum(An.*Bn,1);
    %     disp(['translation #' num2str(i)]);
%     temp(:,:,i) = reshape(C,[im_width,im_height]);
end
im_corr = mean(corr_store,3);

function im_corr = Flat_corr_2(tstack)

translation = [-1,-1;-1,0;-1,1;0,-1;0,1;1,-1;1,0;1,1]; % 8 directions
[im_height,im_width,time_num] = size(tstack);

%generate a map with 8-neighbor coordinates (in linear indices) for all the
%positions in the image


%allocate memory for the correlation directions
corr_store = zeros(size(tstack,1)*size(tstack,2),size(tstack,3));
A = reshape(tstack,[im_width*im_height,time_num])';
for i = 1:8
    sstack = imtranslate(tstack,translation(i,:));
    B = reshape(sstack,[im_width*im_height,time_num])';

%     % Remove means
%     az = bsxfun(@minus, tstack, mean(tstack,3));
%     bz = bsxfun(@minus, sstack, mean(sstack,3));
%     % Standard Pearson correlation coefficient formula
%     a2 = az .^ 2;
%     b2 = bz .^ 2;
%     ab = az .* bz;
%     corr_store(:,:,i) = sum(ab, 3) ./ sqrt(sum(a2, 3) .* sum(b2, 3));
    corr_store = corr2(A,B) + corr_store;

    % below is the faster ~equivalent (super close) of diag(corr(A,B))
%     An=bsxfun(@minus,A,mean(tstack(:),1));
%     Bn=bsxfun(@minus,B,mean(sstack(:),1));
%     An=bsxfun(@times,An,1./sqrt(sum(An.^2,1)));
%     Bn=bsxfun(@times,Bn,1./sqrt(sum(Bn.^2,1)));
%     C=sum(An.*Bn,1);
    %     disp(['translation #' num2str(i)]);
%     temp(:,:,i) = reshape(C,[im_width,im_height]);
end

im_corr = reshape(corr_store,[size(tstack,1),size(tstack,2)]);