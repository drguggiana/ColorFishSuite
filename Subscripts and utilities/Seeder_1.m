function [seed_list,seed_im] = Seeder_1(corr_pic,thres_seed,thres_nb,thres_area,thres_minarea)

%get the dimensions of the image
im_height = size(corr_pic,1);
im_width = size(corr_pic,2);
          
           
ROI = [];
ROI(300).pxlist = []; % 300 not capping
ROI(300).area = [];
ROI(300).centroid = [];

BW = zeros(im_width,im_height,'uint8'); % mask to store new ROI
ROIcount = 0;

% find 1st seed
m = corr_pic;
[~,IX] = max(m(:));
[I,J] = ind2sub([im_width,im_height],IX);

while m(I,J)>thres_seed, % each seed
    if (I==1 || I==im_width || J==1 || J==im_height),% edge of image
        % find next seed and skip to next loop
        m(I,J) = NaN;
        [~,IX] = max(m(:));
        [I,J] = ind2sub([im_width,im_height],IX);
        continue;
    end
    % add seed position to mask
    BW(I,J) = 1;
    
    % grow ROI: expand mask and examine next neighbors
    IX2 = 0; 
    area = 0; % initialize to dummy
    while ~isempty(IX2) && area<=(thres_area-4), % 4 to grow this round
%         %get x and y coordinates of the seed
%         [x,y] = find(BW);
%         BW = BW(min(x)-1:max(x)+1,min(y)-1:max(y)+1);
        
        BWp = bwpack(BW);
        BW0 = imdilate(BWp,[0 1 0; 1 1 1; 0 1 0],'ispacked');% 4 neighbors
        BW0 = bwunpack(BW0,im_height);
        
        BW0(logical(BW))=0;
        IX1 = find(BW0);
        IX2 = find(m(IX1)>thres_nb);
        if ~isempty(IX2),
            BW(IX1(IX2)) = 1;
            area = length(find(BW));
        end
    end % finished expanding this seed
    
    % save this ROI if bigger than thres
    area = length(find(BW));
    if area >= thres_minarea,% = 3 for ds 0.2
        ROIcount = ROIcount+1;
        ROI(ROIcount).pxlist = find(BW);
        center = regionprops(BW,'centroid');
        ROI(ROIcount).centroid = center.Centroid;
        %                          ROI(ROIcount).pxlist = find(BW)';
        ROI(ROIcount).area = area;
    end
    % reset
    m(logical(BW)) = NaN;
    BW = zeros(im_width,im_height);
    % look for next seed
    [~,IX] = max(m(:));
    [I,J] = ind2sub([im_width,im_height],IX);
end
% no seeds above thres_nb found.
ROI(ROIcount+1:end) = [];

%image with the seeds found
cmap = rand(ROIcount,3);
clrim = zeros(im_width,im_height,3);
for i = 1:ROIcount,
    im = zeros(im_width,im_height);
    im(ROI(i).pxlist) = 1;
    for k = 1:3,
        clrim(:,:,k) = clrim(:,:,k)+im * cmap(i,k);
    end
end
%send outputs
seed_im = clrim;
seed_list = ROI;