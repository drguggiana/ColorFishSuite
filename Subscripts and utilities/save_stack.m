function save_stack(path,stack,varargin)
% function to save a stack to tif file

% if it's a 4D stack (i.e. in color)
if length(size(stack)) == 4
    % get the number of z sections (assumed as 3rd dimension)
    z_num = size(stack,3);
    % create the file and write the first frame
    imwrite(squeeze(stack(:,:,1,:)),path,'tif','WriteMode','Overwrite');
    if z_num > 1
        % write the rest
        for z = 2:z_num
            imwrite(squeeze(stack(:,:,z,:)),path,'tif','WriteMode','append')
        end
    end
    % set the alpha flag
    if size(stack,4) == 4
        alpha_channel = true;
    else
        alpha_channel = false;
    end
else
    % get the number of z sections (assumed as 3rd dimension)
    z_num = size(stack,3);
    % create the file and write the first frame
    imwrite(stack(:,:,1),path,'tif','WriteMode','Overwrite')
    if z_num > 1
        % write the rest
        for z = 2:z_num
            imwrite(stack(:,:,z),path,'tif','WriteMode','append')
        end
    end
    % set the alpha channel flag
    alpha_channel = false;
end

% if ref_info was provided
if length(varargin) > 1
    ref_info = varargin{1};
    % create the image object
    im_object = Tiff(path,'r+');
    % add the resolution
    setTag(im_object,'ImageDescription',ref_info(1).ImageDescription);
    setTag(im_object,'XResolution',ref_info(1).XResolution);
    setTag(im_object,'YResolution',ref_info(1).YResolution);
    % if the alpha was specified, add it
    if alpha_channel
        setTag(im_object,'ExtraSamples',Tiff.ExtraSamples.Unspecified);
    end

    rewriteDirectory(im_object);
    close(im_object);
end