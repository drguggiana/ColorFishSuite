function h = projection_plot(input_mat,varargin)

h = figure;
%plot the projections
for dims = 1:3
    %average across the desired dimension
%     ave_z = squeeze(sum(input_mat,dims));
    ave_z = squeeze(max(input_mat,[],dims));
    %if there is a ref stack provided
    if nargin > 1
        ref_stack = varargin{1};
        ref_stack = squeeze(max(ref_stack,[],dims));
        
        ave_z = imfuse(ave_z,ref_stack,'falsecolor');
    end
    %select the subplot (projection)
    switch dims
        case 1
            subplot(2,5,1:3)
            imagesc(permute(ave_z,[2 1 3]))
            set(gca,'YDir','normal')
        case 2
            subplot(2,5,6:8)
            imagesc(permute(ave_z,[2 1 3]))
            set(gca,'YDir','normal')
        case 3
            subplot(2,5,[4:5 9:10])
            imagesc(ave_z)
    end
%     axis equal
end