function plot_trajectory(pca_mat,plot_col,options)

% if the marker option is not there, set o as default
if ~isfield(options,'marker')
    options.marker = 'o';
end

if options.line == 1
    if options.threeD == 1
        %also plot the lines
        plot3(pca_mat(:,1),pca_mat(:,2),pca_mat(:,3),...
            'Color',plot_col,'LineWidth',1)
        hold on
    else
        subplot(1,2,1)
        plot(pca_mat(:,2),pca_mat(:,1),...
            'Color',plot_col,'LineWidth',1)
        hold on
        subplot(1,2,2)
        plot(pca_mat(:,3),pca_mat(:,1),...
            'Color',plot_col,'LineWidth',1)
        hold on
    end
end

% get the point size vector
psize = (1:size(pca_mat,1))./3;
%for all the points
% for points = 1:size(pca_mat,1)
    if options.threeD == 1
        scatter3(pca_mat(:,1),pca_mat(:,2),pca_mat(:,3),psize,...
            'Marker',options.marker,'MarkerFaceColor',plot_col,...
            'MarkerEdgeColor',[0 0 0])
        grid off

    else
        subplot(1,2,1)
        scatter(pca_mat(:,2),pca_mat(:,1),psize,...
            'Marker',options.marker,'MarkerFaceColor',plot_col,...
            'MarkerEdgeColor',[0 0 0])
        hold('on')
        subplot(1,2,2)
        scatter(pca_mat(:,3),pca_mat(:,1),psize,...
            'Marker',options.marker,'MarkerFaceColor',plot_col,...
            'MarkerEdgeColor',[0 0 0])
        hold('on')
    end
    
% end

