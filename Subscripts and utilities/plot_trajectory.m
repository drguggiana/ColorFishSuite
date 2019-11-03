function plot_trajectory(pca_mat,plot_col,options)

%for all the points
for points = 1:size(pca_mat,1)
    if options.threeD == 1
        plot3(pca_mat(points,1),pca_mat(points,2),pca_mat(points,3),...
            'Marker','o','MarkerSize',points/5,'MarkerFaceColor',plot_col,...
            'MarkerEdgeColor',plot_col)
        hold('on')

    else
        subplot(1,2,1)
        plot(pca_mat(points,2),pca_mat(points,1),...
            'Marker','o','MarkerSize',points/5,'MarkerFaceColor',plot_col,...
            'MarkerEdgeColor',plot_col)
        hold('on')
        subplot(1,2,2)
        plot(pca_mat(points,3),pca_mat(points,1),...
            'Marker','o','MarkerSize',points/5,'MarkerFaceColor',plot_col,...
            'MarkerEdgeColor',plot_col)
        hold('on')
    end
    
end

if options.line == 1
    if options.threeD == 1
        %also plot the lines
        plot3(pca_mat(:,1),pca_mat(:,2),pca_mat(:,3),...
            'Color',plot_col)
    else
        subplot(1,2,1)
        plot(pca_mat(:,2),pca_mat(:,1),...
            'Color',plot_col)
        subplot(1,2,2)
        plot(pca_mat(:,3),pca_mat(:,1),...
            'Color',plot_col)
    end
end