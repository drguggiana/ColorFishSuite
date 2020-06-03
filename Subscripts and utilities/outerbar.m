function h = outerbar(file_path,ax_handle,bar_label)

% create the figure and colorbar
figure
% set(gca,'Visible','off')
h = colorbar;
% axis off

% set the bar properties
set(h,'TickLength',0,'LineWidth',2)
ylabel(h,bar_label)

% get the colormap and ranges from the figure
cmap = get(ax_handle,'colormap');

colormap(h,cmap)
% add the name
bar_path = strrep(file_path,'.eps','bar.eps');
% save the bar as a separate image
% print(fullfile(fig_path,bar_path),'-dpng','-r600')
export_fig(bar_path,h)

