function handle_out = style_figure(fig_handle,fig_set)
% Style the provided figure based on the given settings

%% Set the defaults

% define the default values
default_settings = struct([]);
% default colormap
default_settings(1).cmap = magma;
% default font size. the idea is not to have too many font sizes in a
% paper, so I went with 3. you define them down there and this just deines
% which one to use with the figure as default. Check journal requirements,
% though usually it's 7pt
default_settings(1).font_size = 'small';
% default plot X label
default_settings(1).XLabel = '';
% defaults plot y label
default_settings(1).YLabel = '';
% default title, though better not use cause export_fig can't handle it
default_settings(1).Title = '';
% default colorbar label
default_settings(1).colorbar_label = '';
% default for the squared plot box (up and right edges)
default_settings(1).box = 'off';
% whether to include a colorbar, default is off
default_settings(1).colorbar = 0;
% default size (width, height) in cm
default_settings(1).fig_size = [10 10];
% whether to crop, most instances it is better to do so
default_settings(1).crop = 1;
% default font, Arial is a good go to that's also sans serif
default_settings(1).FontName = 'Arial';
% disregard, mostly an old feature
default_settings(1).font_factor = 1;
% default linewidth for the axes, don't go smaller
default_settings(1).LineWidth = 1;
% whether to use painter rendering, activate if the plot resolution is bad
default_settings(1).painters = 0;
% whether to skip an axis, for multiaxes plots/subplots
default_settings(1).skip = 0;

% get the axis handle
ax_list = get(fig_handle,'Children');

% if there are more objects, get the axes
if length(ax_list)>1
    ax_list = findobj(ax_list,'Type','Axes');
end

% get the number of axes
ax_number = length(ax_list);

% for all the axes
for ax_count = 1:ax_number
    
    % get the handle
    ax = ax_list(ax_count);
    % apply them to the structure if the field doesn't exist
    field_list = fields(default_settings);
    % for all the fields
    for field = field_list'
        if ~isfield(fig_set,field{1}) || isempty(fig_set(ax_count).(field{1}))
            % filter for the labeling and ticks fields that could've been
            % specified before
            if any(contains({'XLabel','YLabel','XTick','YTick','XTickLabels','YTickLabels','Title'},field{1})) && ...
                  ~isempty(ax.(field{1}))  
                fig_set(ax_count).(field{1}) = ax.(field{1}).String;
            else
                fig_set(ax_count).(field{1}) = default_settings.(field{1});
            end
        end
    end
    
    % if the skip flag is there, skip the axes
    if fig_set(ax_count).skip
        continue
    end

    % get the font factor
    font_factor = fig_set(ax_count).font_factor;

    % define the font sizes
    switch fig_set(ax_count).font_size
        case 'large'
            font_size_target = 12;
        case 'medium'
            font_size_target = 10;
        case 'small'
            font_size_target = 7;
    end


    % get rid of the inner tick marks
    set(ax,'TickLength',[0 0])
    % set the linewidth
    set(ax,'LineWidth',fig_set(ax_count).LineWidth*font_factor)
    % apply the box setting
    box(ax,fig_set(ax_count).box)
    % apply the labels
    xlabel(ax,fig_set(ax_count).XLabel,'FontName',fig_set(ax_count).FontName)
    ylabel(ax,fig_set(ax_count).YLabel,'FontName',fig_set(ax_count).FontName)
    title(ax,fig_set(ax_count).Title,'FontName',fig_set(ax_count).FontName)

    %% Colorbar

    if fig_set(ax_count).colorbar
        cba = colorbar;
        set(cba,'LineWidth',fig_set(ax_count).LineWidth*font_factor,'TickLength',0)
        ylabel(cba,fig_set(ax_count).colorbar_label,'FontSize',font_size_target*font_factor)
    end

    %% Apply the font size
    set(ax,'FontSize',font_size_target*font_factor,'FontName',fig_set(ax_count).FontName)
    %% Apply colormap

    colormap(fig_set(ax_count).cmap)
end
%% Scaling

% set the units to centimeters
set(fig_handle,'Units','centimeters')

% get the required font size and final figure dimensions
fig_size_target = fig_set(1).fig_size.*font_factor;

% get the current figure size
current_size = fig_handle.Position;

% if only one dimension was given, assume it's the width
if length(fig_size_target) == 1
    % get the proportional height
    fig_size_target = [fig_size_target,current_size(4)*fig_size_target/current_size(3)];
end

% modify the current size
current_size([3 4]) = fig_size_target;

% set the new figure dimensions
set(fig_handle,'Position',current_size)
%% Saving

set(gcf,'Color','w')
% set up the saving
fig_path = fig_set(1).fig_path;
fig_name = fig_set(1).fig_name;

file_path = fullfile(fig_path,fig_name);
if fig_set(1).crop && ~fig_set(1).painters
    handle_out = export_fig(file_path,'-r600');
elseif fig_set(1).painters
    handle_out = export_fig(file_path,'-r600','-painters');
else
    handle_out = export_fig(file_path,'-r600','-nocrop');
end



