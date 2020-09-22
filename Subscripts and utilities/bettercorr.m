function fig_handle = bettercorr(rho,cmap,varargin)
% function to plot symmetric matrices with written values in one half

% if a figure handle is given, use that, otherwise create a new one
if length(varargin) >= 1 && ~isempty(varargin{1})
    fig_handle = varargin{1};
    figure(fig_handle)
else
    fig_handle = figure;
end

% if a range for the colors is given, use that, otherwise use -1 to 1
if length(varargin) >= 2 && ~isempty(varargin{2})
    range = varargin{2};
else
    range = [-1 1];
end

% if significance is provided, add asterisks on the significant squares
if length(varargin) >= 3 && ~isempty(varargin{3})
    significance_matrix = varargin{3};
else
    significance_matrix = [];
end

% take the upper half of the matrix
upper = triu(rho,1);
% make the other values NaN
upper(upper==0) = NaN;
% plot the matrix
imagesc(upper)
hold on
% get the size of the correlation matrix
corr_size = size(rho);
% get the coordinates for text as a meshgrid
[x,y] = meshgrid(1:corr_size(1));
% get the lower triangle only (for the grid and the correlation)
x = tril(x,-1);
y = tril(y,-1);
rho_num = tril(rho,-1);
% take only the non-zero components
x_plot = x(x~=0);
y_plot = y(x~=0);
rho_plot = rho_num(x~=0);


% for all the nonzeros points
for points = 1:length(rho_plot)
    % interpolate the color
    color = interp1(linspace(range(1),range(2),256),cmap,rho_plot(points));
    % get the value to plot
    rho_val = num2str(round(rho_plot(points),2));
%     % trim the decimals (this needs to be cleaner)
%     rho_val = rho_val(1:4);
    % "plot" the correlation or whatever values
    text(x_plot(points),y_plot(points),rho_val,...
        'HorizontalAlignment','center','Color',color,...
        'FontSize',5,'FontWeight','bold')
    
    % if significance_matrix is not empty
    if ~isempty(significance_matrix)
        % if it's a significant point
        if significance_matrix(x_plot(points),y_plot(points)) < 0.05
%             text(y_plot(points),x_plot(points),'*',...
%                 'HorizontalAlignment','center','Color','k',...
%                 'FontSize',7,'FontWeight','bold','VerticalAlignment','cap')
            scatter(y_plot(points),x_plot(points),5,'k*')
        end
    end
end
