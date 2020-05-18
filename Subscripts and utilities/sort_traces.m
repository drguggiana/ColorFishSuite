function trace_out = sort_traces(trace_in,varargin)
% use hierarchical clustering to sort traces

% if a cluster number was specified
if length(varargin) >= 1
    max_clust = varargin{1};
else
    max_clust = 10;
end

% get the linkage tree
T = linkage(trace_in,'complete','euclidean');

% get the indexes
idx_clu = cluster(T,'maxclust',max_clust);

% output the sorted traces
[~,sort_idx] = sort(idx_clu);
trace_out = trace_in(sort_idx,:);

% if the plotting option is on
if length(varargin) >= 2
    if varargin{2}
        figure
        dendrogram(T)
    end
end
    