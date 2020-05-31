function cmap_out = tint_colormap(cmap_in,tint_factor)
% tint the value in a colormap by the given tint factor

% allocate memory for the output
cmap_out = zeros(size(cmap_in));

% for all the values
for val = 1:size(cmap_in,1)
    cmap_out(val,:) = cmap_in(val,:) + [tint_factor*(1-cmap_in(val,1)),...
        tint_factor*(1-cmap_in(val,2)),tint_factor*(1-cmap_in(val,3))];
    
end