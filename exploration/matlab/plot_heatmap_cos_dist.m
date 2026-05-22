% -----------------------------------------------------------------
% Legacy exploratory analysis by a member of the Wilson Lab.
% Historical record only — not part of the bancpipeline release.
% -----------------------------------------------------------------

function plot_heatmap_cos_dist(input_dist_matrix, clusters_rows, zscale, figname, figdirectory)
% plot_heatmap_cos_dist: function to plot heatmaps cosine distances
%
% Usage:
%   plot_heatmap_cos_dist(input_dist_matrix, figname, figdirectory)
%
% Args:
%   x_bins: histogram bins
%   input_dist_matrix: input vectors
%   zscale: scale to use for input distance matrix
%   figname: figure name
%   figdirectory: target directory

% from the distribution it looks like any value below 0.001 is ~= 0 or the
%   lowest distance
dist_cap = 10^-3;
input_dist_matrix(input_dist_matrix < dist_cap) = dist_cap;

if strcmp(zscale, 'log10')
    input_dist_matrix = log10(input_dist_matrix);
    % set default log settings
    icaxis = {[], log10([10^-3 1]), [],  [log10([10^-0.1 1])]};
    colorbarstr = {[], 'log10(cosine distance)', [], 'log10(cosine distance)'};
    colorbarticks = {[], log10([10^-3 10^-2 10^-1 1]), [], log10([10^-0.1 10^-0.05 1])};
elseif strcmp(zscale, 'linear')
    % set default log settings
    icaxis = {[], [0 1], [],  [0.8 1]};
    colorbarstr = {[], 'cosine distance', [], 'cosine distance'};
    colorbarticks = {[], 0:0.2:1, [], 0.8:0.05:1};
end

% generate YTicks for cluster heatmap
[unique_clusters, idx_vector] = unique(clusters_rows);

[unique_clusters, sort_idx] = sort(unique_clusters, 'descend');
idx_vector = idx_vector(sort_idx);

% generate colormaps
colormap_cos_dist = colorGradient(rgb('blue'), rgb('white'), 254);
colormap_clus = prism(length(unique_clusters));

imatrix = {clusters_rows, input_dist_matrix, clusters_rows, input_dist_matrix};
yvars = {1:length(clusters_rows), [], 1:length(clusters_rows), []};
icolormap = {colormap_clus, colormap_cos_dist, colormap_clus, colormap_cos_dist};

figH = figure('Position', [100 100 800 800]);

for i = 1:numel(imatrix)
    
    axH(i) = subplot(2, 2, i);

    if ~isempty(yvars{i})
        imagesc(1, yvars{i}, imatrix{i})
    else
        imagesc(imatrix{i})
    end

    % set colormap
    colormap(axH(i), icolormap{i})
    if ~isempty(icaxis{i})
        caxis(axH(i), icaxis{i});
    end
    
    % add colorbar
    if ~isempty(colorbarstr{i})
        cbarH(i) = colorbar;
        cbarH(i).Label.String = colorbarstr{i};
        cbarH(i).Ticks = colorbarticks{i};
    end

    if ~isempty(yvars{i})
        axH(i).YTick = idx_vector;
        axH(i).YTickLabel = chunk2cell(unique_clusters, 1);
        axH(i).XTick = [];
        axH(i).YLabel.String = ['anatomical clusters (1-', num2str(max(unique_clusters)), ')'];
    else
        axH(i).XTick = [];
        axH(i).YTick = [];
        axH(i).XLabel.String = ['anatomical clusters (1-', num2str(max(unique_clusters)), ')'];
    end

end

% edit panel sizes
axH(1).Position = [0.1 0.55 0.03 0.4];
axH(2).Position = [0.15 0.55 0.4 0.4];
axH(3).Position = [0.1 0.1 0.03 0.4];
axH(4).Position = [0.15 0.1 0.4 0.4];

% edit size of colorbar
cbarH(2).Position = [0.56 0.75 0.02 0.15];
cbarH(4).Position = [0.56 0.25 0.02 0.15];

figname = [figdirectory, filesep, figname];
saveas(figH, [figname, '.fig'])
saveas(figH, [figname, '.png'])
close(figH)

end