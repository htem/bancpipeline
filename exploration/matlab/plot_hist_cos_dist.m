% -----------------------------------------------------------------
% Legacy exploratory analysis by a member of the Wilson Lab.
% Historical record only — not part of the bancpipeline release.
% -----------------------------------------------------------------

function plot_hist_cos_dist(x_bins, input_dist_matrix, figname, figdirectory)
% plot_hist_cos_dist: function to plot histogram of cosine distances
%
% Usage:
%   plot_vel_hist_panels(x_bins, input_dist_matrix, figname, figdirectory)
%
% Args:
%   x_bins: histogram bins
%   input_dist_matrix: input vectors
%   figname: figure name
%   figdirectory: target directory

diag_matrix = diag(input_dist_matrix);

% generate histograms & normalize histograms to unity
all_hist = hist(input_dist_matrix(:), x_bins);
diag_hist = hist(diag_matrix(:), x_bins);

all_hist = all_hist./sum(all_hist, 2);
diag_hist = diag_hist./sum(diag_hist, 2);

xval = {x_bins, x_bins, x_bins, x_bins};
yval = {all_hist, diag_hist, all_hist, diag_hist};

% default parameters
% display y axis
vel_yscale = 'linear';
vel_yscale_log = 'log';
xlim_per_vel = {[-0.1 1.1], [-0.1 1.1], [-0.1 1.1], [-0.1 0.1]};
ylim_per_vel = {[], [], [0 10^-4], [0 0.01]};
ylim_per_vel_log = {[10^-5 1], [10^-5 1], [10^-5 10^-4], [10^-5 1]};
xlabel_str = {'cosine distance', 'cosine distance', 'cosine distance', 'cosine distance'};
ylabel_str_a = {'probability 0-1 (linear)', 'probability 0-1 (linear)', ...
    'probability 0-1 (linear)', 'probability 0-1 (linear)'};
ylabel_str_b = {'probability 0-1 (log)', 'probability 0-1 (log)', ...
    'probability 0-1 (log)', 'probability 0-1 (log)'};
title_str = {'full matrix', 'diagonal only', 'full matrix zoom in y', 'diagonal only zoom in x'};
line_color_left = [1 0 0];
line_color_right = [.5 .5 .5];

figH(1) = figure('Position', [100 100 1000 500]);

% plot yaw/side/forward velocity histograms
for i = 1:numel(yval)
    axH(i) = subplot(2, 2, i);
    plot(xval{i}, yval{i}, '-', 'Color', line_color_left)
    hold(axH(i), 'on')
    yyaxis(axH(i), 'right')
    plot(xval{i}, yval{i}, '-', 'Color', line_color_right)
end

% edit axis
for i = 1:numel(yval)
    % edit X axis
    axH(i).XLabel.String = xlabel_str{i};
    axH(i).XLim = xlim_per_vel{i};
    % edit Y axis scale
    axH(i).YAxis(1).Scale = vel_yscale;
    axH(i).YAxis(1).Label.String = ylabel_str_a{i};  
    axH(i).YAxis(2).Label.String = ylabel_str_b{i};  
    if ~isempty(ylim_per_vel{i})
        axH(i).YAxis(1).Limits = ylim_per_vel{i};
    end
    axH(i).YAxis(2).Scale = vel_yscale_log;
    axH(i).YAxis(2).Limits = ylim_per_vel_log{i};
    % edit color
    axH(i).YAxis(1).Color = line_color_left;
    axH(i).YAxis(2).Color = line_color_right;
    axH(i).Title.String = title_str{i};
end

if ~exist(figdirectory, 'dir')
    mkdir(figdirectory)
end

figname = [figdirectory, filesep, figname];
saveas(figH, [figname, '.fig'])
saveas(figH, [figname, '.png'])
close(figH)

end