% -----------------------------------------------------------------
% Legacy exploratory analysis by a member of the Wilson Lab.
% Historical record only — not part of the bancpipeline release.
% -----------------------------------------------------------------

%% Script to load matrices with connectivity and anatomical cluster information
% 1) get data directory and load data
% 2) reshape vector

% get data directory
datadir = [strrep( which("dn_anatomy_vs_connectivity.m"), ...
    ['analysis', filesep, 'matlab', filesep, 'dn_anatomy_vs_connectivity.m'], ''), ...
    'data', filesep, 'dn', filesep];

% define figure directory
figdir = [strrep( which("dn_anatomy_vs_connectivity.m"), ...
    ['analysis', filesep, 'matlab', filesep, 'dn_anatomy_vs_connectivity.m'], ''), ...
    'figures', filesep, 'dn_connectivity_dist'];

% load data
files2load = {'dn_input_cos_dist', 'dn_joint_cos_dist', 'dn_output_cos_dist'};
cos_dist_all = [];
dn_types_rows_all = [];
dn_types_columns_all = [];
inputfiletype = 2;

for i = 1:numel(files2load)

    % load file
    if inputfiletype == 1
        load([datadir, files2load{i}])
        matrix_siz = length(cos_dist)^0.5*[1 1];
        cos_dist_all{i} = flip(reshape(cos_dist, matrix_siz), 1);
    else
        load([datadir, files2load{i}], 'dn_types_rows', 'dn_types_columns')
        cos_dist = readtable([datadir, files2load{i}, '.csv']);
        cos_dist = table2array(cos_dist(2:end, 2:end));
        cos_dist_all{i} = flip(cos_dist, 1);
        matrix_siz = size(cos_dist);
    end
    % reshape from vector to matrix
    dn_types_rows_all{i} = flip(reshape(dn_types_rows, matrix_siz), 1);
    dn_types_columns_all{i} = flip(reshape(dn_types_columns, matrix_siz));
    clear cos_dist dn_types_rows dn_types_columns

end

%% 1) check the distribution of distances
% Note: basically values in the diagonal are < 10^-3

x_bins = -0.1:0.00001:1;
figname = {'input_cos_dist_hist', 'joint_cos_dist_hist', 'output_cos_dist_hist'};

for i = 1:numel(cos_dist_all)
    plot_hist_cos_dist(x_bins, cos_dist_all{i}, figname{i}, figdir)
end

clear figname x_bins

%% 2.1) plot heatmap of dn_types and cosine distances from 
%   vectorized matrix in different colormaps

% plot log scale
figname = {'input_cos_dist_heatmap_log', 'joint_cos_dist_heatmap_log', ...
    'output_cos_dist_heatmap_log'};

for i = 1:numel(cos_dist_all)
    plot_heatmap_cos_dist(cos_dist_all{i}, dn_types_rows_all{i}(:, 1), ...
        'log10', figname{i}, figdir)
end

% plot linear scale
figname = {'input_cos_dist_heatmap_linear', 'joint_cos_dist_heatmap_linear', ...
    'output_cos_dist_heatmap_linear'};

for i = 1:numel(cos_dist_all)
    plot_heatmap_cos_dist(cos_dist_all{i}, dn_types_rows_all{i}(:, 1), ...
        'linear', figname{i}, figdir)
end

clear figname

%% 2.2) plot heatmap of dn_types and cosine distances from 
%   matrix in different colormaps

% plot linear scale
figname = {'input_cos_dist_heatmap_linear_mat', 'joint_cos_dist_heatmap_linear_mat', ...
    'output_cos_dist_heatmap_linear_mat'};

for i = 1:numel(cos_dist_all)
    plot_heatmap_cos_dist(cos_dist_all{i}, dn_types_rows_all{i}(:, 1), ...
        'linear', figname{i}, figdir)
end

clear figname

%% 2.3) plot heatmap of dn_types and cosine distances from 