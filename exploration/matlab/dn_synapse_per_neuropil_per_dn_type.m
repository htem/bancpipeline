% -----------------------------------------------------------------
% Legacy exploratory analysis by a member of the Wilson Lab.
% Historical record only — not part of the bancpipeline release.
% -----------------------------------------------------------------

%% Script to load matrices with synapses per dn type in vnc and brain
% neuropils
% 1) get data directory and load data
% 2) reshape vector

% get data directory
datadir = [strrep( which("dn_anatomy_vs_connectivity.m"), ...
    ['analysis', filesep, 'matlab', filesep, 'dn_anatomy_vs_connectivity.m'], ''), ...
    'data', filesep, 'dn_synapses', filesep];

% define figure directory
figdir = [strrep( which("dn_anatomy_vs_connectivity.m"), ...
    ['analysis', filesep, 'matlab', filesep, 'dn_anatomy_vs_connectivity.m'], ''), ...
    'figures', filesep, 'dn_synapse_neuropil_per_dn_type'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

% load data
files2load = {'combine_by_dn_type', 'combine_by_dn_type_norm', 'combine_by_dn_type_num'};
syn_per_np = [];
dn_types = [];
np_names = [];

for i = 1:numel(files2load)

    % load file
    temp_matrix = readtable([datadir, files2load{i}, '.csv']);
    dn_types{i} = table2array(temp_matrix(1:end, 1));
    np_names{i} = temp_matrix.Properties.VariableNames(2:end);
    syn_per_np{i} = table2array(temp_matrix(1:end, 2:end));
    clear temp_matrix

end

%% 1) plot heatmap of synapse count per neuropil for each dn_type

% plot log scale with Court np
figname = {'dn_syn_per_dn_type_heatmap_log_court_mat', 'dn_normsyn_per_dn_type_heatmap_linear_court_mat', ...
    'dn_neuronn_per_dn_type_heatmap_linear_court_mat'};
scalestr = {'log10', 'linear', 'linear'};

for i = 1:numel(syn_per_np)
    plot_heatmap_syn_per_np(syn_per_np{i}, dn_types{i}, ...
        np_names{i}, scalestr{i}, 'court', figname{i}, figdir)
end

% plot log scale with MANC np
figname = {'dn_syn_per_dn_type_heatmap_log_manc_mat', 'dn_normsyn_per_dn_type_heatmap_linear_manc_mat', ...
    'dn_neuronn_per_dn_type_heatmap_linear_manc_mat'};
scalestr = {'log10', 'linear', 'linear'};

for i = 1:numel(syn_per_np)
    plot_heatmap_syn_per_np(syn_per_np{i}, dn_types{i}, ...
        np_names{i}, scalestr{i}, 'manc', figname{i}, figdir)
end

clear figname

% plot types as heatmap
% plot each heatmap (vnc/brain)
% add labels
% and keep   the same range for both