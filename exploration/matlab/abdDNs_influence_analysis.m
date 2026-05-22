% -----------------------------------------------------------------
% Legacy exploratory analysis by a member of the Wilson Lab.
% Historical record only — not part of the bancpipeline release.
% -----------------------------------------------------------------

%% Script to load and plot influence scores for selected cell types

% 1) define main directory
influence_dir = 'C:\Users\papers\BANC-project\data\influence\frankenbrain_v1.5\Use_this';

% 2) get list of csv files
% pull sensory influence by body part
bodypartsen_files = rdir([influence_dir, filesep, 'Rachel_request', filesep, '*']);
[filename_bodypart, foldername_bodypart] =  split_path({bodypartsen_files.name});

% pull sensory influence by modality
modality_files = rdir([influence_dir, filesep, 'Modality', filesep, '*']);
[filename_mod, foldername_mod] =  split_path({modality_files.name});

% cell types of interest
neuron_celltype = {'oviDNa_a', 'oviDNa_a', ...
    'oviDNa_b', 'oviDNa_b', 'oviDNa_b', 'oviDNa_b', ...
    'oviDNb', 'oviDNb', ...
    'DNpe034', 'DNpe034', ...
    'DNpe044', 'DNpe044', ...
    'DNpe046', 'DNpe046', ...
    'DNpe047', 'DNpe047'};

neuron_id = {'720575941651301909', '720575941643911032', ...
    '720575941599689257', '720575941563681799', '720575941539873917', '720575941654228180', ...
    '720575941633146043', '720575941687810828', ...
    '720575941584744926', '720575941503840306', ...
    '720575941678731197', '720575941542904325', ...
    '720575941512235075', '720575941515241465', ...
    '720575941641360722', '720575941553339034'};

% load and pull cell type of interest
[column_name_bodypart, scores2use_bodypart, idsfound_bodypart] = extract_influence_bodyparts(filename_bodypart, foldername_bodypart, neuron_id);

% load and pull cell type of interest
[column_name_mod, scores2use_mod, idsfound_mod] = extract_influence_modality(filename_mod, foldername_mod, neuron_id);

% get cell types and sort ids based on cell type
for i = 1:numel(idsfound_bodypart)
    idx2use = contains(neuron_id, idsfound_bodypart(i));
    if sum(idx2use) > 0
        row_name(i) = neuron_celltype(idx2use);
    end
end

% plot scores bodyparts
figH = figure('Position', [100 100 1300 450]);
   
axH(1) = subplot(1, 1, 1);
%imagesc(1 - log(scores2use_bodypart))
imagesc(scores2use_bodypart)

cbarH = colorbar;

axH(1).YTick = 1:size(scores2use_bodypart, 1);
axH(1).YTickLabel = row_name;
axH(1).YLabel.String = 'Cell types';
axH(1).XTick = 1:size(scores2use_bodypart, 2);
axH(1).XTickLabel = column_name_bodypart;
axH(1).XTickLabelRotation = 45;
axH(1).XLabel.String = 'Sensory influence';
axH(1).XAxisLocation = 'top';

cbarH.Position = [0.93 0.4 0.02 0.15];

% plot scores bodyparts
figH = figure('Position', [100 100 1300 450]);
   
axH(1) = subplot(1, 1, 1);
%imagesc(1 - log(scores2use_mod))
imagesc(scores2use_mod)

cbarH = colorbar;

axH(1).YTick = 1:size(scores2use_mod, 1);
axH(1).YTickLabel = row_name;
axH(1).YLabel.String = 'Cell types';
axH(1).XTick = 1:size(scores2use_mod, 2);
axH(1).XTickLabel = column_name_mod;
axH(1).XTickLabelRotation = 45;
axH(1).XLabel.String = 'Sensory influence';
axH(1).XAxisLocation = 'top';

cbarH.Position = [0.93 0.4 0.02 0.15];