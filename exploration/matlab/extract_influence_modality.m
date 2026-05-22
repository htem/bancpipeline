% -----------------------------------------------------------------
% Legacy exploratory analysis by a member of the Wilson Lab.
% Historical record only — not part of the bancpipeline release.
% -----------------------------------------------------------------

function [column_name, scores2use, idsfound] = extract_influence_modality(filename, foldername, neuron_id)
% extract_influence_modality: function to read csv files containing
%   influence scores for all modalities
%
% Usage:
%   extract_influence_modality(filename, foldername, neuron_id)
%
% Args:
%   filename: csv file names
%   foldername: csv file directories
%   neuron_id: ids to find

% Initialize variables
column_name = [];
scores2use = [];
idsfound = {};

for i = 1:numel(filename)

    % Read the table
    table_temp = readtable([foldername{i}, filesep, filename{i}], 'Format','%f%s%f');

    % Process the first file for index extraction
    if i == 1
        % Find matching indices using ismember
        [isMatch, orderIdx] = ismember(neuron_id, table_temp.id);
        
        % Collect idsfound in the order of neuron_id
        idsfound = table_temp.id(orderIdx(isMatch));

        % Correctly order by the provided neuron_id, and ensure valid matches only
        index2use = isMatch;
    end
    
    % Store column names
    column_name{i, 1} = table_temp.Properties.VariableNames{3};
    
    % Sort scores according to the order of neuron_id
    scores2use(:, i) = table_temp{orderIdx(index2use), column_name{i, 1}};

end

% simplify column names
column_name = strrep(column_name, '_influence_norm_unsigned_forward_steady_state', '');
column_name = strrep(column_name, '_influence_norm_unsigned_forward_steady_', '');
column_name = strrep(column_name, '_influence_norm_unsigned_forward_s', '');
column_name = strrep(column_name, '_influence_norm_unsigned_forward_', '');

% fix some names:
column_name = strrep(column_name, '_neuronst', '_neuron');
column_name = strrep(column_name, 'proprioceptivete', 'proprioceptive');
column_name = strrep(column_name, 'tactilesta', 'tactile');
column_name = strrep(column_name, 'otherst', 'other');
column_name = strrep(column_name, 'ingestion_motor_neurons', 'ingestion_motor_neuron');
column_name = strrep(column_name, 'proprioceptive_tactiles', 'proprioceptive_tactile');
column_name = strrep(column_name, 'proboscis_motor_neurons', 'proboscis_motor_neuron');

% set deafult order of sensory influences:
modality_ordered = {'visual_achromatic', ...
    'visual_chromatic', ...
    'ocellar', ...
    'olfactory', ...
    'thermosensory', ...
    'hygrosensory', ...
    'chemosensory', ...
    'chemosensory_proprioceptive', ...
    'chemosensory_tactile', ...
    'tactile', ...
    'gustatory', ...
    'enteric_gustatory', ...
    'auditory', ...
    'wind_gravity', ...
    'mechanosensory', ...
    'johnstons_organ_other', ...
    'ciberial_mechanosensory_neuron', ...
    'proprioceptive', ...
    'proprioceptive_tactile', ...
    'sensory_ascending', ...
    'endocrine', ...
    'eye_motor_neuron', ...
    'antennal_motor_neuron', ...
    'neck_motor_neuron', ...
    'proboscis_motor_neuron', ...
    'grooming', ...
    'haltere', ...
    'ingestion_motor_neuron', ...
    'ascending', ...
    'descending', ...
    'mixed'};

[~, idx] = ismember(modality_ordered, column_name);

scores2use = scores2use(:, idx);
column_name = column_name(idx);
column_name = strrep(column_name, '_', '-');

end