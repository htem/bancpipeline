% -----------------------------------------------------------------
% Legacy exploratory analysis by a member of the Wilson Lab.
% Historical record only — not part of the bancpipeline release.
% -----------------------------------------------------------------

function [column_name, scores2use, idsfound] = extract_influence_bodyparts(filename, foldername, neuron_id)
% extract_influence_bodyparts: function to read csv files containing
%   influence scores for body parts
%
% Usage:
%   extract_influence_bodyparts(filename, foldername, neuron_id)
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

% set deafult order of sensory influences:
bodyparts_ordered = {'JO', ...
    'eye_bristle', ...
    'head_bristle', ...
    'ocellar', ...
    'neck', ...
    'brain_endocrine', ...
    'gustatory_brain', ...
    'ciberial_mechanosensory', ...
    'pharynx', ...
    'enteric_gustatory_brain', ...
    'prothoracic_chordotonal_organ', ...
    'front_leg', ...
    'middle_leg', ...
    'hind_leg', ...
    'wing', ...
    'haltere', ...
    'notum', ...
    'abdomen'};

[~, idx] = ismember(bodyparts_ordered, column_name);

scores2use = scores2use(:, idx);
column_name = column_name(idx);
column_name = strrep(column_name, '_', '-');

end