%% 2026-03 APT - Source Reconstruction %%

% Source reconstruction of EEG data
% Inputs are REFTEP headmodels from realign_elctrodes_create_headmodels.m
% and cleaned EEG data (pre and post TMS)

%% Initial Settings 

clear
close all
clc

eeglab_path = 'C:/Program Files/MATLAB/R2025a/eeglab_current';
addpath(genpath(eeglab_path));
eeglab nogui;

addpath('\\storage.neurologie.uni-tuebingen.de\bbnp_lab\Projects\2024-09 HOUGLANDPHD\Toolboxes')
addpath('Z:\Projects\2024-09 HOUGLANDPHD\Toolboxes\Matti functions')
addpath('Z:\Projects\2024-09 HOUGLANDPHD\Toolboxes\hbf_distribution_open_v170624')
addpath('C:\Program Files\MATLAB\R2025a\fieldtrip-20260617')
addpath('Z:\Projects\2024-09 HOUGLANDPHD\Toolboxes\plotroutines_v170706')
addpath('Z:\Projects\2024-09 HOUGLANDPHD\Toolboxes\plotroutines_v180921')
addpath(genpath('Z:\Projects\2024-09 HOUGLANDPHD\Toolboxes\Source_JM'))

ft_defaults;

formatSpec = '%03.0f';

% path to data directory
dataset = 'M1';
headmodel_path.load = ['Z:\Projects\2026-03 APT\Analysis\' dataset '\source_paolo'];
eeg_path.load = ['E:\APT\' dataset];
source_path.save = ['Z:\Projects\2026-03 APT\Analysis\' dataset '\source_paolo'];

% Load your binary ROI mask (15684x1, ones = vertices of interest)
% C must be loaded/defined before the loop — adjust path as needed
load('Z:\Projects\2026-03 APT\Analysis\Scripts\handknob.mat');  % loads variable C
roi_indices = find(C == 1);   % vertex indices belonging to your ROI
fprintf('ROI contains %d vertices\n', length(roi_indices))

% subject list
subj = [18,19,20,21,22,23,24,25,26,27,28,29,31,32,34,35,36,38,39,40,41,42,43,44,45,46,48,49,50]; 


%%

for id = subj
    
    fprintf('Subject %d ', id)

    % load headmodel
    load([headmodel_path.load '\REFTEP_' num2str(id, formatSpec) '_headmodel'])

    % load cleaned EEG (pre and post TMS)
    EEG_pre = pop_loadset([eeg_path.load '\REFTEP_' num2str(id, formatSpec) '\REFTEP_' ...
        num2str(id, formatSpec) '_EEG_pre_processed.set']);

    EEG_post = pop_loadset([eeg_path.load '\REFTEP_' num2str(id, formatSpec) '\REFTEP_' ...
        num2str(id, formatSpec) '_EEG_post_processed.set']);
    

    %% LCMV Beamforming 

    % Match EEG channels to headmodel
    eeg_labels = {EEG_pre.chanlocs.labels}';
    [~, elec_indx] = ismember(eeg_labels, headmodel.label);
    
    if any(elec_indx == 0)
        warning('Removing %d channels not in headmodel', sum(elec_indx==0))
        EEG_pre.data(elec_indx == 0, :, :)  = [];
        EEG_post.data(elec_indx == 0, :, :) = [];
        eeg_labels(elec_indx == 0)           = [];
        [~, elec_indx] = ismember(eeg_labels, headmodel.label);
    end

    % Reorder the leadfield matrix to match the electrode order in the data
    L = headmodel.leadfield(elec_indx, :);
    headmodel.leadfield = L;
    headmodel.label = eeg_labels;

    % Average reference leadfield
    L = L - mean(L, 1);

    % Reshape EEG to [channels x all_samples] for covariance
    chN  = size(EEG_pre.data, 1);
    Xst  = reshape(EEG_pre.data, chN, []);

    % Covariance matrix
    Cov    = Xst * Xst' / size(Xst, 2);
    lambda = 10 * max(eig(Cov));
    invCy  = pinv(Cov + lambda * eye(size(Cov)));

    % Compute beamformer weights for every dipole
    n_dipoles  = size(L, 2);
    n_channels = size(L, 1);
    weights    = zeros(n_dipoles, n_channels);

    fprintf('Computing beamformer weights (%d dipoles)...\n', n_dipoles)
    for i = 1:n_dipoles
        lf = L(:, i);
        weights(i, :) = pinv(lf' * invCy * lf) * lf' * invCy;
    end

    % Subset weights to ROI vertices only
    weights_roi = weights(roi_indices, :);   % [n_roi_vertices x n_channels]


    %% Apply weights and extract ROI timecourse (mean over ROI vertices)

    n_times_pre   = size(EEG_pre.data, 2);
    n_trials_pre  = size(EEG_pre.data, 3);
    roi_tc_pre    = zeros(n_times_pre, n_trials_pre);   % [times x trials]

    n_times_post  = size(EEG_post.data, 2);
    n_trials_post = size(EEG_post.data, 3);
    roi_tc_post   = zeros(n_times_post, n_trials_post);

    % Pre TMS
    fprintf('Extracting pre-TMS ROI timecourse...\n')
    for trial = 1:n_trials_pre
        trial_data      = squeeze(EEG_pre.data(:, :, trial));       % [channels x times]
        source_tc_roi   = weights_roi * trial_data;                  % [n_roi_vertices x times]
        roi_tc_pre(:, trial) = mean(source_tc_roi, 1)';             % mean over vertices
    end

    % Post TMS
    fprintf('Extracting post-TMS ROI timecourse...\n')
    for trial = 1:n_trials_post
        trial_data      = squeeze(EEG_post.data(:, :, trial));       % [channels x times]
        source_tc_roi   = weights_roi * trial_data;                  % [n_roi_vertices x times]
        roi_tc_post(:, trial) = mean(source_tc_roi, 1)';            % mean over vertices
    end


    %% Save

    save([source_path.save '\REFTEP_' num2str(id, formatSpec) '_handknob_source_pre'], ...
        'roi_tc_pre', 'roi_indices', '-v7.3')

    save([source_path.save '\REFTEP_' num2str(id, formatSpec) '_handknob_source_post'], ...
        'roi_tc_post', 'roi_indices', '-v7.3')

    fprintf('Saved source reconstruction for Subject %d\n', id)

    clear Xst Cov invCy weights weights_roi L elec_indx source_tc_roi

end

fprintf('All subjects processed.\n')