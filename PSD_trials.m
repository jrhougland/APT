
%% 2026-03 APT - PSD for each trial %%

% calculate PSD for each trial from cleaned/preprocessed pre-stim EEG

%% Initial Settings 

clear
close all
clc

% eeglab path
eeglab_path = 'C:/Program Files/MATLAB/R2025a/eeglab_current';
addpath(genpath(eeglab_path));
eeglab nogui;

addpath('\\storage.neurologie.uni-tuebingen.de\bbnp_lab\Projects\2024-09 HOUGLANDPHD\Toolboxes\Neurone')

addpath('\\storage.neurologie.uni-tuebingen.de\bbnp_lab\Projects\2024-09 HOUGLANDPHD\Toolboxes')

% path to data directory
dataset = 'DLPFC';     % M1 or DLPFC or DLPFC Pilot
path.load = ['E:\APT\' dataset '\'];
path.save = ['E:\APT\' dataset '\'];

% identify subject id lists and electrodes of interest
if strcmp(dataset, 'M1')

    subj = [18,19,20,21,22,23,24,25,26,27,28,29,31,32,34,35,36,38,39,40,41,42,43,44,45,46,48,49,50];    % REFTEP
    subj_filename = 'REFTEP_';
    test.electrodes = {'C3', 'C1', 'FC3', 'FC1'};  % M1

elseif strcmp(dataset, 'DLPFC')

    subj = [1,2,3,4,6,7,8,9,10];  %FRONTEP
    subj_filename = 'FRONTEP_';
    test.electrodes = {'F3', 'F1', 'FC3', 'F5', 'FC1', 'FC5'}; % DLPFC

else

    subj = [2,3];  % FRONTEP_pilot
    subj_filename = 'FRONTEP_pilot';
    test.electrodes = {'F3', 'F1', 'FC3', 'F5', 'FC1', 'FC5'}; % DLPFC
    path.load = ['E:\APT\DLPFC\'];
    path.save = ['E:\APT\DLPFC\'];

end



% to get PSD for all electrodes
get_all_chans = 0;  % 1 for yes, 0 for no

formatSpec = '%03.0f';

%%

for id = subj

  % start notification
  fprintf('\nStarting participant %03.0f. \n', id)
  disp('Start loading...')

%% Check / create directories

 if ~exist([path.save '\' num2str(id,formatSpec)], 'dir')
      mkdir([path.save '\' subj_filename num2str(id,formatSpec)])
 end

 %% Load data

 EEG = pop_loadset('filename', [subj_filename num2str(id, formatSpec) '_EEG_pre_processed.set'], ...
                      'filepath', [path.load subj_filename num2str(id, formatSpec) '\']);

 

%%  Find electrodes of interest

    if get_all_chans == 1

        % Get all electrode labels
        test.electrodes = {EEG.chanlocs.labels};    % Get all channels

    else 

    end

    % Convert struct to cell (labels are in first row)
    electrode_labs = struct2cell(EEG.chanlocs);

    % Find indices
    for i = 1:length(test.electrodes)
        test.electrode_idx(i) = find(strcmpi(test.electrodes{i}, electrode_labs(1,:)));
    end


%% PSD Calculation

    Fs      = EEG.srate;
    nTrials = EEG.trials;
    nChan   = length(test.electrode_idx);
    
    f_range = [4 50];
    
    win_len  = round(1 * Fs);
    noverlap = round(0.50 * win_len);
    window   = hanning(win_len);
    nFFT     = win_len;
    
    % Get frequency vector once
    x0 = double(squeeze(EEG.data(test.electrode_idx(1),:,1)));
    [Pxx0, f] = pwelch(x0, window, noverlap, nFFT, Fs);
    
    idx   = f >= f_range(1) & f <= f_range(2);
    freqs = f(idx);
    nFreq = length(freqs);
    
    % Preallocate: trials × selected_channels × freqs
    psd = zeros(nTrials, nChan, nFreq);
    
    for tr = 1:nTrials
        for ci = 1:nChan
    
            ch = test.electrode_idx(ci);
    
            x = double(squeeze(EEG.data(ch,:,tr)));
            [Pxx, ~] = pwelch(x, window, noverlap, nFFT, Fs);
    
            psd(tr, ci, :) = Pxx(idx);
    
        end
    end
    
    save([path.save '\' subj_filename num2str(id,formatSpec) '\' subj_filename num2str(id,formatSpec) '_eeg_trials_psd.mat'], 'psd', 'freqs');
    
    
    
   % clear EEG 

end

  
 
    
