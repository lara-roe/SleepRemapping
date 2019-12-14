%% Create scan protocol for fMRI analysis so that each TR is associated with tone frequency per run

%--------------------------------------------------------------------------
clear all;
addpath('/home/lr499/bin/eeglab13_5_4b/');
eeglab;

outpath = '/derivatives/regressors/';

% run types

runs = {'Ascending','Descending'}

%--------------------------------------------------------------------------
for r_idx = 1:length(runs)
    
    run = runs{r_idx};
    
    % cd to where the imported files are stored
    files_dir = '/derivatives/imported_eeg/';
    files_dir_run = [files_dir run];
    
    cd(files_dir_run)
    
    % list files
    filenames = dir(fullfile(cd, '*.set'));
    filenames = {filenames.name};
    
    %---------------------------------------------------------------------
    % read in file
    
    for idx = 1:length(filenames)
        filename = filenames{idx};
        EEG = pop_loadset('filename', filename, 'filepath', cd);
        EEG = eeg_checkset(EEG);
     
    %---------------------------------------------------------------------    
        % find indexes of TRs
        tr_idx =find(strcmp({EEG.event.type},'TR')==1);
        tr_nr = length(tr_idx);
        
        stim_idx = find(not(cellfun('isempty', strfind({EEG.event.type},'Hz'))));
        stim_nr = length(stim_idx);
        
        first_stim = min(stim_idx);
        last_stim = max(stim_idx);
        
        first_rel_tr = max(find(tr_idx<first_stim));
        last_rel_tr = find(tr_idx>last_stim); 
        
         
    %---------------------------------------------------------------------
        % create vector with hz values
        if strcmp(run,'Descending')==1
        all_tones_cell = {'Hz4000' 'Hz2828' 'Hz2000' 'Hz1414' 'Hz1000' ...
                          'Hz707' 'Hz500' 'Hz354' 'Hz250' 'Hz177' 'Hz125'};
        all_tones = [4000 2828 2000 1414 1000 707 500 354 250 177 125];
        else
        all_tones_cell = {'Hz125' 'Hz177' 'Hz250' 'Hz354' 'Hz500' 'Hz707'...
                            'Hz1000' 'Hz1414'  'Hz2000'  'Hz2828' 'Hz4000'};
        all_tones = [125 177 250 354 500 707 1000 1414 2000 2828 4000];
        end
        
    %---------------------------------------------------------------------
        % create vector that contains stimulations 
        %(30 repetitions of 11 tones, starting at onset)
        stim_vec = NaN(1,tr_nr);
        nr_repetitions = 30;
        increments = 13;
        onset_markers = first_rel_tr + (0:nr_repetitions-1)*increments;
        for onset_idx = 1:length(onset_markers)
            onset = onset_markers(onset_idx);
            stim_vec(onset:(onset+10)) = all_tones;
        end
    %---------------------------------------------------------------------
        % save TRs per tone in structure
        for i = 1:length(all_tones)
           out.(char(all_tones_cell(i)))= find(ismember(stim_vec,all_tones(i)));
        end
        
        % save output
        out_name = [outpath filename(1:end-4) '_reg.mat'];
        save(out_name,'out');
        
    end
end