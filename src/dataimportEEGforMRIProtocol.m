%% Reads EEG .mff files, downsamples to 125 Hz and renames files so they have same names as BIDS fmri dataset
%
%   Lara Rösler

%--------------------------------------------------------------------------
%add eeglab path
addpath(genpath('/src/toolboxes/eeglab13_5_4b/'));
eeglab;

% define input
subjects = 1:32;

% define
upperdir = '/EEG_data/';

%--------------------------------------------------------------------------
% loop through subjects
for sub = 1:length(subjects)
    
    sub_id = ['S' num2str(subjects(sub), '%02.f')]; % e.g. S01
    
    % loop through sessions
    for session = [1:2]
        
        % enter subject directory
        filepath = [upperdir sub_id];
        cd(filepath);
        
        % list all files in directory
        filenames = dir(['*']);
        
        %-----------------------------------------------------------------
        %convert and resample data
        
        for file_idx = 3:size({filenames.name},2)
            filename = filenames(file_idx).name;
            
            fprintf('\nProcessing %s.\n\n', filename);
            
            EEG = pop_readegimff(sprintf('%s%s', filepath, filename));
            EEG = eeg_checkset(EEG);
            
            % re-sampling data to 125 Hz
            fprintf('\nre-sampling data to %s.\n\n', freq);
            
            freq = 125;
            EEG  = pop_resample(EEG, freq);
            EEG  = eeg_checkset(EEG);
            
            % renaming markers
            fprintf('Renaming markers.\n');
            for e = 1:length(EEG.event)
                if strcmp(EEG.event(e).type, 'TONE')
                    EEG.event(e).type = ['Hz' num2str(EEG.event(e).codes{2,2})];
                end
            end % e
            %--------------------------------------------------------------
           
            % rename output file
            % read out participant nr and session
            index  = strfind(filepath,'sourcedata/S')+12; %+2 because that's length of strfind
            sub_nr = filepath(index:index+1);
            
            index_ses = strfind(filepath,'session')+7;
            sess_nr   = filepath(index_ses);
            
            % sort data according to run
            % adopt same naming per file as BIDS fmri data
            
            if      ~isempty(strfind(filename,'a_sleep'))
                filename_new=sprintf('sub-S%s_ses-%s_task-SleepAscending',...
                                        sub_nr, sess_nr);
                out_path='/derivatives/imported_eeg/ascending/';
            elseif  ~isempty(strfind(filename,'d_sleep'))
                filename_new=sprintf('sub-S%s_ses-%s_task-SleepDescending',...
                                        sub_nr, sess_nr);
                out_path='/derivatives/imported_eeg/descending/';
            elseif  ~isempty(strfind(filename,'RSER_sleep'))
                filename_new=sprintf('sub-S%s_ses-%s_task-SleepRSER',...
                                        sub_nr, sess_nr);
                out_path='/derivatives/imported_eeg/rser/';
            elseif  ~isempty(strfind(filename,'a_wake'))
                filename_new=sprintf('sub-S%s_ses-%s_task-WakeAscending',...
                                        sub_nr, sess_nr);
                out_path='//derivatives/imported_eeg/ascending/';
            elseif  ~isempty(strfind(filename,'d_wake'))
                filename_new=sprintf('sub-S%s_ses-%s_task-WakeDescending',...
                                        sub_nr, sess_nr);
                out_path='/derivatives/imported_eeg/descending/';
            elseif  ~isempty(strfind(filename,'RSER_wake'))
                filename_new=sprintf('sub-S%s_ses-%s_task-WakeRSER',...
                                        sub_nr, sess_nr);
                out_path='/derivatives/imported_eeg/rser/';
            else
                filename_new=sprintf('sub-S%s_ses-%s_task-MRS',...
                                        sub_nr, sess_nr);
                out_path='/derivatives/imported_eeg/mrs/';
            end
            
            EEG.setname  = filename_new;
            EEG.filename = filename_new;
            EEG.filepath = out_path;
            
             %-------------------------------------------------------------
            fprintf('Saving %s%s.\n', EEG.filepath, EEG.filename);
            pop_saveset(EEG,'filename', EEG.setname, 'filepath', out_path,'version','7.3');
            
        end
    end
end
