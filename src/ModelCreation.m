%% Here we create a model based on the convolution of stimulus timing and spm HRF function

clear all;

%% Prepare all paths and specify which runs should be analyzed
% path to readcsv.m

addpath(genpath('/src/toolboxes/'));

reg_dir = '/derivatives/regressors/';
fmridir =  '/fmriprep/'; 
outpath = '/derivatives/models/';

sublist=dir(fullfile(fmridir,'sub*'));
dirFile   = [sublist.isdir];
sublist = {sublist(dirFile).name};

% define number of sessions and task
scode = 'ses-1';
taskids = {'WakeAscending','WakeDescending','SleepDescending','SleepAscending'};  

% loop through subjectss
for sub = [7:32]
    
    % loop through tasks
    for taskidx = 1:length(taskids)
        taskid = taskids{taskidx};
        
        % read in events file
        reg_mat = [reg_dir sublist{sub},'_',scode,'_task-',taskid, '_reg.mat'];
        load(reg_mat);
        labels = fieldnames(design);
        
        % find largest TR value in events file
        maxTR = max(structfun(@(x)max(x(:)),design));  
      
        % create model structure in which we save loop output (copied from
        % design because that's desired length)
        model_struct = design;
        
        % loop through all labels (different tone conditions)
        for label_idx = 1:length(labels)
          
        % create vector of zeros which has same lengths as TRs
        time_vec = zeros(maxTR,1);
      
        design_vec = rmmissing(design.(labels{label_idx})); % remove NAs
        time_vec(design_vec) = 1; %replace Tone TRs with 1s (+1 because first TR in design vec labelled as 0 but in other scripts will be 1)
        conv_vec = conv(spm_hrf(2),time_vec); % convolve HRF with time vec
        model_struct.(labels{label_idx}) = conv_vec; % save in model struct
        
        end
        
        outputname = [outpath sublist{sub},'_',scode,'_task-',taskid, 'ModelStruct' ];
        save(outputname, '-struct', 'model_struct')
    end
end
