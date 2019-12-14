%% Here we cross-correlate voxel timeseries and tone models 
% Highest correlation coefficient per voxel received label of tone
%
% Timeseries from fmriprep MNI normalized functional runs
%

clear all;
close all;
addpath(genpath('/src/toolboxes/'));

% read in subject list
fmridir =  'derivatives/fmriprep/'; 
sublist = dir(fullfile(fmridir,'sub*'));
dirFile = [sublist.isdir];
sublist = {sublist(dirFile).name};

for sub = 1:length(sublist)
  %-------------------------------------------------------------------------
  % initial definitions of paths, runs, slices ..
  
    model_path = '/derivatives/models/';
    func_path  = ['derivatives/fmriprep/' sublist{sub} '/ses-1/func/'];
    anat_path  = ['/derivatives/fmriprep/' sublist{sub} '/anat/'];
    roi_path   = '/derivatives/rois/';
    out_path   = '/derivatives/maps/';
    
    % define tasks to be analyzed
    taskids    = {'WakeAscending','SleepAscending', ...
                  'WakeDescending','SleepDescending'};
    
    %initialize empty output matrix
    all_coeffs = [];
    
    %define which slices to loop through (decided based on visual
    %inspection)
    slices = [19:26];
    
    %define start and end TR
    start_TR=25;
    end_TR=219;
    
    %----------------------------------------------------------------------
      
    for tidx = 1:length(taskids) % loop through tasks
        
        task = taskids{tidx};
        
        % read in nifti file
        nii_nam   = [func_path sublist{sub} '_ses-1_task-' task ...
                       '_space-MNI152NLin2009cAsym_desc-preproc_bold.nii'];
        niii_strc = spm_vol(nii_nam);
        nii_4D    = spm_read_vols(niii_strc);
        nii_4D    = nii_4D(:,:,:,start_TR:end_TR); 
        
        % read in (original) anatomical file
        anat_nam = [anat_path sublist{sub} ...
                       '_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii'];
        
        % ---------------- only do this once -----------------------------
         %    adjust dimensions of anat nifti so they fit functional
        if tidx == 1
            [bb,VX] = spm_get_bbox(niii_strc);
            resize_img(anat_name_mni, VX, bb);
        end
        % -----------------------------------------------------------------
        
        % read in ROI file
        roi_nam   = [roi_path 'rAll_ROIs.nii'];
        roi_strc  = spm_vol(roi_nam);
        roi_3d    = spm_read_vols(roi_strc);
        
        % get coordinates of PAC area (defined by roi-mask)
        PAC_vox = [];
        for i = 1:size(roi_3d,1)
            for j = 1:size(roi_3d,2)
                for k = 1:size(roi_3d,3)
                    if roi_3d(i,j,k) > 0
                        PAC_vox= [PAC_vox; [i,j,k]];
                    end
                end
            end
        end
        
        
        % read in (resized) anatomical file
        anat_nam  = [anat_path 'r' sublist{sub} ...
                      '_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii'];
        anat_strc = spm_vol(anat_nam);
        anat_3d   = spm_read_vols(anat_strc);
        
        
        % convert anat nifti to normalized values from 0 to 1
        anat_norm   = mat2gray(anat_3d);
        
        % read in model
        modelname  = [model_path '/' sublist{sub} '_ses-1_task-' ...
                      task 'ModelStruct.mat'];
        model      = load(modelname);
        
        %------------------------------------------------------------------
        % cut model so that only chosen timeframe is considered
        rows_to_remove = [1:(start_TR-1) (end_TR+1):size(model.Hz1000,1)];
        model_clean    = structfun(@(x) (removerows(x, 'ind',     ...
                                         rows_to_remove)), model, ...
                                         'UniformOutput', false);
        labels         = fieldnames(model);
        
        % create coeff structure for output
        coeff.(task) = zeros(length(PAC_vox),(length(fieldnames(model))+3));
        coeff.(task)(:,1:3) = PAC_vox;
        
         
        %------------------------------------------------------------------
        % now loop through PAC voxels to determine correlation coefficient
        for pac_idx = 1:length(PAC_vox)
            
            % read out indices for nifti image
            xdim = PAC_vox(pac_idx,1);
            ydim = PAC_vox(pac_idx,2);
            zdim = PAC_vox(pac_idx,3);
            
            % read out voxel timeseries per voxel
            x = [1:size(nii_4D,4)];
            xmat = vec2mat(x,1);
            y = squeeze(nii_4D(xdim,ydim,zdim,:));
            fit_y = smooth(xmat,y,0.99,'loess'); % fit timeseries points
            detrend_y = y-fit_y; %detrend timeseries
            
            % loop through different model types
            for label_idx = 1:length(labels)
                current_design = model_clean.(labels{label_idx});
                crosscorr = xcorr(current_design,detrend_y,0,'normalized');
                coeff.(task)(pac_idx,label_idx+3)=crosscorr;
            end
        end
        
        %------------------------------------------------------------------
        % create average matrix across tasks
        all_coeffs = coeff.(task);
        
        % find highest value per row
        all_idx = [];
        for row  = 1:length(all_coeffs)
            [val,idx] = max(all_coeffs(row,4:14));
            best_freq = [val, idx];
            all_idx = [all_idx;best_freq];
        end
        
        all_idx = [all_idx,PAC_vox]; % add respective voxel index again
      
         %-----------------------------------------------------------------
        namemat =  [out_path sublist{sub} task '_ContrastingTRs.mat']
        save(namemat,'all_coeffs_ordered')
        
    end
end