clear all;
close all;

%% Prepare all paths and specify which runs should be analyzed

addpath(genpath('/home/lr499/rds/rds-tb419-bekinschtein/Lara/Programs/toolbox/'));

reg_dir = '/home/lr499/rds/rds-tb419-bekinschtein/Lara/SleepRemapping/derivatives/regressors/PhaseEncoding/';
fmridir =  '/home/lr499/rds/rds-tb419-bekinschtein/Lara/SleepRemapping/derivatives/fmriprep/'; % root inputdir for sublist
outpath = '/home/lr499/rds/rds-tb419-bekinschtein/Lara/Scripts/phaseEncoding/Models/';

sublist=dir(fullfile(fmridir,'sub*'));
isFile   = [sublist.isdir];
sublist = {sublist(isFile).name};

% define number of sessions and task
scode = 'ses-1';
taskids = {'WakeAscending'};

% choose subject (only necessary if you want to underlie voxel timeseries)
for sub = [2]
    
    figHandle = figure;
    
    model_path = '/home/lr499/rds/rds-tb419-bekinschtein/Lara/Scripts/phaseEncoding/Models/';
    func_path  = ['/home/lr499/rds/rds-tb419-bekinschtein/Lara/SleepRemapping/derivatives/fmriprep/' sublist{sub} '/ses-1/func/'];
    anat_path  = ['/home/lr499/rds/rds-tb419-bekinschtein/Lara/SleepRemapping/derivatives/fmriprep/' sublist{sub} '/anat/'];
    roi_path   = '/home/lr499/rds/rds-tb419-bekinschtein/Lara/SleepRemapping/derivatives/spm_2ndlevel/ROIs/anatlabels_surf_mni/mni152_posterior2anterior/';
    out_path   = '/home/lr499/rds/rds-tb419-bekinschtein/Lara/SleepRemapping/derivatives/phaseencoding/plots/';
    
    % loop through tasks
    for taskidx = 1:length(taskids)
        task = taskids{taskidx};
        
        % read in nifti file
        nifti_name = [func_path sublist{sub} '_ses-1_task-' task '_space-MNI152NLin2009cAsym_desc-preproc_bold.nii'];
        nifti_struct =  spm_vol(nifti_name);
        nifti_4D = spm_read_vols(nifti_struct);
        %nifti_4D = nifti_4D(:,:,:,8:size(nifti_4D,4));
        
        %define start and end TR
        start_TR=1;
        end_TR=417;
        
        % read in model
        modelname = [model_path sublist{sub} '_ses-1_task-' task 'ModelStruct.mat'];
        model = load(modelname);
        rows_to_remove = [(end_TR+1):size(model.Hz1000,1)];
        model_clean = structfun(@(x) (removerows(x, 'ind', rows_to_remove)), model, 'UniformOutput', false);
        labels = fieldnames(model);
        
        % read in initial design file to check timing of tones
        events_csv = [reg_dir sublist{sub},'_',scode,'_task-',task, '_reg.csv'];
        design = readcsv(events_csv,';',',',NaN);
        
        % read in ROI file
        roi_name = [roi_path 'rAll_ROIs.nii'];
        roi_struct = spm_vol(roi_name);
        roi_3d = spm_read_vols(roi_struct);
        
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
        
        % read out indices for nifti image
        xdim = PAC_vox(295,1);   % pick voxel here, arbitrary choice of ten
        ydim = PAC_vox(295,2);
        zdim = PAC_vox(295,3);
        
        % create figure holder
        figure1=figure('Position', [50, 50, 860, 1300]);
        
        % create colormap for plotting
        colormap = jet(11);
        colormap = flipud(colormap); % flipped so first row is red
        
        col_idx = 0; % start color count at 0 to later index colormap
        for l_idx = [2,4,6,8,10,1,3,5,7,9] % ordered here so tones will be presented in ascending order
            
            col_idx  = col_idx + 1; % go one up in color vector
            rgb_val = colormap(col_idx,:);
            
            x = [1:size(nifti_4D,4)];
            xmat = vec2mat(x,1);set(gca,'YTick',[]); %which will get rid of all the markings for the y axis
            
            y = squeeze(nifti_4D(xdim,ydim,zdim,:));
            fit_y = smooth(xmat,y,0.99,'loess');
            detrend_y = y-fit_y;
            detrend_y_norm = normalize(detrend_y);
            subplot(11,1,col_idx)
            
            % ----------------------------------------------------------
            % add this if you want to see voxel timeseries underneath
            
            %plot(x, detrend_y_norm,'Color', [0 0 0]);
            %hold on
            
            % ----------------------------------------------------------
            plot(x, model_clean.(labels{l_idx}),'Color', rgb_val)
            xticks([5:5:410])
            set(gca,'Xticklabel',[])
            
            
            % ----------------------------------------------------------
            % add  if you want model rescaled for better visualization
            
            %rescaled_model = rescale(model_clean.(labels{l_idx}),-2,2);
            % plot(x, rescaled_model,'Color', rgb_val)
            
            % ----------------------------------------------------------
            
            title(labels{l_idx}) 
        end
        
        name = [out_path sublist{sub} task 'BestModelFit_withTicks.png'];
        print(figure1,'-dpng','-r0',name)
    end
end
