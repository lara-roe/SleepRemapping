
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define paths and filenames
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


datapath <- "derivatives/fmriprep/"
infopath <- "derivatives/sleepscoring/"

# read in sub info
sub_info <- read.csv2(paste0(infopath,"SleepScoring_UsableTRs_N2.csv"))
subs <- sub_info$Subject

wake_motion <- numeric()
sleep_motion <- numeric()

for (sub in subs){
  start_tr <- sub_info$Start_TR[sub_info$Subject==sub]
  end_tr   <- sub_info$End_TR[sub_info$Subject==sub]
  run_sleep <- sub_info$Run_Sleep[sub_info$Subject==sub]
  run_wake <- sub_info$Run_Wake[sub_info$Subject==sub]
  
  wake_file <- paste0("sub-S",ifelse(sub<10,"0",""),sub,"_ses-1_task-",run_wake,"_desc-confounds_regressors.tsv")
  wake_confounds <- read.table(paste0(datapath,wake_file), sep = '\t', header = TRUE)
  wake_confounds <- wake_confounds[start_tr:end_tr,]
  wake_confounds$framewise_displacement <- as.numeric(as.character(wake_confounds$framewise_displacement))
  fd_sub_wake <- mean(wake_confounds$framewise_displacement)
  
  sleep_file <- paste0("sub-S",ifelse(sub<10,"0",""),sub,"_ses-1_task-",run_sleep,"_desc-confounds_regressors.tsv")
  sleep_confounds <- read.table(paste0(datapath,sleep_file), sep = '\t', header = TRUE)
  sleep_confounds <- sleep_confounds[start_tr:end_tr,]
  sleep_confounds$framewise_displacement <- as.numeric(as.character(sleep_confounds$framewise_displacement))
  fd_sub_sleep <- mean(sleep_confounds$framewise_displacement)
  
  wake_motion <- c(wake_motion,fd_sub_wake)
  sleep_motion <- c(sleep_motion,fd_sub_sleep)
}


t.test(wake_motion,sleep_motion)
