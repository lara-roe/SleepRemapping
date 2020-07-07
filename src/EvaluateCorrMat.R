library(rmatio)
library(tidyverse)
library(gridExtra)
library(rgl)
library(corrplot)
library(nlme)
library(lme4)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define paths and filenames
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


datapath <- "/derivatives/maps/"
infopath <- "/derivatives/sleepscoring/"
sub_info <- read.csv2(paste0(infopath,"SleepScoring_UsableTRs_N2.csv"))

subs <- sub_info$Subject
all_subs <- numeric()
all_subs_all <- numeric()

wak_mats <- sub_info$Run_Wake
slp_mats <- sub_info$Run_Sleep

wak_all <- numeric()
slp_all <- numeric()


for (sub in 1:length(subs)){
  sub_cod <- paste0("sub-S",ifelse(subs[sub]<10,"0",""),subs[sub])
  wak_mat <- paste0(sub_cod,wak_mats[sub],"_ContrastingTRs_N2.mat")
  slp_mat <- paste0(sub_cod,slp_mats[sub],"_ContrastingTRs_N2.mat")
  
  wak_all <- cbind(wak_all,wak_mat)
  slp_all <- cbind(slp_all,slp_mat)
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in general model 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mod_dir <- "/Volumes/Seagate/Experiments/SleepRemapping/Data/Models/"
mod_file <- list.files(mod_dir,full.names = TRUE)

mod_mat <- read.mat(mod_file[1])
mod_df  <- as.data.frame(mod_mat)
mod_cor <- cor(mod_df)

#png(filename = "model_correlations.png", width = 900/72*600, height = 450/72*600, res = 600, bg = "transparent")
#par(mfrow=c(1, 1))
#corrplot(mod_cor, method = "color", tl.col = "black", is.corr = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare data placeholders and labels
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_coeffs_all <- numeric()

# Create empty output matrices
# Group matrix
Wak_Cor_group <- matrix(0,11,11)
Slp_Cor_group <- matrix(0,11,11)

# Thresholded group matrix
Wak_thresh_Cor_group <- matrix(0,11,11)
Slp_thresh_Cor_group <- matrix(0,11,11)

# Empty df for mixed model analysis mixed model analysis
mixmod_df <- numeric()

# Prepare labels of tones and tone pairs
tones <-  c("Hz125","Hz177","Hz250","Hz354","Hz500","Hz707","Hz1000","Hz1414","Hz2000","Hz2828","Hz4000")
tonepairs <- expand.grid(tones, tones) %>% filter(Var1 != Var2) %>% apply(1, paste, collapse = "_") %>% .[1:55]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (i in 1:length(wak_all)){
  wakfile <- wak_all[i]
  slpfile <- slp_all[i]
  sub_cod <- subs[i]
  
  Wak <- read.mat(paste0(datapath,wakfile))
  Slp <- read.mat(paste0(datapath,slpfile))
  
  Wak <- as.data.frame(Wak$all_coeffs)
  Slp <- as.data.frame(Slp$all_coeffs)
  
  colnames(Wak) <- c("x","y","z","Hz125","Hz177","Hz250","Hz354","Hz500","Hz707","Hz1000","Hz1414","Hz2000","Hz2828","Hz4000")
  colnames(Slp) <- c("x","y","z","Hz125","Hz177","Hz250","Hz354","Hz500","Hz707","Hz1000","Hz1414","Hz2000","Hz2828","Hz4000")
  
  Wak <- Wak %>% select(-c("x","y","z"))
  Slp <- Slp %>% select(-c("x","y","z"))
  
  Wak_thresh <- Wak %>% filter_all(any_vars(.> 0.1))
  Slp_thresh <- Slp %>% filter_all(any_vars(.> 0.1))
  
  Wak_Cor <- cor(Wak)
  Slp_Cor <- cor(Slp)
  
  Wak_thresh_cor <- cor(Wak_thresh)
  Slp_thresh_cor <- cor(Slp_thresh)
  
  Wak_Cor_group <- Wak_Cor_group + Wak_Cor
  Slp_Cor_group <- Slp_Cor_group + Slp_Cor
  
  Wak_thresh_Cor_group <- Wak_thresh_Cor_group + Wak_thresh_cor
  Slp_thresh_Cor_group <- Slp_thresh_Cor_group + Slp_thresh_cor
    
  # prepare df for later mixed model analysis
   Wak_sub <- Wak_thresh_cor[lower.tri(Wak_thresh_cor)] %>% 
    unlist() %>% 
    bind_cols(
      id = rep(sub_cod, each = (11*12/2-11)), 
      tonepair = tonepairs,
      vigilance = rep("wake", (11*12/2-11)),
      value = .
    )
   
   Slp_sub <- Slp_thresh_cor[lower.tri(Slp_thresh_cor)] %>% 
     unlist() %>% 
     bind_cols(
       id = rep(sub_cod, each = (11*12/2-11)), 
       tonepair = tonepairs,
       vigilance = rep("sleep", (11*12/2-11)),
       value = .
     )
   
   sub_df <- rbind(Wak_sub,Slp_sub)
   mixmod_df <- rbind(mixmod_df,sub_df)
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Correct output dataframe
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Divide by 17 (coefficients were simply added above)
Wak_Cor_group <- Wak_Cor_group/17
Slp_Cor_group <- Slp_Cor_group/17

# Divide by 17 (coefficients were simply added above)
Wak_thresh_Cor_group <- Wak_thresh_Cor_group/17
Slp_thresh_Cor_group <- Slp_thresh_Cor_group/17

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Plot data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# without model correlation subtraction
png(filename = "simmaps_perVox.png", width = 900/72*600, height = 450/72*600, res = 600, bg = "transparent")
par(mfrow=c(1, 2))
corrplot(Wak_Cor_group, method = "color", tl.col = "black", is.corr = FALSE)
corrplot(Slp_Cor_group, method = "color", tl.col = "black", is.corr = FALSE)
dev.off()

# with model correlation subtraction
png(filename = "simmaps_modelsubtracted.png", width = 900/72*600, height = 450/72*600, res = 600, bg = "transparent")
par(mfrow=c(1, 2))
corrplot(Wak_thresh_Cor_group, method = "color", tl.col = "black", is.corr = FALSE,cl.lim = c(-0.37,0.4))
corrplot(Slp_thresh_Cor_group, method = "color", tl.col = "black", is.corr = FALSE, cl.lim = c(-0.37,0.4))
dev.off()

# difference between wake and sleep
png(filename = "diffmap_wakevssleep.png", width = 900/72*600, height = 450/72*600, res = 600, bg = "transparent")
par(mfrow=c(1, 1))
corrplot(Wak_thresh_Cor_group-Slp_thresh_Cor_group, method = "color", tl.col = "black", is.corr = FALSE)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Mixed model analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot the distributions
mixmod_df %>% ggplot(aes(x = value, fill = as_factor(vigilance))) + geom_density()

# plot distribution of log-transform
mixmod_df$log_val <- log(mixmod_df$value+1)

mixmod_df %>% ggplot(aes(x = log_val, fill = as_factor(vigilance))) + geom_density()

# perform different models
m1 <- lmer(value ~  tonepair*vigilance + (1 | id), data = mixmod_df)
plot(resid(m1))

m2 <- lme(value ~ tonepair*vigilance, mixmod_df,method = 'ML',
          ~1|id,weights = varIdent(form=~1|tonepair*vigilance),
          control=lmeControl(maxIter = 50, opt=("optim")))
plot(resid(m2))

m3 <- lme(log_val ~ tonepair*vigilance, mixmod_df,method = 'ML',
          ~1|id,weights = varIdent(form=~1|tonepair*vigilance),
          control=lmeControl(maxIter = 50, opt=("optim")))
plot(resid(m3))