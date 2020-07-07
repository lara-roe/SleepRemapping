#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SleepRemapping Project
#
# Evaluate whether correlation coefficients differ per state
# (wake versus sleep)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


library(rmatio)
library(tidyverse)
library(gridExtra)
library(rgl)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define paths and filenames
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_coeffs_all <- numeric()

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
  
  # from wide to long
  
  Wak_long <- Wak %>% gather(Tone,Coeff,Hz125:Hz4000)
  Slp_long <- Slp %>% gather(Tone,Coeff,Hz125:Hz4000)
  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Only look at highest coefficient per voxel
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  Wak[,"max"]  <- apply(Wak[,4:14],1,max) # read out max value per row
  Wak[,"name"] <- apply(Wak[,4:14],1,function(x) names(Wak[,4:14]) [which(x==max(x))]) # read out column name of the max val
  
  Slp[,"max"]  <- apply(Slp[,4:14],1,max) # read out max value per row
  Slp[,"name"] <- apply(Slp[,4:14],1,function(x) names(Slp[,4:14]) [which(x==max(x))]) # read out column name of the max val
  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Throw out voxels with coefficients < 0.1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Wak_clean <- Wak[,c(1:3,15:16)]
  Slp_clean <- Slp[,c(1:3,15:16)]
  
  Wak_clean <- Wak_clean[Wak_clean$max>0.1,]
  Slp_clean <- Slp_clean[Slp_clean$max>0.1,]
  
  nrow(Wak_clean)
  nrow(Slp_clean)
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Relabel and group tones
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# summarize per tone 

Wake_clean_all <- Wak_clean %>% select(max,name)
Wake_clean_all <- Wake_clean_all %>% group_by(name) %>% dplyr::summarize(mean(max))
Wake_clean_all$sub <- rep(sub_cod,nrow(Wake_clean_all))
Wake_clean_all$drowsy <- rep("wake",nrow(Wake_clean_all))

Sleep_clean_all <- Slp_clean %>% select(max,name)
Sleep_clean_all <- Sleep_clean_all %>% group_by(name) %>% dplyr::summarize(mean(max))
Sleep_clean_all$sub <- rep(sub_cod,nrow(Sleep_clean_all))
Sleep_clean_all$drowsy <- rep("sleep",nrow(Sleep_clean_all))

sub_coeffs_all <- rbind(Wake_clean_all,Sleep_clean_all)
all_coeffs_all <- rbind(all_coeffs_all,sub_coeffs_all)

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Prepare df for anova
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
all_coeffs_all <- as.data.frame(all_coeffs_all)
colnames(all_coeffs_all) <- c("tone","coeff","sub","drowsy")
all_coeffs_all$drowsy <- as.factor(all_coeffs_all$drowsy)
all_coeffs_all$sub <- as.factor(all_coeffs_all$sub)
all_coeffs_all$tone <- as.factor(all_coeffs_all$tone)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculat anova
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(afex)

sleep_anova_all <- aov_ez(data=all_coeffs_all, dv = "coeff", id = "sub", within =  c("tone","drowsy")) 
summary(sleep_anova_all)

#Post-hoc tests
ls1 <- lsmeans(sleep_anova_all, c("drowsy"),by="tone")
diffs <- update(pairs(ls1), by=NULL, adjust = "none")
write.table(diffs, file = "diffs_thresh02.txt", sep = ";", quote = FALSE, row.names = F)

# create numerical (then factor) frequency column of tones to get rid of annoying Hz string 
all_coeffs_all$freq <- gsub("Hz","",all_coeffs_all$tone)
all_coeffs_all$freq  <- as.numeric(all_coeffs_all$freq)
all_coeffs_all$freq <- as.factor(all_coeffs_all$freq)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# =========================== function definition ===========================

## define function of geom_flat_volin
# devtools::install_github(repo = "IndrajeetPatil/ggstatsplot")

"%||%" <- function(a, b) {
  if (!is.null(a))
    a
  else
    b
}

geom_flat_violin <-
  function(mapping = NULL,
           data = NULL,
           stat = "ydensity",
           position = "dodge",
           trim = TRUE,
           scale = "area",
           show.legend = NA,
           inherit.aes = TRUE,
           ...) {
    ggplot2::layer(
      data = data,
      mapping = mapping,
      stat = stat,
      geom = GeomFlatViolin,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(trim = trim,
                    scale = scale,
                    ...)
    )
  }

GeomFlatViolin <-
  ggproto(
    "GeomFlatViolin",
    Geom,
    setup_data = function(data, params) {
      data$width <- data$width %||%
        params$width %||% (resolution(data$x, FALSE) * 0.9)
      
      # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
      data %>%
        dplyr::group_by(.data = ., group) %>%
        dplyr::mutate(
          .data = .,
          ymin = min(y),
          ymax = max(y),
          xmin = x,
          xmax = x + width / 2
        )
    },
    
    draw_group = function(data, panel_scales, coord)
    {
      # Find the points for the line to go all the way around
      data <- base::transform(data,
                              xminv = x,
                              xmaxv = x + violinwidth * (xmax - x))
      
      # Make sure it's sorted properly to draw the outline
      newdata <-
        base::rbind(
          dplyr::arrange(.data = base::transform(data, x = xminv), y),
          dplyr::arrange(.data = base::transform(data, x = xmaxv), -y)
        )
      
      # Close the polygon: set first and last point the same
      # Needed for coord_polar and such
      newdata <- rbind(newdata, newdata[1,])
      
      ggplot2:::ggname("geom_flat_violin",
                       GeomPolygon$draw_panel(newdata, panel_scales, coord))
    },
    
    draw_key = draw_key_polygon,
    
    default_aes = ggplot2::aes(
      weight = 1,
      colour = "grey20",
      fill = "white",
      size = 0.5,
      alpha = NA,
      linetype = "solid"
    ),
    
    required_aes = c("x", "y")
  )
# =========================== function definition ===========================


raincloud_plot  <- ggplot(data = all_coeffs_all, aes(y = coeff, x = freq, fill = drowsy)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  xlab("Tone (in Hz)") + ylab("Strength of correlation coefficient") +
  geom_point(aes(y = coeff, col = sub), position = position_jitter(width = .15), size = 1, alpha = 0.8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  scale_color_viridis_d(option="B",end = 0.8) +
  scale_fill_viridis_d(begin=0.5) + 
  theme_bw() +
  raincloud_theme +
  guides(fill=guide_legend(title="Wakefulness"),color=guide_legend(title="Participant ID"))

ggsave(filename=paste0("CoefficientsPerTone.jpg"),raincloud_plot,
       width=11,height=7)

plot <- 
  ggplot(all_coeffs,aes(x=tone,y=coeff)) + geom_boxplot(aes(fill=drowsy))+
  facet_grid(.~ sub) 
ggsave("CorrCoeffs_ControlWake.jpg",height=8,width=8,units="in")


               