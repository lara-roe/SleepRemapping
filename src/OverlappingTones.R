#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SleepRemapping Project
#
# Evaluate overlapping voxels across states
# (wake versus sleep)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


library(rmatio)
library(tidyverse)
library(gridExtra)
library(rgl)
library(corrplot)
library(scmamp)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define paths and filenames
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

datapath <- "/derivatives/maps/"
infopath <- "/derivatives/sleepscoring/"

sub_info <- read.csv2(paste0(infopath,"SleepScoring_UsableTRs_N2.csv"))
sub_info$Total_TR <- sub_info$End_TR - sub_info$Start_TR

total_tr_plot <- ggplot(sub_info, aes(x=Subject, y = Total_TR)) + geom_bar(stat="identity", fill = "steelblue") + 
                 ggtitle("Number of usable TRs per participant") + ylab("Total number of TRs per run") +  theme_modern() + 
                theme(plot.title = element_text(hjust = 0.5)) 


subs <- sub_info$Subject

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
# Read in files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
overlap_all <- numeric()

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
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Only look at highest coefficient per voxel
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  Wak[,"max"]  <- apply(Wak[,4:14],1,max) # read out max value per row
  Wak[,"name"] <- apply(Wak[,4:14],1,function(x) names(Wak[,4:14]) [which(x==max(x))]) # read out column name of the max val
  
  Slp[,"max"]  <- apply(Slp[,4:14],1,max) # read out max value per row
  Slp[,"name"] <- apply(Slp[,4:14],1,function(x) names(Slp[,4:14]) [which(x==max(x))]) # read out column name of the max val
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Throw out voxels with coefficients < 0.1
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  Wak_clean <- Wak[,c(1:3,15:16)]
  Slp_clean <- Slp[,c(1:3,15:16)]
  
  Wak_clean <- Wak_clean[Wak_clean$max>0.1,]
  Slp_clean <- Slp_clean[Slp_clean$max>0.1,]
  
  Wak_clean <- Wak_clean[,c(1:3,5)]
  Slp_clean <- Slp_clean[,c(1:3,5)]
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Find matching rows across states
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #transposte dataframe for match() function to work (is applied column wise)
  
  Wak_t <- transpose(Wak_clean)
  Slp_t <- transpose(Slp_clean)  
  
  matching_rows <- match(Wak_t,Slp_t)
  matching_rows <- matching_rows[!is.na(matching_rows)]
  
  overlap <- Slp_clean[matching_rows,]
  sub_overlap <- overlap %>% group_by(name) %>% summarize(n())
  
  # figure out if there's no overlap for some tones
  if (length(sub_overlap$name) != 11){
  missing_tones <- setdiff(Slp_clean$name, sub_overlap$name)
  miss_mat <- cbind(missing_tones,0)
  colnames(miss_mat) <- c("name","n()")
  miss_mat <- as.data.frame(miss_mat)
  sub_overlap <- rbind(sub_overlap,miss_mat)
  }
  
  sub_overlap <- sub_overlap[order(sub_overlap$name),]

  sub_overlap_t <- sub_overlap %>% as.data.frame() %>% t()
  overlap_all <- rbind(overlap_all, sub_overlap_t[2,])
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Group statistics
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

overlap_all <- as.data.frame(overlap_all)

# add subject id
colnames(overlap_all) <- sub_overlap$name
overlap_all$name <- subs

# only from 707 onwards

overlap_all_707 <- overlap_all %>% select(Hz707,Hz1000,Hz1414,Hz2000,Hz2828,Hz4000)
overlap_all_707 <- mutate_if(overlap_all_707, is.factor, ~ as.numeric(levels(.x))[.x])

#convert to long format
overlap_all_long <-  overlap_all %>% gather("Hz1000":"Hz707",key="tone",value="NrOfVox")

# function to compute p-mat
p.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      
      tmp <- t.test(mat[, i], mat[, j], paired = TRUE,method ="fdr",, ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      p.mat[i, j] <- p.adjust(pmat[i,j],method="fdr",n=15)
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat 
}

# matrix of t statistics
t.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  t.mat<- matrix(NA, n, n)
  diag(t.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      
      tmp <- t.test(mat[, i], mat[, j], paired = TRUE, method ="fdr",...)
      t.mat[i, j] <- t.mat[j, i] <- tmp$statistic
    }
  }
  colnames(t.mat) <- rownames(t.mat) <- colnames(mat)
  t.mat
}


pmat <- p.mtest(overlap_all_707)
tmat <- t.mtest(overlap_all_707)


plot1 <- corrplot(tmat,type="upper",  tl.col="black", tl.srt=45,
                  p.mat = pmat, sig.level = 0.05, is.corr = FALSE, title="Post hoc comparisons - FDR corrected",  mar=c(0,0,1,0))

plot2 <- corrplot(tmat,type="upper",  tl.col="black", tl.srt=45,
                  p.mat = pmat, insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "white", is.corr = FALSE)

# create numerical (then factor) frequency column of tones to get rid of annoying Hz string 
overlap_all_long$freq <- gsub("Hz","",overlap_all_long$tone)
overlap_all_long$freq  <- as.numeric(overlap_all_long$freq)
overlap_all_long$freq <- as.factor(overlap_all_long$freq)

overlap_all_long$NrOfVox <- as.numeric(overlap_all_long$NrOfVox)
overlap_all_long$name <- as.factor(overlap_all_long$name)
overlap_all_long$tone <- as.factor(overlap_all_long$tone)

overlap_from707 <- overlap_all_long[overlap_all_long$freq == 707 | overlap_all_long$freq == 1000 | overlap_all_long$freq == 1414 |
                                    overlap_all_long$freq == 2000 | overlap_all_long$freq == 2828 | overlap_all_long$freq == 4000 ,]

condcomp <- function(cond1,cond2,method="none") {
  p <- numeric()
  for (i in 1:ncol(cond1)) {
    p <- c(p,t.test(cond1[,i],cond2[,i],paired=TRUE)$p.value)
  }
  
  p <- p.adjust(p,method)
  return(p)
}



library(afex)
aov1 <- aov_ez(data=overlap_all_long, dv = "NrOfVox", id = "name", within =  c("tone"), factorize = FALSE) 
summary(aov1)

aov2 <- aov_ez(data=overlap_from707, dv = "NrOfVox", id = "name", within =  c("tone"), factorize = FALSE) 
summary(aov2)

library(see)
library("viridis") 

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


plot_overlap <-  ggplot(overlap_all_long, aes(x=freq,y=NrOfVox, fill= freq)) + geom_violindot(size_dots = 20) +
                  theme_modern() + ylab("Nr of overlapping voxels") + xlab("Tone in Hz") +
                  scale_fill_viridis_d() + theme(legend.position = "none")


raincloud_plot  <- ggplot(data = overlap_all_long, aes(y = NrOfVox, x = freq, fill = freq)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  xlab("Tone (in Hz)") + ylab("Number of voxels") + ggtitle("Overlap between sleep and wakefulness per tone") +
  geom_point(aes(y = NrOfVox, col = name), position = position_jitter(width = .15), size = 1, alpha = 0.8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  scale_color_viridis_d(option="B",end = 0.8) +
  scale_fill_viridis_d() + 
#  scale_fill_brewer(palette = "RdYlBu") +
  theme_classic() +
 theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

ggsave(filename=paste0("OverlapPerTone.jpg"),raincloud_plot,
       width=7,height=5)

  

