#Quantifying Scleractinian Coral Growth Forms Code

# R version 3.5.0 (2018-04-23)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows >= 8 x64 (build 9200)
#
# Matrix products: default
#
# locale:
#   [1] LC_COLLATE=English_United Kingdom.1252  LC_CTYPE=English_United Kingdom.1252    LC_MONETARY=English_United Kingdom.1252 LC_NUMERIC=C
# [5] LC_TIME=English_United Kingdom.1252
#
# attached base packages:
#   [1] grid      parallel  stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
# [1] e1071_1.7-1       GGally_1.4.0      rgl_0.100.19      magick_2.0        pca3d_0.10        nnet_7.3-12       data.table_1.12.2 gtable_0.3.0
# [9] gridExtra_2.3     png_0.1-7         ggthemes_4.1.1    ade4_1.7-13       ellipse_0.4.1     reshape2_1.4.3    factoextra_1.0.5  lme4_1.1-21
# [17] Matrix_1.2-14     caret_6.0-84      lattice_0.20-35   multcomp_1.4-10   TH.data_1.0-10    MASS_7.3-49       survival_2.41-3   mvtnorm_1.0-10
# [25] ggrepel_0.8.0     ggpubr_0.2        magrittr_1.5      forcats_0.4.0     stringr_1.4.0     dplyr_0.8.0.1     purrr_0.3.2       readr_1.3.1
# [33] tidyr_0.8.3       tibble_2.1.1      ggplot2_3.1.1     tidyverse_1.2.1
#
# loaded via a namespace (and not attached):
# [1] nlme_3.1-137            lubridate_1.7.4         progress_1.2.0          RColorBrewer_1.1-2      webshot_0.5.1           httr_1.4.0
# [7] tools_3.5.0             backports_1.1.4         R6_2.4.0                rpart_4.1-13            lazyeval_0.2.2          colorspace_1.4-1
# [13] manipulateWidget_0.10.0 withr_2.1.2             prettyunits_1.0.2       tidyselect_0.2.5        compiler_3.5.0          cli_1.1.0
# [19] rvest_0.3.3             xml2_1.2.0              sandwich_2.5-1          labeling_0.3            scales_1.0.0            digest_0.6.18
# [25] minqa_1.2.4             pkgconfig_2.0.2         htmltools_0.3.6         htmlwidgets_1.3         rlang_0.3.4             readxl_1.3.1
# [31] rstudioapi_0.10         shiny_1.3.2             generics_0.0.2          zoo_1.8-5               jsonlite_1.6            crosstalk_1.0.0
# [37] ModelMetrics_1.2.2      Rcpp_1.0.1              munsell_0.5.0           stringi_1.4.3           yaml_2.2.0              plyr_1.8.4
# [43] recipes_0.1.5           promises_1.0.1          crayon_1.3.4            miniUI_0.1.1.1          haven_2.1.0             splines_3.5.0
# [49] hms_0.4.2               knitr_1.22              pillar_1.3.1            boot_1.3-20             codetools_0.2-15        stats4_3.5.0
# [55] glue_1.3.1              modelr_0.1.4            nloptr_1.2.1            httpuv_1.5.1            foreach_1.4.4           cellranger_1.1.0
# [61] reshape_0.8.8           assertthat_0.2.1        xfun_0.6                gower_0.2.0             mime_0.6                prodlim_2018.04.18
# [67] xtable_1.8-4            broom_0.5.2             later_0.8.0             class_7.3-14            timeDate_3043.102       iterators_1.0.10
# [73] lava_1.6.5              ipred_0.9-9

#Set up ----

rm(list=ls())

outputdir <- ("figures")

#.....Libraries ----

library(tidyverse)
library("ggpubr")
library(ggrepel)
library("multcomp")
library("caret")
library("ellipse")
library("grid")
library("ggthemes")
library("png")
library("gridExtra")
library("gtable")
library("nnet")
library(GGally)

#.....Functions ----

#ggGetLegend:====
#extract a figure legend from a ggplot object.

ggGetLegend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#McFaddens psuedo-R2 ====
#calculates McFaddens psuedo-R2 for estimating goodness of fit for multinomial models

McFR2 <- function(model,nullmodel) { return(1-(logLik(model)/logLik(nullmodel)))}

#GGally custom correlation function ====
my_custom_cor <- function(data, mapping, color = I("black"), sizeRange = c(1, 5), ...) {

  # get the x and y data to use the other code
  x <- eval(mapping$x, data)
  y <- eval(mapping$y, data)
  ct <- cor.test(x,y)
  sig <- symnum(
    ct$p.value, corr = FALSE, na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("***", "**", "*", ".", " ")
  )
  r <- unname(ct$estimate)
  rt <- format(r, digits=2)[1]
  # since we can't print it to get the strsize, just use the max size range
  cex <- max(sizeRange)
  # helper function to calculate a useable size
  percent_of_range <- function(percent, range) {
    percent * diff(range) + min(range, na.rm = TRUE)
  }
  # plot the cor value
  ggally_text(label = as.character(rt),
              mapping = aes(),
              xP = 0.5, yP = 0.5,
              size = I(percent_of_range(cex * abs(r), sizeRange)),
              color = color,
              ...) +
    # add the sig stars
    geom_text(aes_string(x = 0.8, y = 0.8),
              label = sig,
              size = I(cex),
              color = color,
              ...) +
    # remove all the background stuff and wrap it with a dashed line
    theme_classic() +
    theme(panel.background = element_rect(colour = "black",linetype = "solid"),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank()
    )
}

#GGally custom smooth function ====
my_custom_smooth <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(color = "black", shape = 21, fill = "grey", alpha = 0.8) +
    geom_smooth(method = "loess", color = I("black"), se = F, alpha = 0.6, ...) +
    theme_classic() +
    theme(
      panel.background = element_rect(
        color = "black",
        linetype = "solid"
      ),
      axis.line = element_blank()
      #axis.ticks = element_blank(),
      #axis.text.y = element_blank(),
      #axis.text.x = element_blank()
    )

}
#GGally custom density function ====
my_custom_density <- function(data, mapping, ...) {
  ggplot(data = data, mapping=mapping) +
    geom_histogram(bins = 10, fill = "white",colour = "black") +
    #scale_x_continuous(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    #expand_limits(y=0) +
    theme_classic() +
    theme(
      panel.background = element_rect(
        color = "black",
        linetype = "solid"
      ),
      axis.line = element_blank()
      #axis.ticks = element_blank(),
      #axis.text.y = element_blank(),
      #axis.text.x = element_blank()
    )
}

#Zuur Functions ====
#source(paste0(funcdir,"ZuurHighstatLibV10.R"))

#.....ggplot settings ----

base_size <- 32

PlotTheme <- theme_foundation(base_size=base_size) +
  theme(plot.title=element_text(size=rel(1.2), face="bold", hjust = 0.5),
        text = element_text(),
        rect = element_rect(colour = "black", linetype = "solid"),
        panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.border = element_rect(colour = NA),
        axis.title = element_text(size = rel(1)),
        axis.title.y = element_text(angle=90,vjust =2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text = element_text(),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(),
        panel.grid.major = element_line(colour="#f0f0f0"),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(colour = NA),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.size= unit(0.5, "cm"),
        legend.key.height=unit(1.5,"cm"),
        legend.key.width=unit(1.5,"cm"),
        legend.spacing = unit(1, "cm"),
        legend.title = element_text(),
        legend.title.align = 0.5,
        #plot.margin=unit(c(10,5,5,5),"mm"),
        panel.spacing = unit(5, "mm"),
        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
        strip.text = element_text(face="bold"))

GFNinePallete <- c("#5DA5DA",
                   "#FAA43A",
                   "#B276B2",
                   "#F15854",
                   "#DECF3F",
                   "#60BD68",
                   "#F17CB0",
                   "#B2912F",
                   "#4D4D4D"
                   ) #pallete for nine growth forms

#.....Data ====

# Download data to 'data' folder from https://figshare.com/articles/CoralColonyMorphologyData_WholeColonies_MediumToHighQuality_LaserScanned/8115674
# Change code here to point to data
# DOI: 10.6084/m9.figshare.2067414.v1

AllDat <- read.csv("data/3DLaserScannedColonies_MedToHighQuality_updatedFractalDimension.csv")


#..........Extra variable calculation ====

AllDat$SizeClass <- factor(round(log10(AllDat$ColonyVol),0))

#..........data subsetting/sorting ====

#removed as too damaged
AllDat <- AllDat[-which(AllDat$ScanID == "MTQ62MeLhem"),]

#removed branching-closed due to odd assignment of categories/too few replicates
AllDat <- AllDat[-which(AllDat$GrowthForm == "BranchingClosed"),]

#remove growth forms with less than six replicates.
GFTotals <- table(AllDat$GrowthForm) #N per growth form category
LZDat <- AllDat[AllDat$GrowthForm %in% names(GFTotals[GFTotals >=6]),] #if less than 5, drop.
LZDat$GrowthForm <- droplevels(LZDat$GrowthForm) #drop empty levels.

#order data by decending colony volume in each growth form
GFs <- levels(LZDat$GrowthForm) #Existing growth form levels.
LZDat <- do.call("rbind", #bind list of data frames
                 lapply(1:length(GFs), function(x) #for each growth form
                   subset(LZDat[order(LZDat$ColonyVol), ], #subset the dataframe, ordered by colony volume
                          LZDat$GrowthForm[order(LZDat$ColonyVol)] %in% GFs[x]) #by growth form, ordered by colony volume.
                   ))

LZDat <- LZDat %>%
  arrange(GrowthForm, ColonyVol)


#Summarise data per growth form
GFSummary <- LZDat %>%
  group_by(GrowthForm) %>%
  summarise(n = n(),
            MinSize = range(ColonyVol)[1],
            MaxSize = range(ColonyVol)[2],
            SpeciesNumber = length(unique(Species)),
            SpeciesList = paste(unique(Species), collapse = ", "))

#Analysis ====
#.....PCA Analysis ====
#..........set-up ====

PCADat <- data.frame('Sphericity' = log(LZDat$Sphericity/(1-LZDat$Sphericity)), #logit Sphericity
                     'Packing'= log10(LZDat$Packing), #log10 packing
                     'Convexity' = log(LZDat$Convexity/(1-LZDat$Convexity)), #logit convexity
                     "Vertical2ndMomentAreaScaled" = log10(LZDat$ColonySA2ndMomentVertScaled), #log10 Varea
                     "Vertical2ndMomentVolumeScaled" = log10(LZDat$ColonyVol2ndMomentVertScaled), #log10 Vvol
                     'FractalDimension' = (LZDat$FractalDimension), #fractal dimension
                     'Volume' = log10(LZDat$ColonyVol), #log10 volume
                     "GrowthForm" = LZDat$GrowthForm #growth form
                     )


#..........analysis ====

PCA2 <- prcomp(PCADat[,c("Sphericity","Packing","Convexity","FractalDimension",
                         "Vertical2ndMomentVolumeScaled","Vertical2ndMomentAreaScaled")], center = T, scale. = T)

PCAggregate <- aggregate(PCA2$x[,1]~LZDat$GrowthForm, FUN = mean) #average PC1 score for sorting

#..........plotting ====

#sort LZDat and PCAdat by average PC1 Score
LZDat$GrowthForm <- factor(LZDat$GrowthForm,
                           levels = levels(LZDat$GrowthForm)[order(PCAggregate$`PCA2$x[, 1]`)])
PCADat$GrowthForm <- factor(PCADat$GrowthForm, levels = levels(LZDat$GrowthForm))


Loadings <- as.data.frame(cor(PCADat[,c("Sphericity","Packing","Convexity","FractalDimension",
                                    "Vertical2ndMomentVolumeScaled","Vertical2ndMomentAreaScaled")], PCA2$x)) #factor loadings
Arrows <- data.frame(xS = rep(0,nrow(Loadings)), #loading arrow x center
                     yS = rep(0,nrow(Loadings)), #loading arrow y center
                     xE = Loadings$PC1, #loading arrow x coordinate
                     yE = Loadings$PC2) #loading arrow y coordinate
Arrows <- Arrows *4 #scaling of arrows

Text <- data.frame(X = Arrows$xE *1.1, # x coordinate for variable labels, scaled away from arrow end
                  Y = Arrows$yE *1.1) # y coordinate for variable lables, scale away from arrow end
row.names(Text) <- c("S","P","C","FD","V[VOL]","V[AREA]") #variable labels
Vars <- c("S","P","C","FD","V[VOL]","V[AREA]") #variable labels

#Plotting coral images - not available
#Corals <- list.files(paste0(getwd(),"/data/","TraitExtremesSpecimens"), full.names = T)
#Max <- Corals[which(grepl(Corals, pattern = "Max") == T)] #variable maximum specimens
#Min <- Corals[which(grepl(Corals, pattern = "Min") == T)] #variable minumim specimens
#MaxPng <- lapply(1:length(Max), function(c) rasterGrob(readPNG(Max[c]),interpolate=TRUE)) #load max pngs into R
#MinPng <- lapply(1:length(Min), function(c) rasterGrob(readPNG(Min[c]),interpolate=TRUE)) #load minpngs in R
#PngXCoords <- c(Text$X[2],Text$X[4],Text$X[1],Text$X[6],Text$X[5],Text$X[3]) #coral image x coords
#PngYCoords <- c(Text$Y[2],Text$Y[4],Text$Y[1],Text$Y[6],Text$Y[5],Text$Y[3]) #coral image y coords
#DistanceFactor <- 1.25 #coralimage distance factor multiplier
#PngXCoords <- PngXCoords * DistanceFactor #scale image x coords
#PngYCoords <- PngYCoords * DistanceFactor #scale image y coords
                       #S,P,C,FD,2V,2A
#ScalingFactor <- 0.7 #image size scaling distance
#PngXMin <- PngXCoords - ScalingFactor #image min x coords
#PngXMax <- PngXCoords + ScalingFactor #image max x coords
#PngYMin <- PngYCoords - ScalingFactor #image min y coords
#PngYMax <- PngYCoords + ScalingFactor #image max y coords


PlotDat <- data.frame(PCA2$x,
                      GF = factor(PCADat$GrowthForm,
                                  levels = (c("Massive",
                                              "Submassive",
                                              "Digitate",
                                              "Arborescent",
                                              "Laminar",
                                              "Tabular",
                                              "Corymbose")))) #PCA data and growth form dataframe

#Biplot

Biplot <- ggplot(data = PlotDat, aes(x =PC1,y=PC2,colour= GF,fill = GF)) + #ggplot set up
  PlotTheme + #plotting theme
  stat_conf_ellipse(alpha = 0.7, #95% confidence ellipse around the mean
                    linetype = "dashed",
                    geom = "polygon",
                    level= 0.95) +
  stat_conf_ellipse(colour = "black", #95% confidence ellipse outline
                    alpha = 0.0,
                    linetype = "dashed",
                    geom = "polygon",
                    level= 0.95,
                    show.legend = F) +
  geom_segment(data = Arrows, #Outwards loading arrows
               aes(x= xS,y=yS,xend=xE,yend=yE),
               colour = "grey21",
               size = 0.85,
               arrow = arrow(length = unit(0.25,"cm"),
                             type = "closed"),
               inherit.aes = F) +
  geom_segment(data = Arrows, #Inwards loading arrows
               aes(x= xS,y=yS,xend=-xE,yend=-yE),
               linetype = "dashed",
               colour = "grey21",
               size = 0.85,
               inherit.aes = F) +
  geom_point(shape = 21, colour = "black" , size = 5,alpha = 0.7,show.legend = F) + #Scatter plot of observed PC1 and 2 values
  geom_text(data = as.data.frame(Vars), #Loading variable text labels
            aes(x = Text$X, y = Text$Y, label = Vars),
            size = 9,
            colour = "black" ,
            inherit.aes = F,
            parse = T) +
  #geom_text(aes(label = Labels), colour = "black") +
  scale_fill_manual(values= GFNinePallete) + #set pallete colours
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey31") + #horizontal reference line
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey31") + #verticalreference line
  #lapply(1:length(MaxPng), function(p) annotation_custom(MaxPng[[p]],PngXMin[p],PngXMax[p],PngYMin[p],PngYMax[p])) + #add maximum variable specimens
  #lapply(1:length(MinPng), function(p) annotation_custom(MinPng[[p]],-PngXMin[p],-PngXMax[p],-PngYMin[p],-PngYMax[p])) + #add minimum variable speciemens
  scale_color_manual(values = GFNinePallete) + #set pallete colours
  guides(colour = F, fill = guide_legend(nrow = 1)) + #remove colour guide, force fill legend to  single row
  theme(legend.title = element_blank(), #remove title
        legend.background = element_blank(), #remove background
        legend.justification = "center", #justify legend to center
        legend.position = c(0.5, 0.05), #legend position adjustment
        axis.ticks = element_blank(), #remove axis ticks
        panel.grid.major = element_blank()) + #remove panel grid lines
  coord_cartesian(xlim = c(-6,6), ylim = c(-5.5,4)) + #setplot x and y limits
  labs(fill = "Growth Form:")

#save as png
png(paste0(outputdir,"/BiplotFinal.png"), 1920, 1080, type='cairo')
Biplot
dev.off()

#save as pdf
ggsave(file=paste0(outputdir,"/BiplotFinal.pdf"), plot=Biplot, width=3.45, height=2, scale = 8)



#.....All shape + size variable correlations ====
#..........Set up ====
CorDat <- data.frame(
                     V = LZDat$ColonyVol/10000,
                     SA = LZDat$ColonySA/10000,
                     S = LZDat$Sphericity,
                     C = LZDat$Convexity,
                     P = LZDat$Packing,
                     FD = LZDat$FractalDimension,
                     VVOL = LZDat$ColonyVol2ndMomentVertScaled,
                     VAREA = LZDat$ColonySA2ndMomentVertScaled
                     )
CorLabels <- c("V[m^3]",
               "SA[m^2]",
               "S","C","P","FD",
               "V[VOL~mm^4]",
               "V[AREA~mm^3]")
UpperTheme <- PlotTheme

#......... Plot
PairPlot <- ggpairs(data = CorDat,
                    lower = list(continuous = wrap(my_custom_smooth, alpha = 0.7), PlotTheme),
                    upper = list(continuous = wrap(my_custom_cor, sizeRange = c(2,4)), UpperTheme),
                    diag = list(continuous = wrap(my_custom_density, labels = CorLabels)),
                    labeller = "label_parsed",
                    columnLabels = CorLabels) +
  theme(panel.spacing = unit(0, "mm"),
        strip.background = element_rect(fill = NULL, colour = "black"),
        strip.text = element_text(size = 22))
ggsave(file=paste0(outputdir,"/PairPlot.png"), plot=PairPlot, width=2, height=2,scale = 6)
ggsave(file=paste0(outputdir,"/PairPlot.pdf"), plot=PairPlot, width=2, height=2,scale = 6)


#.....Size by shape variables ----
#..........analysis ====

NewData <- lapply(1:length(levels(PCADat$GrowthForm)), function(g) { #for each growth form
  dat <- PCADat[which(PCADat$GrowthForm == levels(PCADat$GrowthForm)[g]),] #subset data by growth form
  data.frame(Volume = seq(from = range(dat$Volume)[1],
                          to = range(dat$Volume)[2],
                          by = 0.1), #create vector covering size range of growth form
             GrowthForm = rep(levels(dat$GrowthForm)[g],
                                 times = length(seq(from = range(dat$Volume)[1],
                                                    to = range(dat$Volume)[2],
                                                    by = 0.1)))) #create a vector of growth form category
})

NewData <- rbindlist(NewData) #create single dataframe from list

Vars <- c("Sphericity","Convexity", "Packing", "FractalDimension",  "Vertical2ndMomentVolumeScaled","Vertical2ndMomentAreaScaled") #selected shape variables

Mods <- lapply(1:length(Vars), function(x) { #for each variable
  Mod <- lm(get(Vars[x]) ~ Volume * GrowthForm, data = PCADat) #run a linear regression
  return(Mod)
})

GFLevels <- length(levels(LZDat$GrowthForm)) #length of the number of levels in the growth form data

CoefsList <- lapply(1:length(Vars), function(x) { #for each variables
  Coefs<- data.frame(Slope = c(coef(Mods[[x]])[2], #create a data frame of the model coefficients for slopes
                       coef(Mods[[x]])[(GFLevels+2):(GFLevels+GFLevels)] +
                         coef(Mods[[x]])[2]),
             rbind(confint(Mods[[x]])[2,], #with  95% confidence intervals
                   confint(Mods[[x]])[(GFLevels+2):(GFLevels+GFLevels),] +
                     confint(Mods[[x]])[2,]))
  Ints <- data.frame(Intercept = c(coef(Mods[[x]])[1],  #create a data frame of the model coefficients for intercepts
                                   coef(Mods[[x]])[3:(GFLevels+1)] +
                                     coef(Mods[[x]])[1]),
                        rbind(confint(Mods[[x]])[1,],  #with  95% confidence intervals
                              confint(Mods[[x]])[3:(GFLevels+1),] +
                                confint(Mods[[x]])[1,]))
  Coefs$SE <- c(summary(Mods[[x]])$coefficients[2,2],summary(Mods[[x]])$coefficients[(GFLevels+2):(GFLevels+GFLevels),2]) #Standard errrors
  Coefs$T <- Coefs$Slope/Coefs$SE #T scores
  Coefs$N <- aggregate(Sphericity ~ GrowthForm, PCADat,length)[,2] #Growth form counts
  Coefs$P <- 2*pt(-abs(Coefs$T), df = 124) #P values
  Coefs$UPR <- Coefs$Slope + (1.96*Coefs$SE) #Upper slope confidence interval
  Coefs$LWR <- Coefs$Slope - (1.96*Coefs$SE) #Lower slope confidence interval
  Coefs$DirectionConfint <- factor(ifelse(Coefs$LWR>0, #If the lower bound is greater than zero
                                          "positive", #mark as positive
                                          ifelse(Coefs$UPR< 0, #else, if the upper bound is less than zero
                                                            "negative","zero")), #masr as negative, otherwise mark as zero
                                   levels = c("positive","negative","zero")) #set levels for factor
  Coefs$DirectionP <- factor(ifelse(Coefs$P >0.05, #If P > 0.05
                                    "P >0.05", #mark as p > 0.05
                                    ifelse(Coefs$Slope > 0, #else, if the sloep is greater than zero
                                                     "positive", #mark as positive
                                           "negative")), #else, mark as negative.
                             levels = c("P >0.05","negative","positive")) #set levels for factor
  Coefs$GrowthForm <- levels(PCADat$GrowthForm) #set levels for growth form.
  return(Coefs) #Return coefficients data frame
})

PredsList <- lapply(1:length(Vars), function(x) { #for each variable
  Preds <- predict(Mods[[x]], newdata = NewData, interval = "confidence") #predict values using model and new data
  Preds <- data.frame(NewData,Preds) #create combined dataframe
  Preds <- merge.data.frame(Preds,CoefsList[[x]], by = "GrowthForm") #merge predicted values and coefficients data via growth form
})

#..........plots ====

YLabs <- c(expression("logit S"),
           expression(log^"10"*" C"),
           expression(log^"10"*" P"),
           "FD",
           expression(atop(log^"10"*" V"["VOL"]," mm"^"4")),
           expression(atop(log^"10"*" V"["AREA"]," mm"^"3")))

VarVolPlots <- lapply(1:length(Vars), function(x) { #for each variable
  ggplot(PCADat, aes(y = get(Vars[x]), x = Volume)) + #ggplot set up
    PlotTheme + #set plot theme
    scale_fill_manual(values = c("blue","red","black"), drop = F) + #set colours
    scale_colour_manual(values = c("black","red","blue"), drop = F) + #set colours DEPRECATED
    geom_point(size = 4, colour = "black", fill = "grey", alpha =0.8,shape = 21) + #data points
    geom_ribbon(data = PredsList[[x]], #confidence bands
                aes(x = Volume,
                    y = fit,
                    ymin = lwr,
                    ymax = upr,
                    fill = DirectionConfint),
                alpha = 0.5) + #, colour = SpherCoefs$Direction)) +
    geom_line(data = PredsList[[x]], aes(x = Volume, y = fit), size = 2) + #regression line
    guides(colour = F, fill = F, size = F) + #remove legends
    theme(legend.title = element_blank(),
          plot.margin = unit(c(0,2,0,2), "cm"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          #axis.ticks.length = unit(-0.25,"cm"),
          #panel.grid.major.x = element_line( size=.1, color="black" ),
          strip.text = element_text(face = "plain"),
          strip.background = element_rect(colour = "black"),
          panel.grid = element_blank(),
          panel.spacing = unit(0,"cm"),
          panel.border = element_rect(colour = "black")) +
    scale_y_continuous(expand = c(0.125,0))+
    labs(y = YLabs[x], x = NULL) + #set y axis labels
    coord_cartesian(xlim = c(3.5,7.5)) +
    facet_wrap(~GrowthForm, ncol=length(GFs)) #facet by growth form
})

VarVolPlots[[length(Vars)]] <- VarVolPlots[[length(Vars)]] +
  guides(fill = "legend") #add legend to bottom panel

VarVolPlots[[length(Vars)]] <- VarVolPlots[[length(Vars)]] +
  xlab(expression(log^"10"*" volume "*"mm"^"3")) +
  theme(axis.text.x = element_text(),
        axis.ticks.x = element_line())  #add x axis labels to bottom panel

for(x in 2:length(Vars)){
  VarVolPlots[[x]] <- VarVolPlots[[x]] +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())
}

VarVolAll <- grid.arrange(rbind(ggplotGrob(VarVolPlots[[1]]),
                             ggplotGrob(VarVolPlots[[2]]),
                             ggplotGrob(VarVolPlots[[3]]),
                             ggplotGrob(VarVolPlots[[4]]),
                             ggplotGrob(VarVolPlots[[5]]),
                             ggplotGrob(VarVolPlots[[6]]),
                             size = "last"))

ggsave(file=paste0(outputdir,"/VarVolAll.png"), plot=VarVolAll, width=3.45, height=3.45,scale = 6)
ggsave(file=paste0(outputdir,"/VarVolAll.pdf"), plot=VarVolAll, width=3.45, height=3.45,scale = 6)


#multinomial models ----
#......analysis ====

#............Models ====

PredDat <- PCADat #dataset for predicting growth form from shape data
PredDat$GrowthForm <- droplevels(PredDat$GrowthForm)
ModNull <-  multinom(formula = GrowthForm ~ 1, data = PredDat, maxit = 10000) #null model

#Full model
ModSelVars <- multinom(formula = GrowthForm ~ #model that balances accuracy with variable reduction
                               #Volume : Sphericity +
                               #Sphericity +
                               #Volume : Convexity +
                               Convexity +
                               Volume : Packing +
                               Packing +
                               Volume : FractalDimension +
                               FractalDimension +
                               #Volume : Vertical2ndMomentVolumeScaled +
                               Vertical2ndMomentVolumeScaled,
                               #Volume : Vertical2ndMomentAreaScaled +
                               #Vertical2ndMomentAreaScaled,# +
                             data = PredDat,
                             maxit = 10000)
AIC(ModAllBestAIC,ModSelVars)

#Leave one observation out models
LeaveOneOutModels <- lapply(1:nrow(PredDat), function(x) {
  multinom(formula = GrowthForm ~ #model with best AIC using all six available variables
            #Volume : Sphericity +
            #Sphericity +
            #Volume : Convexity +
            Convexity +
            Volume : Packing +
            Packing +
            Volume : FractalDimension +
            FractalDimension +
            #Volume : Vertical2ndMomentVolumeScaled +
            Vertical2ndMomentVolumeScaled
            #Volume : Vertical2ndMomentAreaScaled +
            #Vertical2ndMomentAreaScaled,# +
            ,
           data = PredDat[-x,],
           maxit = 10000)
})

# Predicted probabilites
PredVals <- lapply(1:(length(LeaveOneOutModels)), function(x)
  predict(LeaveOneOutModels[[x]],
          newdata = PredDat[x,],
          type = "prob",
          se.fit = T,
          confint = "confidence"))
PredVals <- unlist(PredVals)

PredValsDF <- data.frame(ActualGF = PredDat$GrowthForm, matrix(round(PredVals,3), nrow=nrow(PredDat), byrow=T))

GFLevels <- length(levels(LZDat$GrowthForm)) #length of the number of levels in the growth form data

colnames(PredValsDF)[2:(1+GFLevels)] <- levels(droplevels(PredDat$GrowthForm))

PredValsLong <- melt(PredValsDF, id.vars= "ActualGF",
                     measure.vars = levels(PredDat$GrowthForm),
                     variable.name = "PredictedGF")

PredValsLongMeans <- aggregate(value ~ ActualGF + PredictedGF, PredValsLong, mean)
PredValsLongMeans$n <- aggregate(value ~ ActualGF + PredictedGF, PredValsLong, length)[,3]
PredValsLongMeans$sd <- aggregate(value ~ ActualGF + PredictedGF, PredValsLong, sd)[,3]
PredValsLongMeans$se <- PredValsLongMeans$sd/sqrt(PredValsLongMeans$n)
PredValsLongMeans$confint <- 1.96*PredValsLongMeans$se


# predicted classes
PredClass <- lapply(1:(length(LeaveOneOutModels)), function(x)
  predict(LeaveOneOutModels[[x]],
          newdata = PredDat[x,]))
PredClass <- unlist(PredClass)

# Confusion matrix
confusionMatrix(PredClass, PredDat$GrowthForm)

ProbPlot <- ggplot(data = PredValsLongMeans, aes(x = ActualGF, y = value, fill = PredictedGF, ymin = value-se,ymax = value+se)) +
  geom_bar(position = position_dodge(.9),colour = "black", stat = "identity", alpha = 0.9) +
  geom_errorbar(position = position_dodge(.9), width = .5) +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Observed growth form",
       y = "Probability",
       fill = element_blank()) +
  geom_abline(slope = 0, intercept = 0.14, col = "black", linetype = "dashed") +
  scale_color_manual(values = GFNinePallete) +
  scale_fill_manual(values = GFNinePallete) +
  scale_y_continuous(expand = c(0,0))+
  PlotTheme +
  theme(axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(title.position="left",
                             title.hjust = 0.5,
                             nrow = 1,
                             title.theme = element_text(size = 25, angle = 0),
                             override.aes = list(colour = NA))
         ) +
  labs(fill = "Predicted growth form:")

ggsave(file=paste0(outputdir,"/ProbPlot.png"), plot=ProbPlot, width=3.45, height=2,scale = 6)
ggsave(file=paste0(outputdir,"/ProbPlot.pdf"), plot=ProbPlot, width=3.45, height=2,scale = 6)




