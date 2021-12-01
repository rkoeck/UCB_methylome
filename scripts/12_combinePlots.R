###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: script to compose (multi-panel) figures, make them uniform and save them in an appropriate format for submission

# source code origin: https://rdrr.io/github/jorvlan/openvis/src/R/R_rainclouds.R

###########################################################################################################################

# load packages

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(cowplot))

source("rainCloudPlot.R")

# set directory

fig.dir = "figures/"

# function to generalise formatting (fonts etc)

plot.format = theme(text = element_text(family = "Helvetica", size = 7),
                    plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
                    axis.text.x = element_text(size = 7),
                    legend.title = element_text(size = 7, face = "bold"),
                    legend.text = element_text(size = 7),
                    axis.text.y = element_text(size = 7),
                    axis.title.x = element_text(size = 7),
                    axis.title.y = element_text(size = 7))
  
####################### global methylation #############################

# input: PCA plot from 04PCACovariates.R

p1 = readRDS(file.path(fig.dir, "globalPCA.rds")) + plot.format +
  theme(legend.position = "bottom",
        legend.margin = margin(t = -10))

pdf(width = 3.34, height = 3, file = file.path(fig.dir, "Figure2.pdf"))
p1
dev.off()

jpeg(file.path(fig.dir, "Figure1.jpg"),
     width = 3.34, height = 3, units = "in",
     res = 900)
p1
dev.off()

##################### DNA methylation at individual CpGs #############

# input: plots generated with 07VolcanoPlot.R

p2 = readRDS(file.path(fig.dir, "dmpComplexImprinted.rds")) + plot.format

p3 = readRDS(file.path(fig.dir, "dmpComplexBwt.rds")) + plot.format

jpeg(file.path(fig.dir, "Figure3.jpg"),
     width = 4.5, height = 3, units = "in",
     res = 900)
p2 + p3 + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 8, face = "bold"))
dev.off()


pdf(width = 4.5, height = 3, file = file.path(fig.dir, "Figure3.pdf"))
p2 + p3 + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 8, face = "bold"))
dev.off()

##################### regional DNA methylation ######################

# input: plots generated with 07VolcanoPlot.R

p4 = readRDS(file.path(fig.dir, "genesComplexImprinted.rds")) + ggtitle("Genes") + plot.format + theme(legend.position = "none")

p5 = readRDS(file.path(fig.dir, "promoterComplex.rds")) + plot.format + theme(legend.position = "none")

p6 = readRDS(file.path(fig.dir, "islandComplicated.rds")) + plot.format + theme(legend.position = "none")

jpeg(file.path(fig.dir, "Figure4.jpg"),
     width = 6, height = 3, units = "in",
     res = 900)
p4 + p5 + p6 + plot_annotation(tag_levels = list(c("C", "D", "E"))) & theme(plot.tag = element_text(size = 8, face = "bold"))
dev.off()

pdf(width = 10, height = 5, file = file.path(fig.dir, "Figure4.pdf"))
p4 + p5 + p6 + plot_annotation(tag_levels = list(c("C", "D", "E"))) & theme(plot.tag = element_text(size = 8, face = "bold"))
dev.off()

################### outliers per sample #######################

# input: plot generated with 08OutliersVisualisation.R

p7 = readRDS(file.path(fig.dir, "outliersFull.rds")) + plot.format

jpeg(file.path(fig.dir, "Figure5.jpg"),
     width = 3.34, height = 2.6, units = "in",
     res = 300)
p7
dev.off()


pdf(width = 3.34, height = 2.6, file = file.path(fig.dir, "Figure5.pdf"))
p7
dev.off()

################### GAA ############################################

p8 = readRDS(file.path(fig.dir, "GAAAllIVF.rds")) + plot.format

pdf(width = 3.34, height = 2.6, file = file.path(fig.dir, "Figure6.pdf"))
p8
dev.off()

################# preprocessing ###############################

# input: plots from 02SexPredictionSest.R, 03CellCompositionVisualisation.R, 04PCACovariates.R

s1 = readRDS(file.path(fig.dir, "sEST.rds")) + plot.format

s2 = readRDS(file.path(fig.dir, "idol.rds")) + plot.format +
  theme(legend.position = "bottom")
  
s3 = readRDS(file.path(fig.dir, "associationsHeatmap.rds")) 

left = plot_grid(s1, s3, ncol = 1, rel_heights = c(1.5,2), 
          labels = c("A", "C"),
          label_size = 8,
          label_fontface = "bold")
  
pdf(width = 6.69, height =6, file = file.path(fig.dir, "SuppFigure1.pdf"))
plot_grid(left, s2, ncol = 2,
          rel_widths = c(2, 1),
          labels = c("", "B"),
          label_size = 8,
          scale = 0.95)
dev.off()
  
  
############# uncomplicated DMPs ###########################

# input: plots from 07VolcanoPlot.R

s4 = readRDS(file.path(fig.dir, "dmpUncomplicatedImprinted.rds")) + plot.format

s5 = readRDS(file.path(fig.dir, "dmpUncomplicatedBwt.rds")) + plot.format

jpeg(file.path(fig.dir, "SuppFigure2.jpg"),
     width = 4.5, height = 3, units = "in",
     res = 300)
s4 + s5 + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 8, face = "bold"))
dev.off()

pdf(width = 4.5, height = 3, file = file.path(fig.dir, "SuppFigure2.pdf"))
s4 + s5 + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 8, face = "bold"))
dev.off()

########## uncomplicated DMRs #############################

# input: plots from 07VolcanoPlot.R

s6 = readRDS(file.path(fig.dir, "genesUncomplicatedImprinted.rds")) + ggtitle("Genes") + plot.format + theme(legend.position = "none")

s7 = readRDS(file.path(fig.dir, "promoterUncomplicated.rds")) + plot.format + theme(legend.position = "none")

s8 = readRDS(file.path(fig.dir, "islandUncomplicated.rds")) + plot.format + theme(legend.position = "none")

jpeg(file.path(fig.dir, "SuppFigure3.jpg"),
     width = 6, height = 3, units = "in",
     res = 300)
s6 + s7 + s8 + plot_annotation(tag_levels = list(c("C", "D", "E"))) & theme(plot.tag = element_text(size = 8, face = "bold"))
dev.off()

pdf(width = 6, height = 3, file = file.path(fig.dir, "SuppFigure3.pdf"))
s6 + s7 + s8 + plot_annotation(tag_levels = list(c("C", "D", "E"))) & theme(plot.tag = element_text(size = 8, face = "bold"))
dev.off()

######### DNA methylation variance ########################

# input: plots from 08OutliersVisualisation.R

s9 = readRDS(file.path(fig.dir, "outliersFullHypoHyper.rds")) + plot.format

s10 = readRDS(file.path(fig.dir, "outliersUncomplicated.rds")) 

s11 = readRDS(file.path(fig.dir, "outliersUncomplicatedHypoHyper.rds")) + plot.format + theme()

u = plot_grid(s10 + theme(axis.text.x = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  axis.title.x = element_blank()) +
            ggtitle("Total outliers"),
          s11, ncol = 1, rel_heights = c(1,2))

pdf(width = 6.69, height = 5.5, file = file.path(fig.dir, "SuppFigure4.pdf"))
s9 + u + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 8, face = "bold"))
dev.off()

###### GAA uncomplicated ##################################

s12 = readRDS(file.path(fig.dir, "GAAUncomplicated.rds")) + plot.format

pdf(width = 3.34, height = 2.6, file = file.path(fig.dir, "SuppFigure5.pdf"))
s12
dev.off()
