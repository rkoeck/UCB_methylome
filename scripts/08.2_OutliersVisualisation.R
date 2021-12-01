###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: visualisation of the number of outliers per sample

# input: output from 08Outliers.R (summarised total number of outliers per sample)
#        conducted to compare G5 and HTF within the IVF cohort (including and excluding pregnancy complications)

# output: multi-faceted plot containing scatter plot of hypo to hyper methylation outliers in the centre with adjacent distribution summaries

# source code origin: https://rdrr.io/github/jorvlan/openvis/src/R/R_rainclouds.R

###########################################################################################################################

# load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))


source("rainCloudPlot.R")

# plot formatting

plot.format = theme(text = element_text(family = "Helvetica", size = 7),
                    plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
                    axis.text.x = element_text(size = 7),
                    legend.title = element_text(size = 7, face = "bold"),
                    legend.text = element_text(size = 7),
                    axis.text.y = element_text(size = 7),
                    axis.title.x = element_text(size = 7),
                    axis.title.y = element_text(size = 7),
                    plot.margin = margin())

# load the data

annotation = fread("annotation.csv") %>%
  mutate(Sample_ID = as.character(Sample_ID))

outliers = fread("sampleThresholdOutliers.csv",
                 header = TRUE) %>% 
  column_to_rownames("V1") %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("Sample_ID") %>%
  mutate(Sample_ID = as.character(Sample_ID))

# merge the outlier values with the annotation data

outliers = full_join(annotation, outliers, by = "Sample_ID")

outliers$MEDIUM = if_else(outliers$MEDIUM == "vg5", "G5", "HTF")

outliers$MEDIUM = factor(outliers$MEDIUM, levels = c("HTF", "G5"))

# visualise the results as a scatter plot with density information summarised at the side 

scatter = outliers %>% ggplot(aes(x = hypo, y = hyper, colour= MEDIUM)) +
  geom_point(size = 0.3) +
  scale_x_log10(breaks = c(100, 1000, 10000), labels = c( "100", "1000", "10000"), limits = c(49, 40000), expand = expansion(add = c(0.1, 0.3))) +
  scale_y_log10(breaks = c(100, 1000, 10000), labels = c( "100", "1000", "10000"), limits = c(49, 40000), expand = expansion(add = c(0.1, 0.3))) +
  theme_classic() +
  scale_fill_manual(values = c(HTF = "deepskyblue2", G5 = "darkgoldenrod3")) +
  scale_colour_manual(values = c(HTF = "deepskyblue2", G5 = "darkgoldenrod3" )) +
  guides(colour = guide_legend(override.aes = list(size = 1))) +
  plot.format +
  theme(legend.title.align = 0.5,
        legend.background = element_rect(color = "black")) +
  ylab("Hypermethylation outliers per sample") +
  xlab("Hypomethylation outliers per sample")

x = outliers %>% ggplot(aes(x = MEDIUM, y = hypo, fill = MEDIUM, colour = MEDIUM)) +
  geom_flat_violin(position = position_nudge(x = -0.2, y = 0), adjust = 0.75, trim = TRUE) +
  geom_boxplot(aes(x = as.numeric(MEDIUM) -0.2, y = hypo), outlier.shape = NA, alpha = 0.5, width = 0.1, colour = "black", lwd = 0.3) +
  theme_classic() + 
  scale_y_log10(breaks = c(100, 1000, 10000), labels = c( "100", "1000", "10000"), limits = c(49, 40000), expand = expansion(add = c(0.1, 0.3))) +
  coord_flip() +
  scale_fill_manual(values = c("deepskyblue2", "darkgoldenrod3")) +
  scale_colour_manual(values = c("deepskyblue2", "darkgoldenrod3" )) +
  scale_shape_manual(values = c(19, 15, 1, 17)) +
  xlab(NULL) +
  ylab("log10(number of outliers per sample)") +
  plot.format +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank())

y = outliers %>% ggplot(aes(x = MEDIUM, y = hyper, fill = MEDIUM, colour = MEDIUM)) +
  geom_flat_violin(position = position_nudge(x = -0.2, y = 0), adjust = 0.75, trim = TRUE) +
  geom_boxplot(aes(x = as.numeric(MEDIUM) -0.2), outlier.shape = NA, alpha = 0.5, width = 0.1, colour = "black", lwd = 0.3) +
  theme_classic() +
  scale_y_log10(breaks = c(100, 1000, 10000), labels = c( "100", "1000", "10000"), limits = c(49, 40000), expand = expansion(add = c(0.1, 0.3))) +
  scale_fill_manual(values = c("deepskyblue2", "darkgoldenrod3")) +
  scale_colour_manual(values = c("deepskyblue2", "darkgoldenrod3" )) +
  scale_shape_manual(values = c(19, 15, 1, 17)) +
  xlab(NULL) +
  ylab("log10(number of outliers per sample)") +
  #theme(legend.position = "none") +
  plot.format +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) 

pdf(file = file.path(fig.dir, "altOutliers.pdf"),
    width = 3.34, height = 3.34)
x + guide_area() + scatter + y + 
  plot_layout(ncol = 2,
              heights = c(1,3),
              widths = c(3, 1),
              guides = "collect")
dev.off()





