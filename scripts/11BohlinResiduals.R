# Calculation of the residuals after Bohlin eGA estimation & exploration of the results

# input: annotated annotation file generated with 11BohlinEstimateGA.R (must also include cell composition and covariates)
#         epigenetic gestational age was only calculated for samples without pregnancy complications


# output: statistics relating to the residuals, residuals plotted as a raincloud plot

# method replicated from: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1063-4
# source code origin: https://rdrr.io/github/jorvlan/openvis/src/R/R_rainclouds.R

# load packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(lmerTest))
suppressPackageStartupMessages(library(Metrics))
suppressPackageStartupMessages(library(jmuOutlier))

source("rainCloudPlot.R")

# set the directories

base.dir <- "output"

bohlin.file <- file.path(base.dir, "estimatedGABohlin.csv")

# load the gestational age estimates and bind these to the annotation file

estimateGA <- fread(bohlin.file)

# determine the root mean squared error of the model

rmse(estimateGA$gestational_age, estimateGA$predictedGA)

# calculate the residuals correcting for GA and also cell composition

acc.model <- lm(predictedGA ~ gestational_age + counts.CD8T + counts.CD4T + counts.NK + counts.Mono + counts.Gran + counts.nRBC + counts.Bcell, data = estimateGA)

residuals <- acc.model$residuals

estimateGA$residuals <- residuals

# see if there are any group differences

stats.test = lmer(formula = residuals ~ MEDIUM + gender + maternal_age + (1|center), data = estimateGA)

summary(stats.test)

# calculate average GAA per group
vg5 = estimateGA %>% filter(MEDIUM == "vg5") %>% .$residuals
htf = estimateGA %>% filter(MEDIUM == "htf") %>% .$residuals


mean(htf)
sd(htf)

mean(vg5)
sd(vg5)

# visualise the residuals by medium 

estimateGA$MEDIUM = if_else(estimateGA$MEDIUM == "vg5", "G5", "HTF")

estimateGA$MEDIUM = factor(estimateGA$MEDIUM, levels = c("HTF", "G5"))


plot.gaa = estimateGA %>% ggplot(aes(x = MEDIUM, y = residuals, fill = MEDIUM, colour = MEDIUM)) +
  geom_flat_violin(position = position_nudge(x = 0, y = 0), adjust = 0.75, trim = FALSE) +
  geom_point( aes(x = as.numeric(MEDIUM) - 0.2, y = residuals), position = position_jitter(width = 0.1), size = 0.3) +
  geom_boxplot(aes(x = as.numeric(MEDIUM) + 0.0, y = residuals), outlier.shape = NA, alpha = 0.3, width = 0.1, colour = "black", lwd = 0.3) +
  theme_classic() + 
  coord_flip() +
  scale_fill_manual(values = c("deepskyblue2", "darkgoldenrod3")) +
  scale_colour_manual(values = c("deepskyblue2", "darkgoldenrod3" )) +
  xlab(NULL) +
  ylab("Gestational age acceleration (weeks)") +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16))


######################## look for an association with birth weight ####################

cor(estimateGA$residuals, estimateGA$Birthweight, method = "pearson")

estimateGA %>% ggplot(aes(x = Birthweight, y = residuals)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "lm", se = FALSE) +
  ylab("Gestational age acceleration (weeks)") +
  xlab("Birth weight (g)")

ggsave(filename = file.path(base.dir, "/estimateGA/bohlinWeightResiduals.png"), plot = wt.residuals,
       height = 12, width = 15, unit = "cm")


x = estimateGA$Birthweight

y = estimateGA$residuals

test = perm.cor.test(x = x, y = y, alternative = "two.sided", method = "pearson", num.sim = 10000)

test$p.value

