###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: PCA to explore global methylation & determine associations with covariates

# input: swan normalised beta values (output from 01preprocessingSWAN.R)
#        conducted on each cohort separately to identify batch effects 
#        conducted on the full dataset (IVF & naturally conceived data)

# output: PCA plot
#         statistics for associations between sample characteristics & PCs
#         heatmap showing associations between sample characteristics & PCs

###########################################################################################################################

#load packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(jmuOutlier))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(tidyr))

# set data directories
base.dir <- "output"

annotation.file <- file.path(base.dir, "annotation.csv")

betas.file <- file.path(base.dir, "betas.csv")

# load the data files

annotation <- fread(annotation.file) %>% mutate(Sample_ID = as.character(Sample_ID))

betas <- fread(betas.file) %>% 
  select(ID, all_of(annotation$Sample_ID)) %>% column_to_rownames("ID")

betas.pca <- betas %>% na.omit() %>%
  as.matrix() %>%
  t()

# conduct a principal component analysis
pca <- prcomp(betas.pca, center = T, scale. = F)

summary(pca)$importance[, 1:8]

coords <- as.data.frame(pca$x) %>% rownames_to_column("Sample_ID") %>%
  full_join(annotation, .[ , 1:15], by = "Sample_ID")

coords$`Sample Plate` <- if_else(coords$Sample_Plate == "WG5839007-BCD", "1", "2")

coords$`Culture Medium` <- if_else(coords$MEDIUM == "htf", "HTF", "G5")

pca.plt <- coords %>% ggplot(aes(x = PC1, y = PC2, colour = `Culture Medium`)) + 
  geom_point(size = 0.6) +
  scale_colour_manual(values = c("darkgoldenrod3", "deepskyblue2")) +
  theme_classic() +
  xlab("PC1 (24%)") +
  ylab("PC2 (9%)") 

######################### look for associations between the PCs (1 to 8) and technical / demographic features #################

# sample plate (binary categorical = Wilcoxon rank)

plate = list()

for(n in 1:8){
  
  pc = paste0("PC", n)
  
  x = coords %>% filter(Sample_Plate == "WG5839007-BCD") %>% select(all_of(pc))
  
  x = x[,1]
  
  y = coords %>% filter(Sample_Plate != "WG5839007-BCD") %>% select(all_of(pc)) %>% as.vector()
  
  y = y[,1]
  
  test = wilcox.test(x = x, y = y, alternative = "two.sided")
  
  plate[[pc]] = test$p.value
  
}

# gender (binary categorical = Wilcoxon rank)

gender = list()

for(n in 1:8){
  
  pc = paste0("PC", n)
  
  x = coords %>% filter(gender == "M") %>% select(all_of(pc))
  
  x = x[,1]
  
  y = coords %>% filter(gender == "F") %>% select(all_of(pc)) %>% as.vector()
  
  y = y[,1]
  
  test = wilcox.test(x = x, y = y, alternative = "two.sided")
  
  gender[[pc]] = test$p.value
  
}

# culture medium (binary categorical = Wilcoxon rank)

medium = list()

for(n in 1:8){
  
  pc = paste0("PC", n)
  
  x = coords %>% filter(MEDIUM == "vg5") %>% select(all_of(pc))
  
  x = x[,1]
  
  y = coords %>% filter(MEDIUM == "htf") %>% select(all_of(pc)) %>% as.vector()
  
  y = y[,1]
  
  test = wilcox.test(x = x, y = y, alternative = "two.sided")
  
  medium[[pc]] = test$p.value
  
}

# gestational age (continuous - correlation permutation test)

age = list()
cor.age = list()

for(n in 1:8){
  
  pc = paste0("PC", n)
  
  data = coords[ , c("gestational_age", pc)]
  
  test = perm.cor.test(x = data[,1], y = data[,2], alternative = "two.sided", method = "pearson", num.sim = 10000)
  
  age[[pc]] = test$p.value
  
  cor.age[[pc]] = cor(data[,1], data[,2], method = "pearson")
  
}

# maternal age (continuous - correlation permutation test)

mat_age = list()
mat_cor.age = list()

for(n in 1:8){
  
  pc = paste0("PC", n)
  
  data = coords[ , c("maternal_age", pc)]
  
  test = perm.cor.test(x = data[,1], y = data[,2], alternative = "two.sided", method = "pearson", num.sim = 10000)
  
  mat_age[[pc]] = test$p.value
  
  mat_cor.age[[pc]] = cor(data[,1], data[,2], method = "pearson")
  
}

# sample cell composition (continuous - correlation permutation test)

cd8t = list()

cor.cd8t = list()

for(n in 1:8){
  
  pc = paste0("PC", n)
  
  data = coords[ , c("counts.CD8T", pc)]
  
  test = perm.cor.test(x = data[,1], y = data[,2], alternative = "two.sided", method = "pearson", num.sim = 10000)
  
  cd8t[[pc]] = test$p.value
  
  cor.cd8t[[pc]] = cor(data[,1], data[,2], method = "pearson")
  
}

cd4t = list()

cor.cd4t = list()

for(n in 1:8){
  
  pc = paste0("PC", n)
  
  data = coords[ , c("counts.CD4T", pc)]
  
  test = perm.cor.test(x = data[,1], y = data[,2], alternative = "two.sided", method = "pearson", num.sim = 10000)
  
  cd4t[[pc]] = test$p.value
  
  cor.cd4t[[pc]] = cor(data[,1], data[,2], method = "pearson")
  
}

nk = list()

cor.nk = list()

for(n in 1:8){
  
  pc = paste0("PC", n)
  
  data = coords[ , c("counts.NK", pc)]
  
  test = perm.cor.test(x = data[,1], y = data[,2], alternative = "two.sided", method = "pearson", num.sim = 10000)
  
  nk[[pc]] = test$p.value
  
  cor.nk[[pc]] = cor(data[,1], data[,2], method = "pearson")
  
}

bcell = list()

cor.bcell = list()

for(n in 1:8){
  
  pc = paste0("PC", n)
  
  data = coords[ , c("counts.Bcell", pc)]
  
  test = perm.cor.test(x = data[,1], y = data[,2], alternative = "two.sided", method = "pearson", num.sim = 10000)
  
  bcell[[pc]] = test$p.value
  
  cor.bcell[[pc]] = cor(data[,1], data[,2], method = "pearson")
  
}

mono = list()

cor.mono = list()

for(n in 1:8){
  
  pc = paste0("PC", n)
  
  data = coords[ , c("counts.Mono", pc)]
  
  test = perm.cor.test(x = data[,1], y = data[,2], alternative = "two.sided", method = "pearson", num.sim = 10000)
  
  mono[[pc]] = test$p.value
  
  cor.mono[[pc]] = cor(data[,1], data[,2], method = "pearson")
  
}

gran = list()

cor.gran = list()

for(n in 1:8){
  
  pc = paste0("PC", n)
  
  data = coords[ , c("counts.Gran", pc)]
  
  test = perm.cor.test(x = data[,1], y = data[,2], alternative = "two.sided", method = "pearson", num.sim = 10000)
  
  gran[[pc]] = test$p.value
  
  cor.gran[[pc]] = cor(data[,1], data[,2], method = "pearson")
  
}

nrbc = list()

cor.nrbc = list()

for(n in 1:8){
  
  pc = paste0("PC", n)
  
  data = coords[ , c("counts.nRBC", pc)]
  
  test = perm.cor.test(x = data[,1], y = data[,2], alternative = "two.sided", method = "pearson", num.sim = 10000)
  
  nrbc[[pc]] = test$p.value
  
  cor.nrbc[[pc]] = cor(data[,1], data[,2], method = "pearson")
  
}

# sentrix ID (categorical with >2 groups - Kruskall-Wallis)

sentID = list()

for(n in 1:8){
  
  pc = paste0("PC", n)
  
  data = coords[ , c("Sentrix_ID", pc)]
  
  test = kruskal.test(data[,1]~data[,2])
  
  sentID[[pc]] = test$p.value
  
}

# sentrix position (categorical with >2 groups - Kruskall-Wallis)

sentPos = list()

for(n in 1:8){
  
  pc = paste0("PC", n)
  
  data = coords[ , c("Sentrix_Position", pc)]
  
  test = kruskal.test(data[,1]~data[,2])
  
  sentPos[[pc]] = test$p.value
  
}

# combine the statistical tests for associations to allow plotting

associations = data.frame(row.names = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")) %>%
  rownames_to_column("Principal_component") %>%
  mutate(Sample_plate = unlist(plate), 
         Sentrix_ID = unlist(sentID),
         Sentrix_Position = unlist(sentPos),
         Sex = unlist(gender),
         Gestational_age = unlist(age),
         Maternal_age = unlist(mat_age),
         Culture_Medium = unlist(medium),
         CD8T_cells = unlist(cd8t),
         CD4T_cells = unlist(cd4t),
         NK_cells = unlist(nk),
         B_cells = unlist(bcell),
         Monocytes = unlist(mono),
         Granulocytes = unlist(gran),
         nRBCs = unlist(nrbc))

# adjust the data for plotting

to.plot = associations %>% column_to_rownames("Principal_component") %>%
  t() %>% -log10(. + 10^-8)  

t = to.plot %>% as.data.frame() %>% rownames_to_column("V1")

t2 = gather(t, key = "PC", value = "log.pvals", PC1:PC8)

t2$V1 = factor(t2$V1, levels = rev(rownames(to.plot)), ordered = T)

t2 = t2 %>% mutate(p.vals = 10^-log.pvals, text = p.vals <= 0.05, pval.text = formatC(p.vals, format = "e", digits = 1))

# a very small value is added to avoid getting Inf values when calculating the -log10 of p-values that are exactly 0

# plot the results

plot2 = t2 %>% ggplot(aes(x = PC, y = V1, fill = log.pvals)) + 
  geom_tile(colour = "white") + 
  scale_fill_gradientn(colours = c("white", "lightblue", "dodgerblue2"), 
                       values = scales::rescale(c(0, -log10(0.05), 8)), limits = c(0,5), oob = scales::squish) +
  geom_text(data = function(x){filter(x, text)}, aes(label = pval.text), size = 2, angle = 20) +
  theme_classic() +
  labs(fill = "-log10 \n p-value") +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(), 
        text = element_text(family = "Helvetica", size = 7),
        plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 7),
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 7),
        axis.text.y = element_text(size = 7))
