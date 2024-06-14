#Principal component analysis (PCA) on the RNASeq data of Enterococcus faecalis OG1RF strain for determining its bacterial growth-associated transcriptomic changes under the treatment of antimicrobial compounds.
#Beign#
library(knitr)
library(mixOmics)
library(dplyr)
library(stats)


WGS_biofilm_0410 = read.csv(file = file.choose(), sep = ",", header = TRUE, stringsAsFactors = FALSE)    #open file "group_top52_9h_Differential_Expression-PCA_240122"

X = WGS_biofilm_0410[, 1:2947]  
X_mean = apply(X, 2, mean)
X_center = scale(X, center = X_mean, scale = FALSE)
Y = WGS_biofilm_0410$Sample

pca.WGS = pca(X_center, ncomp = 6, center = FALSE, scale = FALSE)

plot(pca.WGS, ylim = c(0, 0.8), main = "PCA on WGS")


ev = pca.WGS$prop_expl_var   #explained_variance
print(ev)
PCA_ev = data.frame(ev)
write.csv(PCA_ev, "PCA_ev_240122.csv")


sco = pca.WGS$variates   #prinicipal component values
print(sco)
sco_PC = data.frame(SampleNo = WGS_biofilm_0410$Sample, sco)
write.csv(sco_PC, "PCA_PC_240122.csv")
