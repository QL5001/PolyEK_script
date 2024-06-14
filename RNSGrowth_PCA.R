#Principal component analysis (PCA) on the RNASeq data of Enterococcus faecalis strain for determining its growth-associated transcriptomic changes tirggered by Klebsiella pneumoniae strain in CAUTI polymicrobial community.
#Beign#

library(knitr)
library(mixOmics)
library(dplyr)
library(stats)

RNS_Growth = read.csv(file = file.choose(), sep = ",", header = TRUE, stringsAsFactors = FALSE)    #open csv file "RNSGrowth_PCA" containing gene differential expression data of E. faecalis strain

X = RNS_Growth[, 1:2947]  
X_mean = apply(X, 2, mean)
X_center = scale(X, center = X_mean, scale = FALSE)
Y = RNS_Growth$Sample

pca.RNS = pca(X_center, ncomp = 6, center = FALSE, scale = FALSE)
plot(pca.RNS, ylim = c(0, 0.8), main = "PCA on WGS")

ev = pca.RNS$prop_expl_var   #explained_variance
print(ev)
PCA_ev = data.frame(ev)
write.csv(PCA_ev, "RNSGrowth_PCA_ev.csv")  #export explained variance for each principal component into a csv file

sco = pca.RNS$variates   #prinicipal component values
print(sco)
sco_PC = data.frame(SampleNo = RNS_Growth$Sample, sco)
write.csv(sco_PC, "RNS_Growth_PCA_PC.csv")  #export the principal component values into a csv file

#Exported csv files will be used for further data analyses and drawing figures in Excel, Prism, and R programming
#End#
