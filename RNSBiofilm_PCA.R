#Principal component analysis (PCA) on the RNASeq data of Enterococcus faecalis strain for determining its biofilm-associated transcriptomic changes tirggered by Klebsiella pneumoniae strain in CAUTI polymicrobial community.
#Beign#

library(knitr)
library(mixOmics)
library(dplyr)
library(stats)

RNS_Bioilm = read.csv(file = file.choose(), sep = ",", header = TRUE, stringsAsFactors = FALSE)    #open csv file "RNSBiofilm_PCA" containing gene differential expression data of E. faecalis strain

X = RNS_Biofilm[, 1:2972]  
X_mean = apply(X, 2, mean)
X_center = scale(X, center = X_mean, scale = FALSE)
Y = RNS_Biofilm$Sample

pca.RNS = pca(X_center, ncomp = 6, center = FALSE, scale = FALSE)
plot(pca.RNS, ylim = c(0, 0.8), main = "PCA on RNSBiofilm")

ev = pca.RNS$prop_expl_var   #explained_variance
print(ev)
PCA_ev = data.frame(ev)
write.csv(PCA_ev, "RNS_Biofilm_PCA_ev.csv")  #export explained variance for each principal component into a csv file

sco = pca.RNS$variates   #prinicipal component values
print(sco)
sco_PC = data.frame(SampleNo = RNS_Biofilm$Sample, sco)
write.csv(sco_PC, "RNS_Biofilm_PCA_PC.csv")  #export the principal component values into a csv file

#Exported csv files will be used for further data analyses and drawing figures in Excel, Prism, and R programming
#End#
