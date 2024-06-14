#Sparse partial least squares discriminant analysis (sPLSDA) on untargeted LC-MS profiling of Enterococcus faecalis exometabolome that is affected by Klebsiella pneumoniae in CAUTI polymicrobial community.
#Beign#

library(knitr)
library(mixOmics)
library(snow)

MS_EF = read.csv(file = file.choose(), sep = ",", header = TRUE, stringsAsFactors = FALSE) #open the file "ZZJAW_EMS230706_3_ForAna.csv"
X = MS_EF[, 1:2536]
X_mean = apply(X, 2, mean)
X_center = scale(X, center = X_mean, scale = FALSE)
Y1 = MS_EF$Group1   #Fraction groups
splsda.MS = splsda(X, Y1, ncomp = 10)

ev = explained_variance(splsda.MS$X, splsda.MS$variates$X, ncomp = 10)  #explained variance
write.csv(ev, "MS_sPLSDA_ExplainedVariance.csv")  #export explained variance for each principal component into a csv file

sco = splsda.MS$variates   #prinicipal component
write.csv(sco, "MS_sPLSDA_PrincipleComponent.csv")  #export the principal component values into a csv file

loa = data.frame(splsda.MS$loadings$X)   #loading values
peak_name = colnames(X)   #peak names
loa_PC1 = data.frame(peak = peak_name, PC1loa = loa[, 1])   #create a dataframe with peak names, and their loading values associated with PC1
write.csv(loa_PC1, "SPLSDA_20230707TPSecC18LCMS_loa_PC1.csv")  #export loading values associated with PC1 into a csv file
loa_PC2 = data.frame(peak = peak_name, PC2loa = loa[, 2])   #create a dataframe with peak names, and their loading values associated with PC2
write.csv(loa_PC2, "SPLSDA_20230707TPSecC18LCMS_loa_PC2.csv")  #export loading values associated with PC2 into a csv file

#Exported csv files will be used for further data analyses and drawing figures in Excel, Prism, and R programming
#End#
