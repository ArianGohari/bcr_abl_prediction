# Install all of the needed dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("XML", repos = "http://www.omegahat.net/R")
BiocManager::install("ALL")
BiocManager::install("MLInterfaces", force=TRUE)

library("ALL")
library("MLInterfaces")

# Load dataset
data("ALL")

# Get subset of only B-Cell ALL data
bALL = ALL[, substr(ALL$BT, 1, 1) == "B"]
fus = bALL[, bALL$mol.biol %in% c("BCR/ABL", "NEG")]
fus$mol.biol = factor(fus$mol.biol)
mads = apply(exprs(fus), 1, mad)
fusk = fus[mads > sort(mads, decr=TRUE)[300],]

# Visualize data using a heatmap
fcol = ifelse(fusk$mol.biol == "NEG", "green", "red")
heatmap(exprs(fusk), ColSideColors = fcol)

# Principal component analysis
PCg = prcomp(t(exprs(fusk)))
plot(PCg)

# Analyse principal components using biplot
pairs(PCg$x[, 1:5], col=fcol, pch=19)
biplot(PCg)

# Train a neural network on the dataset using MLearn
# A helper function for feature selection in cross-validation is used, 
# it picks the top 30 features ranked by absolute t statistic for 
# each cross-validation partition
nnALLFS = MLearn(mol.biol~., fusk, nnetI, xvalSpec("LOG", 5, balKfold.xvspec(5), fs.absT(30)), size=5, decay=.01, MaxNWts=2000)
nnALLFS
nn_eval <- confuMat(nnALLFS)

# Calculate precision, recall and accuracy
prec <- nn_eval[1, 1] / (nn_eval[1, 1] + nn_eval[2, 1])
rec <- nn_eval[1, 1] / (nn_eval[1, 1] + nn_eval[1, 2])
accuracy <- (nn_eval[1, 1] + nn_eval[2, 2]) / sum(nn_eval)
