##########################################################################
#                                                                        #                                
#    Differential expression of the data downloaded by TCGA(type:counts) #
#                   download from fire..（breast）                       #
#                                                                        #
#######################################################################

##########  Extract the counts in the expression data
# rm(list = ls())
setwd("E:/further/8A network/TCGA/lung")

data = read.table("breast.txt", sep = "\t", header = T)
data = as.matrix(data)
type = as.character(data[1,])
countdata = data[,which(type == "raw_counts")]
result = cbind(data[,1], countdata)
write.csv(result, "breast_counts.csv", col.names = F, row.names = F)

##### Put the above files in the order of control case
data = read.csv("breast_counts.csv", header = F)
data = as.matrix(data)
samples = as.character(data[1,-c(1)])
len = length(samples)
control = c(0)
k = 1
case = c(0)
w = 0
for(i in 1:len){
  # i= 1
  substr = strsplit(samples[i], "-")
  state = substr[[1]][4]
  if(state < "10"){
    case[w] = i
    w = w + 1
  }
  if(state >= "10")
  {
    control[k] = i
    k = k + 1
  }
}
result = data[, 1]
result = cbind(result, data[,control+1], data[,case+1])
write.csv(result, "breast_control_case.csv", col.names = F, row.names = F)

######### Do differential expression
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
install.packages("gplots")
foldChange=1.5
padj=0.05

library("edgeR")
library("limma")

rt = read.csv("lung_control_case.csv", header = T)
rt = as.matrix(rt)
rownames(rt) = rt[,1]
exp = rt[,2:ncol(rt)]
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data = avereps(data)  #Genes appearing multiple times are averaged
data = data[rowMeans(data) > 1,] #Remove low-expression data

group = c(rep("normal", 37), rep("tumor", 124))
design = model.matrix(~group)

y = DGEList(counts = data, group = group) #Which are normal and which are cancer, so that edgeR can recognize
y = calcNormFactors(y)   #Factor correction
y = estimateCommonDisp(y)  #Estimate the coefficient of variation, that is, estimate the variance; estimate the degree of internal difference to see whether the difference between groups is greater than the internal difference. If it is large, select the differential gene
y = estimateTagwiseDisp(y)

et = exactTest(y, pair = c("normal", "tumor"))
topTags(et)
ordered_tags = topTags(et, n = 100000)  #Show top 100,000

allDiff = ordered_tags$table
allDiff = allDiff[is.na(allDiff$FDR) == FALSE,]
diff = allDiff
newData = y$pseudo.counts
###all genes
write.csv(diff, "edgerOut.csv", quote = F) #First express all the differences and output them to the edgerOut file
###Significant gene
diffSig = diff[(diff$FDR < padj & (diff$logFC > foldChange | diff$logFC < (-foldChange))),]   #Screen all significant differences
write.csv(diffSig, "diffSig.csv", quote = F) #Significantly different output
###Up- down- gene
diffUp = diff[(diff$FDR < padj & (diff$logFC > foldChange)),] 
write.csv(diffUp, "up.csv", quote = F)
diffDown = diff[(diff$FDR < padj & (diff$logFC < (-foldChange))),] 
write.csv(diffDown, "down.csv", quote = F)

