#####################################
## rank
# rm(list = ls())
 size = c(100,200,300,400,500,600)
for(s in size){
  # s = 96
  f = paste("D:\\work\\breast\\nodelete\\mid\\all_score_sort\\all_score_", s, ".csv", sep = "")
  data = read.csv(f, header = F)
  data = data[-c(1),-c(1)]
#  data = data[,-c(1)]
  data = as.matrix(data)
  len = length(data[8,])
  score = c(0)
  for(i in 1:len){
    score[i] = as.numeric(as.character(data[8,i]))
  }
  data1 = data[,order(-abs(score))]
  f = paste("D:\\work\\breast\\nodelete\\mid\\sort_drugs\\all_normal_abs_sort_", s,".csv", sep = "")
  write.csv(data1, f,row.names = F, col.names = F)
  data1 = data[,order(score)]
  data1 = data1[,which(data1[8,] < 0)]
  print(length(data1[3,]))
  f = paste("D:\\work\\breast\\nodelete\\mid\\sort_drugs\\all_normal_less_sort_", s,".csv", sep = "")
  write.csv(data1, f,row.names = F, col.names = F)
  
}
########### drugs rank #############
# rm(list = ls())
 size = c(100,200,300,400,500,600)
for(s in size){
  # s = 1000
  f = paste("D:\\work\\breast\\nodelete\\mid\\sort_drugs\\all_normal_abs_sort_", s,".csv", sep = "")
  data = read.csv(f, header = F)
  data = data[-c(1),]
  data = as.matrix(data)
  ## all drugs 
  drugs = as.character(data[4,])
  drugs = unique(drugs)
  print(length(drugs))
  
  result = data
  drugs = c("")
  for(i in 1:length(result[4,])){
    drugs[i] = as.character(result[4,i])
  }
  drugs = unique(drugs)
  drugs = tolower(drugs)
  print(length(drugs))
  of = paste("D:\\work\\breast\\nodelete\\mid\\sort_drugs\\all_score_", s, "_abs_sort_drugs.csv", sep = "")

  write.csv(drugs, of, col.names = F, row.names = F)
  
}


# rm(list = ls())
#### ctd 
ctd_data = read.csv("D:\\work\\ctd_alldrugs_id.csv", header = T)
ctd_drugs = as.character(ctd_data$all_drugs)
### breast_ctd
breast_data = read.csv("D:\\work\\breast\\ctd\\ctd_drugs.csv", header = T)
breast_drugs = as.character(breast_data$all_drugs)
print(length(unique(breast_drugs)))
print(length(intersect(ctd_drugs, breast_drugs)))
### pre
 size = c(100,200,300,400,500,600)
for(s in size){
  # s = 400
  f = paste("D:\\work\\breast\\result\\nodelete\\mid\\sort_drugs\\all_score_", s, "_abs_sort_drugs.csv", sep = "")
  pre_data = read.csv(f, header = T)
  pre_drug = tolower(as.character(pre_data$x))
  ## Remove medicines not in CTD
  common = intersect(pre_drug, ctd_drugs)
  common_cure = intersect(common, breast_drugs)
  print("Number of drugs in CTD:intersect(pre_drug, ctd_drugs)£º")
  print(length(common))
  print("Number of drugs related to disease in CTD:intersect(common, breast_drugs)£º")
  print(length(common_cure))
  ## precision
  top_num = seq(10, length(common), by = 10)
  precision = rep(0, length(top_num))
  for(i in 1:length(top_num)){
    top_drug = common[c(1:top_num[i])]
    ctd_num = length(intersect(top_drug, breast_drugs))
    precision[i] = ctd_num / top_num[i]
  }
  precision
  
  par(lty=2, pch=15)
  main = paste(s," abs score", sep = "")
  plot(top_num, precision,type="b",col="black", ylim = c(0,1),xlab = "drug number M",ylab = "precision", main=main)
  
}




