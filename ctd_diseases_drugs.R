##################################################################################################
#                                                                                                #
#                    Extract disease-related drugs from CTD,                                     #
#        use ID to add drugs with aliases for later verification (that is, as a standard)        #
#                                                                                                #
##################################################################################################

########################################################
#     Extract all drugs and IDs in CTD (only once, for subsequent use)
#######################################################
# rm(list = ls())
ctd_drugs = read.csv("E:\\work\\ctd\\CTD_chemicals.csv", header = T)
drugsA = tolower(as.character(ctd_drugs$ChemicalName))
drugB = tolower(as.character(ctd_drugs$Synonyms))
id = as.character(ctd_drugs$ChemicalID)
len = length(drugB)
all_drugs = c("")
drugs_id = c("")
k = 1
for(i in 1:len){
  # i = 1
  all_drugs[k] = drugsA[i]
  drugs_id[k] = id[i]
  k = k + 1
  sub = stringr::str_split(drugB[i], "[|]")
  for(j in 1:length(sub)){
    if(sub[[1]][j] != ""){
      all_drugs[k] = sub[[1]][j]
      drugs_id[k] = id[i]
      k = k + 1
    }
  }
}
result = data.frame(all_drugs, drugs_id)
result = unique(result)
print(length(unique(drugs_id)))
print(length(unique(all_drugs)))
write.csv(result, "E:\\work\\ctd\\ctd_alldrugs_id.csv", row.names = F, col.names = F)

###################################################################
#
#  Use the ID to find out the medicine related to the disease
#
###################################################################

result = read.csv("D:\\work\\ctd_alldrugs_id.csv", header = T)
breast_data = read.csv("merge.csv",header = T)
id = as.character(breast_data$Chemical.ID)
id = unique(id)
result_breast = result[1,]
for(i in 1:length(id)){
  pos = which(result$drugs_id == id[i])
  if(length(pos) > 0){
    result_breast = rbind(result_breast, result[pos,])
  }
}
result_breast = result_breast[-c(1),]
write.csv(result_breast,"ctd_drugs.csv",row.names = F, col.names = F)
