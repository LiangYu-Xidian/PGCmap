##################################################################################################
#                                                                                                #               
#     calculate ks_score for lincs drug(size = 1w+ genes)  [drugs use specific cell line]        #       
#                   disease : TCGA diff_sig   drug: LINCS(diff idoses)                           #
#                                                                                                #           
##################################################################################################

# rm(list = ls())


# function to calculate ks_score
cmap_score <- function(sig_up, sig_down, drug_signature) {
  # calculate interaction score between two sets of top-ranked genes and a
  # drug-induced expression signature.
  # Inputs:
  #   sig_up - numeric of GeneIDs from top end of ranked list (lo expression)
  #   sig_down - numeric GeneIDs from bottom end of ranked list (hi expression)
  #   drug_signature - data.frame GeneIDs and rank
  # Returns:
  #   connectivity score
  
  num_genes <- nrow(drug_signature)
  ks_up <- 0
  ks_down <- 0
  connectivity_score <- 0
  
  # re-rank because the GeneID mapping changed the original rank range
  drug_signature[, "rank"] <- rank(drug_signature[, "value"])           # Ascending order
  
  # Merge the drug signature with the disease signature by GeneID. This becomes
  # the V(j) from the algorithm description.
  up_tags_rank <- merge(drug_signature, sig_up, by.x="GeneID", by.y=1)   
  down_tags_rank <- merge(drug_signature, sig_down, by.x="GeneID", by.y=1)
  
  up_tags_position <- sort(up_tags_rank$rank)
  down_tags_position <- sort(down_tags_rank$rank)

  num_tags_up <- length(up_tags_position)
  num_tags_down <- length(down_tags_position)
  
  if(num_tags_up > 1 && num_tags_down > 1) {
    a_up <- 0
    b_up <- 0
    
    # small speed up by changing from sapply to inline math (~0.5 sec)
    a_up <- max(((1:num_tags_up) / num_tags_up) - (up_tags_position / num_genes))
    b_up <- max((up_tags_position / num_genes) - (((1:num_tags_up)-1) / num_tags_up))
    
    if(a_up > b_up) {
      ks_up <- a_up
    } else {
      ks_up <- -b_up
    }
    
    a_down <- 0
    b_down <- 0
    
    # small speed up by changing from sapply to inline math (~0.5 sec)
    a_down <- max(((1:num_tags_down) / num_tags_down) - (down_tags_position / num_genes))
    b_down <- max((down_tags_position / num_genes) - (((1:num_tags_down)-1) / num_tags_down))
    
    if(a_down > b_down) {
      ks_down <- a_down
    } else {
      ks_down <- -b_down
    }
    
    if (sum(sign(c(ks_down, ks_up))) == 0) {  
      connectivity_score <- ks_up - ks_down # different signs
    }
  }
  return(connectivity_score)
}

######## read data ##############################################################
######## read disease data ########
MS <- read.csv("D:\\work\\breast\\nodelete\\input\\diffSig_id.csv", header = TRUE)

### use for calculate score
### sort(large->small)
d <- MS[order(MS$logFC),]
geneid <- d$ENTREZID
zscore <- d$logFC
ZData1 <- data.frame(geneid,zscore)
colnames(ZData1) <- c("GeneID","value")
### delete NA
ZData1 = na.omit(ZData1)

print(length(which(ZData1$value < 0)))
print(length(which(ZData1$value > 0)))


#### drug`s files
setwd("D:\\work\\breast\\lincs")
temp <- list.files(pattern = '*.csv')

# the parameters of ks, respectively for top size and bottom size
 top_size = c(100,200,300,400,500,600)
 bottom_size = c(100,200,300,400,500,600)

for(size in 1:length(top_size)){
  for(drug_size in 1:length(temp)){
    CMapData <- read.csv(temp[drug_size], header = FALSE)
    CMapData <- CMapData[,-c(2)]  # delete DMSO
    head <- CMapData[c(1:7),]
    CMapData <- CMapData[-c(1:7),]
    colnames(CMapData) = c("GeneID")
    ### pull out matching_genes  
    matching_genes = unique(merge(CMapData, ZData, by.x = "GeneID")[,1])
    print(length(matching_genes)) 
  
	# number of columns in each drug file
    N = ncol(CMapData) - 1
    score=numeric(N)
    score_pvalue=numeric(N)
    
    ### for each drug with different experimental conditions
    for(i in 1:N){
      # i=1
      drug_signature = CMapData[,c(1, i+1)]
      colnames(drug_signature) = c("GeneID","value")
      drug_signature = as.data.frame(drug_signature)
      
      ### select top and bottom of rank and collect enrichment score
      ms_up = ZData1$GeneID[1:top_size[size]]
      gsel = (nrow(ZData1) - bottom_size[size] + 1):nrow(ZData1)
      ms_down = ZData1$GeneID[gsel]
      
      score[i] = cmap_score(ms_up, ms_down, drug_signature)
      
    } # end for
    ### write score 
    file1 = strsplit(temp[drug_size],'.csv')
    ofile = paste("D:\\work\\breast\\nodelete\\mid\\D_",top_size[size],"\\", file1[[1]], "_score.csv", sep = "")
    write.csv(score, file = ofile, quote = F, col.names = F, row.names = F)
    ############ merge and add header information  ###############
     score = t(score) 
    rownames(score) = c("score")
    name = as.character(head[,1])
    name[8] = "score" 
    head <- head[,-c(1)]
    head = as.matrix(head)
    score_fdr = as.matrix(score)
    data <- rbind(head,score_fdr)
    data <- cbind(name,data)
    
    ofile = paste("D:\\work\\breast\\nodelete\\mid\\D_",top_size[size],"\\", file1[[1]], "_head_score.csv", sep = "")
    write.csv(data, ofile,  col.names = F, row.names = F)
  }
}



