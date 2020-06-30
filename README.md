# PGCmap:Cmap method considering Pathogenic Genes

## Dependencies
*MEDTI* is tested to work under R 3.5

## Code and Data
#### Code
- `calculateLOGFC.R`: Calculate the differential expression of logFC for disease data(diseases data download form TCGA)
- `calculateScore.R`: Calculate the correlation score between disease and medicine
- `ctd_diseases_drugs.R`: Extract disease-related drugs from CTD
- `drugRank_precision.R`: Calculate the accuracy of predicted drugs


#### Data: `test_data/` directory
- `diseases differential expression.xlsx` 	    	: Gene differential expression data for multiple diseases
- `A375.csv`       		: Differential expression of drugs in LINCS database under A375 cell line
- `HA1E.csv` 			: Differential expression of drugs in LINCS database under HA1E cell line
- `HT29.csv` 			: Differential expression of drugs in LINCS database under HT29 cell line
**Note**: gene id all use ENTREZ ID


### Tutorial
1. Use `calculateLOGFC.R` to calculate the differential expression value logFC of the genes in the disease data downloaded from TCGA, and calculate the corresponding significant FDR, and screen the genes according to the conditions.
2. The differential expression value of the gene under the disease and the order of the gene under the drug are taken as inputs to calculate the correlation score between the drug and the disease.
3. Download disease-related drug data from the CTD database as a verification standard to calculate and predict the accuracy of the drug.
4. Run `drugRank_precision.R`, Using the drug data of xx as a standard, the accuracy of segment prediction is calculated for every multiple drugs.

### Contacts
If you have any questions or comments, please feel free to email Liang Yu (lyu@xidian.edu.cn) and/or Dan He (hedan@163.com).