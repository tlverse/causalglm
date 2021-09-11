 
load("beta.rcp.Rdata") # This loads the beta matrix and maybe also a "sheet" descriptor file
# This should be the beta/cpg matrix as stored in the RData file.
beta.rcp <- beta.rcp
dim(beta.rcp)

library(data.table) #install.packages("data.table")
dat_cpg <- data.table(CpGName = rownames((beta.rcp)),  (beta.rcp)) #Store the cpg matrix as a data.table with the first column being the cpg names/ids (usually the rownames but maybe not always)

# Download the necessary horvath cpgs from the calculator site
horvathCPG <- data.table(CpGName = read.csv("https://dnamage.genetics.ucla.edu/sites/all/files/tutorials/datMiniAnnotation3.csv")$Name)

# This "trick" will subset dat_cpg to only rows/cpgs that appear in horvathCPG and it will add "NA" values for missing cpgs as required for the horvath clock calculator
horvath_data <- dat_cpg[horvathCPG, on = "CpGName"]
write.csv(horvath_data, "horvathcalculator_beta_input.csv") # To be fed into online calculator
 




pheno <- read.csv("pheno.csv")
# pheno is a phenotype file with ids and sheet is a file that maps the phenotype ids to the cpg data ids. The sheet file may already be merged with the pheno file for you. In which case, this can be simplified.
pheno <- merge(pheno, sheet, by.x = "id", by.y = "id")  # May not be necessary
keep_cols <- intersect(colnames(horvath_data), pheno$Basename) # Find common observations between pheno and beta mat
horvath_data <- horvath_data[,c("CpGName", keep_cols), with = F] # If there are extra observations not found in both then remove them.
pheno <- pheno[na.omit(match(keep_cols, pheno$Basename)), ] # Make sure row order of phenofile matches the column order of beta/cpg matrix (check ids)
pheno$Female <- as.numeric(pheno$sex == 2)  # The pheno file should have a column called exactly "Female" that is 1 for "female" and 0 for "male".
pheno$Age <- pheno$age_9y # The pheno file should have a column called exactly "Age" with the chronological ages
write.csv(horvath_data, "horvathcalculator_pheno_input.csv")# To be fed into online calculator
 