# This script conducts a polygenic risk score (PRS) analysis using PLINK output files and phenotype data.
# It first reads in the phenotype file, principal components (PCs), and covariates. 
# Then, it calculates the null model R-squared and iterates through different p-value thresholds to calculate PRS R-squared,
# coefficients, standard errors, and p-values. Finally, it identifies the best PRS result and saves the results to output files.

p.threshold <- c(0.001,0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5)

# Read the phenotype file 
phenotype <- read.table( "/data/users/itreccani/ADNI/PRS/Target.Final.pruned.fam", header = F)
phen <- data.frame ( phenotype$V1, phenotype$V2, phenotype$V6)
colnames(phen) <- c( "FID",	"IID",	"Phenotype" )
write.table( phen, "/data/users/itreccani/PRS/phenotype.phen", row.names = F, col.names = T, quote = F, sep = "\t") 

# Read the PCs
pcs <- read.table( "/data/users/itreccani/ADNI/PCA.eigenvec", header = T)
pcs <- pcs[,1:4]
colnames(pcs) <- c( "FID", "IID", paste0( "PC", 1:2)) 

# Read the covariate
covariate <- read.table( "/data/users/itreccani/ADNI/PRS/Target.Final.pruned.fam", header=F)
cov<-data.frame( covariate$V1, covariate$V2, covariate$V5)
colnames(cov)<-c("FID",	"IID",	"Sex")

write.table( cov, "/data/users/itreccani/ADNI/PRS/sex.sex", row.names = F, col.names = F, quote = F, sep = "\t") 
pheno <- merge( merge( phen, cov, by = c( "FID", "IID" )), pcs, by = c( "FID", "IID" ))
null.model <- lm( Phenotype~., data = pheno[ , !colnames( pheno ) %in% c( "FID", "IID" )])
null.r2 <- summary(null.model)$r.squared
prs.result <- NULL

for(i in p.threshold){
  
    # Go through each p-value threshold
    prs <- read.table( paste0("/data/users/itreccani/ADNI/PRS/prs.plink.",i,".profile"), header=T)
    pheno.prs <- merge( pheno, prs[ , c("FID","IID", "SCORE")], by = c("FID", "IID"))
    model <- lm( Phenotype~., data = pheno.prs[ , !colnames( pheno.prs ) %in% c("FID","IID")])
    model.r2 <- summary(model)$r.squared
    
    # R2 of PRS is simply calculated as the model R2 minus the null R2
    prs.r2 <- model.r2-null.r2
    
    # Obtain the coeffcient and p-value of association of PRS
    prs.coef <- summary(model)$coeff["SCORE", ]
    prs.beta <- as.numeric(prs.coef[1])
    prs.se <- as.numeric(prs.coef[2])
    prs.p <- as.numeric(prs.coef[4])
    prs.result <- rbind(prs.result, data.frame( Threshold = i, R2 = prs.r2, P = prs.p, BETA = prs.beta, SE = prs.se))
}

# Best result is:
prs.result[ which.max( prs.result$R2 ), ]
write.table( prs.result[ which.max ( prs.result$R2 ), ], "/data/users/itreccani/ADNI/PRS/result.txt",row.names = F, col.names = T, quote = F, sep = "\t") 
write.table( prs.result,  "/data/users/itreccani/ADNI/PRS/plink.prs.result.txt", row.names = F, col.names = T, quote = F,sep ="\t") 