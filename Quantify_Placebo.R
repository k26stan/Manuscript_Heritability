## Quantify Placebo Effect using Time-Series Data ##
## February 26, 2015 ##
## Kristopher Standish ##

library( nlme )

###########################################################
## LOAD DATA ##############################################
###########################################################

## Set Date
DATE <- gsub("-","",Sys.Date())

## Set Paths to Data and to Save
PathToFT <- "/Users/kstandis/Data/Burn/Data/Phenos/Full_Tables/20141229_Full_Table.txt"
PathToRep <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150226_Resp_v_Time.txt"
PathToSave <- paste("/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Plots/",DATE,"_Plac/",sep="")
dir.create( PathToSave )

## Load Real Data
FT <- read.table( PathToFT, sep="\t",header=T )
RP <- read.table( PathToRep, sep="\t",header=T )

## Game Plan
# Repeated Measures
  # Calculate effect of Placebo using all time points from RP table
  # Use simple Linear Regression
  # Test multiple different phenotypes
  # Assume I've already established transformation
  # Calculate using single timepoint as well (t-test)
  #*Calculate Placebo Effect vs Drug Effect
# XXXXXXX THIS WON'T WORK XXXXXXX # Full Table 
  # Use simplified delta-Phenotype calculation at single time-points
  # Use Cochran-Mantel-Hanzel test (via ClTrial) to calculate Plac effect
    # This might not work...will have to 

###########################################################
## QUANTIFY PLACEBO EFFECT ################################
###########################################################

## Calculate Transformations of Phenotypes
 # DAS/lCRP/rSJC/rTJC/rSJC28/rTJC28
PH <- data.frame( lCRP=log10(RP$CRP), rSJC=sqrt(RP$SJC), rTJC=sqrt(RP$TJC), rSJC28=sqrt(RP$SJC28), rTJC28=sqrt(RP$TJC28) )
DAT.1 <- data.frame( RP, PH )
WHICH_PHENOS <- c("DAS","lCRP","rSJC","rTJC")
PH_COLS <- which( colnames(DAT.1) %in% WHICH_PHENOS )# c(16,21:25)
names(PH_COLS) <- colnames(DAT.1)[PH_COLS]

## Remove timepoints w/ DRUG==1
DAT.2 <- DAT.1[ which(DAT.1$DRUG==0), ]

###########################################################
## Simple Linear Regression w/ PLAC as fixed effect
MODS <- COEFS <- list()
PS <- numeric( length(PH_COLS) )
names(PS) <- names(PH_COLS)
for ( c in 1:length(PH_COLS) ) {
	col <- PH_COLS[c]
	pheno <- names(PH_COLS)[c]
	TEMP_DAT <- DAT.2[, c(pheno,"PLAC") ]
	colnames(TEMP_DAT)[1] <- "Pheno"
	MOD <- lm( Pheno ~ PLAC, data=TEMP_DAT )
	COEF <- summary(MOD)$coefficients
	P <- COEF["PLAC","Pr(>|t|)"]
	# Compile
	MODS[[pheno]] <- MOD
	COEFS[[pheno]] <- COEF
	PS[c] <- P
}
LM <- list( MODS, COEFS, PS )
names(LM) <- c("Mods","Coefs","Ps")
unlist(lapply( COEFS, function(x) x[2,1]/x[1,1] ))

###########################################################
## Linear Mixed Models w/ PLAC as fixed effect & IID as Grouping Factor
MODS <- SUMMS <- COEFS <- INTS <- list()
PS <- numeric( length(PH_COLS) )
names(PS) <- names(PH_COLS)
for ( c in 1:length(PH_COLS) ) {
	col <- PH_COLS[c]
	pheno <- names(PH_COLS)[c]
	TEMP_DAT <- DAT.2[, c("IID",pheno,"PLAC") ]
	colnames(TEMP_DAT)[2] <- "Pheno"
	MOD <- lme( fixed = Pheno ~ PLAC, random = ~ 1 | IID, data=TEMP_DAT, subset=which(!is.na(TEMP_DAT$Pheno)) )
	SUMM <- summary(MOD)
	COEF <- SUMM$tTable
	# INT <- intervals(MOD)$fixed
	P <- COEF["PLAC","p-value"]
	# Compile
	MODS[[pheno]] <- MOD
	SUMMS[[pheno]] <- SUMM
	COEFS[[pheno]] <- COEF
	# INTS[[pheno]] <- INT
	PS[c] <- P
}
LME <- list( MODS, COEFS, PS )
names(LME) <- c("Mods","Coefs","Ps")
unlist(lapply( COEFS, function(x) x[2,1]/x[1,1] ))

###########################################################
## QUANTIFY PLACEBO vs DRUG EFFECT ########################
###########################################################

###########################################################
## Simple Linear Regression w/ PLAC as fixed effect
MODS <- COEFS <- list()
PS <- numeric( length(PH_COLS) )
names(PS) <- names(PH_COLS)
for ( c in 1:length(PH_COLS) ) {
	col <- PH_COLS[c]
	pheno <- names(PH_COLS)[c]
	TEMP_DAT <- DAT.1[, c(pheno,"PLAC","DRUG") ]
	colnames(TEMP_DAT)[1] <- "Pheno"
	MOD <- lm( Pheno ~ PLAC+DRUG, data=TEMP_DAT )
	COEF <- summary(MOD)$coefficients
	P.pl <- COEF["PLAC","Pr(>|t|)"]
	P.dr <- COEF["DRUG","Pr(>|t|)"]
	# Compile
	MODS[[pheno]] <- MOD
	COEFS[[pheno]] <- COEF
	PS[c] <- P.pl
}
LM.DP <- list( MODS, COEFS, PS )
names(LM.DP) <- c("Mods","Coefs","Ps")
unlist(lapply( COEFS, function(x) x[3,1]/x[2,1] ))

###########################################################
## Linear Mixed Models w/ PLAC as fixed effect & IID as Grouping Factor
MODS <- SUMMS <- COEFS <- INTS <- list()
PS <- numeric( length(PH_COLS) )
names(PS) <- names(PH_COLS)
for ( c in 1:length(PH_COLS) ) {
	col <- PH_COLS[c]
	pheno <- names(PH_COLS)[c]
	TEMP_DAT <- DAT.1[, c("IID",pheno,"PLAC","DRUG") ]
	colnames(TEMP_DAT)[2] <- "Pheno"
	MOD <- lme( fixed = Pheno ~ PLAC+DRUG, random = ~ 1 | IID, data=TEMP_DAT, subset=which(!is.na(TEMP_DAT$Pheno)) )
	SUMM <- summary(MOD)
	COEF <- SUMM$tTable
	# INT <- intervals(MOD)$fixed
	P.pl <- COEF["PLAC","p-value"]
	P.dr <- COEF["DRUG","p-value"]
	# Compile
	MODS[[pheno]] <- MOD
	SUMMS[[pheno]] <- SUMM
	COEFS[[pheno]] <- COEF
	# INTS[[pheno]] <- INT
	PS[c] <- P.pl
}
LME.DP <- list( MODS, COEFS, PS )
names(LME.DP) <- c("Mods","Coefs","Ps")
unlist(lapply( COEFS, function(x) x[3,1]/x[2,1] ))

###########################################################
## COMPILE & PLOT RESULTS #################################
###########################################################

## Compile Beta Ratios & P-Values
 # Beta
LM.B <- sapply( lapply( LM$Coefs, function(x) x["PLAC",1] ), "[", 1 )
LME.B <- sapply( lapply( LME$Coefs, function(x) x["PLAC",1] ), "[", 1 )
LM.DP.B <- sapply( lapply( LM.DP$Coefs, function(x) x["PLAC",1] ), "[", 1 )
LME.DP.B <- sapply( lapply( LME.DP$Coefs, function(x) x["PLAC",1] ), "[", 1 )
LM.DP.B.dr <- sapply( lapply( LM.DP$Coefs, function(x) x["DRUG",1] ), "[", 1 )
LME.DP.B.dr <- sapply( lapply( LME.DP$Coefs, function(x) x["DRUG",1] ), "[", 1 )
Bs <- rbind( LM.B, LME.B, LM.DP.B, LME.DP.B, LM.DP.B.dr, LME.DP.B.dr )
 # Confidence Intervals (for B)
LM.conf <- sapply( lapply( LM$Mods, function(x) confint(x)["PLAC",] ), "[", c(1,2) )
LME.conf <- sapply( lapply( LME$Mods, function(x) intervals(x)$fixed["PLAC",c(1,3)] ), "[", c(1,2) )
LM.DP.conf <- sapply( lapply( LM.DP$Mods, function(x) confint(x)["PLAC",] ), "[", c(1,2) )
LME.DP.conf <- sapply( lapply( LME.DP$Mods, function(x) intervals(x)$fixed["PLAC",c(1,3)] ), "[", c(1,2) )
LM.DP.conf.dr <- sapply( lapply( LM.DP$Mods, function(x) confint(x)["DRUG",] ), "[", c(1,2) )
LME.DP.conf.dr <- sapply( lapply( LME.DP$Mods, function(x) intervals(x)$fixed["DRUG",c(1,3)] ), "[", c(1,2) )
Confs <- rbind( LM.conf, LME.conf, LM.DP.conf, LME.DP.conf, LM.DP.conf.dr, LME.DP.conf.dr )
 # Beta Ratios
LM.Ratio <- sapply( lapply( LM$Coefs, function(x) x[2,1]/x[1,1] ), "[", 1 )
LME.Ratio <- sapply( lapply( LME$Coefs, function(x) x[2,1]/x[1,1] ), "[", 1 )
LM.DP.Ratio <- sapply( lapply( LM.DP$Coefs, function(x) x[2,1]/x[3,1] ), "[", 1 )
LME.DP.Ratio <- sapply( lapply( LME.DP$Coefs, function(x) x[2,1]/x[3,1] ), "[", 1 )
Ratios <- rbind( LM.Ratio, LME.Ratio, LM.DP.Ratio, LME.DP.Ratio )
 # P-Values
Pvals <- rbind( LM$Ps, LME$Ps, LM.DP$Ps, LME.DP$Ps )
Stars <- array(, dim(Pvals) )
Stars[ which(Pvals<.05) ] <- "*"
Stars[ which(Pvals<.01) ] <- "**"
Stars[ which(Pvals<.001) ] <- "***"

## Plot Beta values/ratios for each Phenotype
COLS <- c("tomato2","slateblue3")
png( paste(PathToSave,"Plac_Effect_Beta.png",sep="/"), height=1200,width=1600, pointsize=26 )
par(mfrow=c(2,1))
# B for PLAC beta & Intercept
YLIM <- c( 1.5*min(Bs[1:2,]), 0 )
X_BARS <- barplot( Bs[1:2,], col=COLS, beside=T, ylim=YLIM, xlab="Phenotype",ylab="Beta",main="Placebo Effect (Beta Value)" )
abline( h=seq(-5,0,.2), lty=2, col="grey50" )
barplot( Bs[1:2,], col=COLS, beside=T, add=T )
X.arrow <- c(t( X_BARS )) # rep( c(1.5,2.5), rep(6,2) ) + rep( seq(0,15,3),2 )
arrows( X.arrow, c(Confs[1,],Confs[3,]), X.arrow, c(Confs[2,],Confs[4,]), code=3,angle=90 )
text( X.arrow, c(Confs[1,],Confs[3,])-.1, label=rep("**",length(X.arrow)) )
# B for PLAC beta & DRUG beta
COLS <- c("tomato2","slateblue3","chartreuse2","deepskyblue2")
YLIM <- c( 1.5*min(Bs[3:6,]), 0 )
X_BARS <- barplot( Bs[3:6,], col=COLS, beside=T, ylim=YLIM, xlab="Phenotype",ylab="Beta",main="Placebo & Drug Effect (Beta Value)" )
abline( h=seq(-5,0,.5), lty=2, col="grey50" )
barplot( Bs[3:6,], col=COLS, beside=T, ylim=YLIM, add=T )
X.arrow <- c(t( X_BARS ))
arrows( X.arrow, c(Confs[5,],Confs[7,],Confs[9,],Confs[11,]), X.arrow, c(Confs[6,],Confs[8,],Confs[10,],Confs[12,]), code=3,angle=90 )
legend( "bottomright", fill=COLS, legend=c("LR: Plac","LMM: Plac","LR: Drug","LMM: Drug"), ncol=2 )
text( X.arrow, c(Confs[5,],Confs[7,],Confs[9,],Confs[11,])-.2, label=rep("**",length(X.arrow)) )
dev.off()

## Plot Beta values/ratios for each Phenotype
COLS <- c("tomato2","slateblue3")
png( paste(PathToSave,"Plac_Effect_Ratio.png",sep="/"), height=1200,width=1600, pointsize=26 )
par(mfrow=c(2,1))
# Ratio b/n PLAC beta & Intercept
YLIM <- c( 0, -1.2*min(Ratios[1:2,]) )
barplot( -Ratios[1:2,], col=COLS, beside=T, ylim=YLIM, xlab="Phenotype",ylab="-Beta(Placebo) / Beta(Intercept)",main="Placebo Effect (Proportion of Intercept)" )
abline( h=seq(0,1,.1), lty=2, col="grey50" )
barplot( -Ratios[1:2,], col=COLS, beside=T, add=T )
# Ratio b/n PLAC beta & DRUG beta
YLIM <- c( 0, 1.2*max(Ratios[3:4,]) )
barplot( Ratios[3:4,], col=COLS, beside=T, ylim=YLIM, xlab="Phenotype",ylab="Beta(Placebo) / Beta(Drug)",main="Placebo Effect (Proportion of Drug Effect)" )
abline( h=seq(0,1,.1), lty=2, col="grey50" )
barplot( Ratios[3:4,], col=COLS, beside=T, ylim=YLIM, add=T )
legend( "topleft", fill=COLS, legend=c("Lin. Regr.","LMM"), ncol=2 )
dev.off()

## Plot P-Values for each Phenotype
png( paste(PathToSave,"Plac_Effect_P.png",sep="/"), height=1200,width=1600, pointsize=26 )
COLS <- c("tomato2","slateblue3","tomato2","slateblue3")
PCHS <- c(1,1,4,4)
YLIM <- c( 0, ceiling(-log10(min(Pvals))) )
XLIM <- c( 1, ncol(Pvals) )
plot( 0,0,type="n", ylim=YLIM,xlim=XLIM, xlab="Phenotype",ylab="-log10(P)",main="Significance of Placebo Effect across Phenotypes",xaxt="n" )
axis( 1, at=1:ncol(Pvals), label=colnames(Pvals) )
abline( h=seq(0,50,5), lty=2,col="grey50" )
abline( v=seq(0,10,1), lty=2,col="grey50" )
for ( r in 1:nrow(Pvals) ) {
	points( 1:ncol(Pvals), -log10(Pvals[r,]), pch=PCHS[r],col=COLS[r],cex=2,lwd=4 )
}
legend( "bottomright", pch=PCHS,col=COLS,legend=c("LR: PLAC","LMM: PLAC","LR: PLAC+DRUG","LMM: PLAC+DRUG"), pt.cex=2,pt.lwd=4 )
dev.off()









###########################################################
## END OF DOC #############################################
###########################################################
