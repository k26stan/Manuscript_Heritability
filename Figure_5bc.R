## Make Figure 5 (Derived Phenotypes) for Resp_Herit Manuscript ##
## July 8, 2015 ##
## Kristopher Standish ##

## Gameplan
 # Show discordance b/n ranks of patients at different timepoints
 # Show correlations b/n phenotypes
 # Show variance around mean/fit (after treatment?)
 # SIMULATIONS???

library( nlme )
library( gplots )
library( lmtest )

###########################################################
## LOAD DATA ##############################################
###########################################################

## Set Date
DATE <- gsub("-","",Sys.Date())

## Set Paths to Data and to Save
PathToDER.313 <- "/Users/kstandis/Data/Burn/Data/Phenos/Full_Tables/20150313_Derived_Pheno_Table.txt"
PathToFT <- "/Users/kstandis/Data/Burn/Data/Phenos/Full_Tables/20150520_Full_Table.txt"
PathToRep <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150530_Resp_v_Time.txt"
PathToWAG <- "/Users/kstandis/Data/Burn/Data/Phenos/Full_Tables/20150520_Single_Pheno_Table.txt"
PathToDER <- "/Users/kstandis/Data/Burn/Data/Phenos/Full_Tables/20150619_Derived_Pheno_Table.txt"
PathToSave <- paste("/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Plots/",DATE,"_Derived",sep="")
dir.create( PathToSave )

## Load Real Data
FT <- read.table( PathToFT, sep="\t",header=T )
RP <- read.table( PathToRep, sep="\t",header=T )
WAG <- read.table( PathToWAG, sep="\t",header=T )
DER <- read.table( PathToDER, sep="\t",header=T )
DER.313 <- read.table( PathToDER.313, sep="\t",header=T )

## Specify Phenotypes of Interest
WHICH_PHENOS <- c("DAS","lCRP","rSJC","rTJC")
COLS <- c("firebrick1","gold1","chartreuse1","dodgerblue1")
COLS.4 <- gsub("1","4",COLS)

###########################################################
## ORGANIZE DATA ##########################################
###########################################################

## Merge Full Table w/ Derived Phenotypes
MG <- merge( FT, DER, by.x="ID_2",by.y="FID" )

## Get set of Weeks
WKS <- as.numeric( unique( RP$WK ) )

## Set number of patients (before Filtering)
Samps.1 <- as.character( MG$ID_2 )
N.samps.1 <- length( Samps.1 )

## Pull weeks participating in trial for each patient
IN <- MG$IN
BL <- c(0,24,16)[factor(MG$GRP)]
names(IN) <- names(BL) <- MG$ID_2

## Remove Patients w/ little or not Treatment Data
RM.drop <- which( IN < 8 )
RM.drop.names <- names(RM.drop)

## Set number of patients
Samps <- setdiff( Samps.1, RM.drop.names )
N.samps <- length( Samps )

###########################################################
## CALCULATE DIFFERENCES/RESIDS ###########################
###########################################################
 # 4, 12, 20, 28 WAG
 # First & Last

###########################################################
## Distributions: Delta

## Create Tables of Delta-Values
# DEL_CATS <- c("MNw","MNwo","MNa","MNcd","Bdr","PRC","Bwk","VARdr","VARwk")
DEL_CATS <- c("MNa","MNcd","PRC","Bwk","VARwk")
DEL <- list()
DEL$DAS <- data.frame( MG[, paste("DEL",DEL_CATS,"DAS",sep="_") ] ) # data.frame( MG[, colnames(DER)[grep("DEL.*DAS",colnames(DER))] ] )
DEL$lCRP <- data.frame( MG[, paste("DEL",DEL_CATS,"lCRP",sep="_") ] )
DEL$rSJC <- data.frame( MG[, paste("DEL",DEL_CATS,"rSJC",sep="_") ] )
DEL$rTJC <- data.frame( MG[, paste("DEL",DEL_CATS,"rTJC",sep="_") ] )
# DEL$lCRP <- data.frame( MG[, colnames(DER)[grep("DEL.*CRP",colnames(DER))] ] )
# DEL$rSJC <- data.frame( MG[, colnames(DER)[setdiff( grep("DEL.*SJC",colnames(DER)),grep("28",colnames(DER)) ) ] ] )
# DEL$rTJC <- data.frame( MG[, colnames(DER)[setdiff( grep("DEL.*TJC",colnames(DER)),grep("28",colnames(DER)) ) ] ] )
for ( i in 1:length(DEL) ) { colnames(DEL[[i]]) <- DEL_CATS }

## Create Tables of Post-Treatment Values
POST_CATS <- c("MNw","MNwo","MNa","Bdr")
POST <- list()
POST$DAS <- data.frame( MG[, paste("POST",POST_CATS,"DAS",sep="_") ] ) # data.frame( MG[, colnames(DER)[grep("POST.*DAS",colnames(DER))] ] )
POST$lCRP <- data.frame( MG[, paste("POST",POST_CATS,"lCRP",sep="_") ] )
POST$rSJC <- data.frame( MG[, paste("POST",POST_CATS,"rSJC",sep="_") ] )
POST$rTJC <- data.frame( MG[, paste("POST",POST_CATS,"rTJC",sep="_") ] )
# POST <- list()
# POST$DAS <- data.frame( MG[, colnames(DER)[grep("POST.*DAS",colnames(DER))] ] )
# POST$lCRP <- data.frame( MG[, colnames(DER)[grep("POST.*CRP",colnames(DER))] ] )
# POST$rSJC <- data.frame( MG[, colnames(DER)[setdiff( grep("POST.*SJC",colnames(DER)),grep("28",colnames(DER)) ) ] ] )
# POST$rTJC <- data.frame( MG[, colnames(DER)[setdiff( grep("POST.*TJC",colnames(DER)),grep("28",colnames(DER)) ) ] ] )
for ( i in 1:length(POST) ) { colnames(POST[[i]]) <- POST_CATS }

###########################################################
## Distributions: Residual vs Initial Values

## Create table to write Pheno/Covs to (for GCTA)
WRITE_PHENO_TAB <- array( ,c(0,2) )

## Categories of Pheno using Pre- as Covariate
# RES_CATS <- c("MNw","MNwo","MNa","Bdr")
RES_CATS <- "MNa"

## Calculate Residuals vs Initial Values
RES <- DEL
BP.test <- array( ,c( length(RES_CATS), length(RES) ) )
colnames(BP.test) <- names(RES)
rownames(BP.test) <- RES_CATS
for ( p in 1:length(RES) ) {
	pheno <- names(RES)[p]
	for ( c in 1:length(RES_CATS) ) {
		DEL_COLNAME.w <- paste( "DEL",RES_CATS[c],pheno,sep="_" )
		DEL_COLNAME <- RES_CATS[c]
		PRE_COLNAME <- paste( "PRE",RES_CATS[c],pheno,sep="_" )
		MOD <- lm( DEL[[p]][,DEL_COLNAME] ~ MG[,PRE_COLNAME] )
		which_rows <- as.numeric( names(resid(MOD)) )
		RES[[p]][which_rows,DEL_COLNAME] <- resid( MOD )
		## BP Test
		BP <- bptest( MOD )$p.value
		BP.test[c,p] <- BP
		WRITE_PHENO_TAB <- rbind( WRITE_PHENO_TAB, c( DEL_COLNAME.w, PRE_COLNAME ) )
	}
	RES[[p]] <- data.frame( RES[[p]][, RES_CATS] ) # grep( paste(RES_CATS,collapse="|"), colnames(RES$DAS) ) ]
	colnames(RES[[p]]) <- RES_CATS
}

###########################################################
## FCT: PLOTS #############################################
###########################################################
COLS <- c("firebrick1","gold1","chartreuse1","dodgerblue1")

## Plot Function
SHIZ <- function(LIST, COLNAME, FILENAME, BIN_SIZES, LABEL) {
	P_VALS <- numeric(length(LIST))
	png( paste(PathToSave,"/4b_Der_Dist-",FILENAME,"_",COLNAME,".png",sep=""), height=700, width=2400, pointsize=30 )
	par(mfrow=c(1,4))
	# Loop Thru Phenotypes
	for ( i in 1:length(LIST) ) {
		# Pheno Name
		PHENONAME <- names(LIST)[i]
		XLIM <- c( floor(min(LIST[[i]][,COLNAME],na.rm=T)), max(LIST[[i]][,COLNAME],na.rm=T) )
		BIN <- BIN_SIZES[i]
		BRKS <- seq( XLIM[1], XLIM[2]+BIN, BIN )
		HIST <- hist( LIST[[i]][,COLNAME], breaks=BRKS, plot=F )
		YLIM <- c( 0, max(HIST$counts) )
		hist( LIST[[i]][,COLNAME], main=paste(LABEL,PHENONAME,"-",COLNAME), xlab=paste(LABEL,PHENONAME,sep="-"), col=COLS[i], xlim=XLIM, breaks=BRKS )
		abline( h=seq(0,100,10), lty=2,col="grey50" )
		hist( LIST[[i]][,COLNAME], main=paste(LABEL,PHENONAME,"-",COLNAME), xlab=paste(LABEL,PHENONAME,sep="-"), col=COLS[i], xlim=XLIM, breaks=BRKS, add=T )			
		P_VAL <- shapiro.test( LIST[[i]][,COLNAME] )$p.value
		text( XLIM[1], quantile(YLIM,.9), label=paste("Shapiro: P=",formatC(P_VAL,format="e",digits=2),sep=""), pos=4 )
		P_VALS[i] <- P_VAL
	}
	dev.off()
	return(P_VALS)
}
# SHIZ( DEL, "WAG_20","1_Delta", c(.25,2,1,1) )
# SHIZ( DEL.t, "WAG_20","1_Delta_t", c(.25,.25,.5,.5) )

###########################################################
## PLOT DISTRIBUTIONS/FITS ################################
###########################################################

SHAP.p <- list()

###########################################################
## Distributions: Delta
SHAP.p[["1_Delta"]] <- array( , c( length(DEL_CATS),length(DEL) ) )
rownames(SHAP.p[["1_Delta"]]) <- DEL_CATS
colnames(SHAP.p[["1_Delta"]]) <- names(DEL)
## Plot it
for ( del in DEL_CATS ) {
	SHAP.p[["1_Delta"]][del,] <- SHIZ( DEL, del,"1_Delta", c(.25,.25,.5,.5), "Delta" )
}

###########################################################
## Distributions: Residuals
SHAP.p[["2_Resid"]] <- array( , c( length(RES_CATS),length(RES) ) )
rownames(SHAP.p[["2_Resid"]]) <- RES_CATS
colnames(SHAP.p[["2_Resid"]]) <- names(RES)
## Plot it
for ( del in RES_CATS ) {
	SHAP.p[["2_Resid"]][del,] <- SHIZ( RES, del,"2_Resid", c(.25,.25,.5,.5), "Residuals" )
}

###########################################################
## PLOT P-VALUES OF ASSUMPTION TESTS ######################
###########################################################

## Plot Results
# P-Values for Shapiro Test for Original and Transformed Data
COLS <- c("firebrick1","gold1","chartreuse1","dodgerblue1")
PCHS <- c(0:4,7)
png( paste(PathToSave,"/4b_Der_Shapiro_Results.png",sep=""), height=1200, width=1200, pointsize=36 )
 # Original
XLIM <- c(1,ncol(SHAP.p[[1]]) ) + c(-.4,.4)
YLIM <- c(0, -log10(Reduce( min, SHAP.p[[1]] )) )
plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, main="Shapiro Test: Derived Stats",xaxt="n",xlab="Phenotype",ylab="-log10(p)")
abline( h=seq(0,YLIM[2],5), lty=2,col="grey50" )
legend( "topright", pch=PCHS, legend=c(rownames(SHAP.p[[1]]),"MNa (Resids)"), pt.cex=1.5,pt.lwd=3 ,ncol=2)
for ( c in 1:ncol(SHAP.p[[1]]) ) {
	points( jitter(rep(c,nrow(SHAP.p[[1]])),amount=.05), -log10(SHAP.p[[1]][,c]), col=COLS[c],pch=PCHS,cex=2,lwd=3 )
}
axis( 1, at=1:ncol(SHAP.p[[1]]), label=names(DEL), las=2 )
 # Transformed
for ( c in 1:ncol(SHAP.p[[2]]) ) {
	points( rep(c,nrow(SHAP.p[[2]])), -log10(SHAP.p[[2]][,c]), col=COLS[c],pch=PCHS[6],cex=2,lwd=3 )
}
dev.off()

###########################################################
## HOMOSCEDASTICITY PLOT ##################################
###########################################################

# P-Values for BP Test for Original and Transformed Data
COLS <- c("firebrick1","gold1","chartreuse1","dodgerblue1")
PCHS <- 0:7
png( paste(PathToSave,"/4b_Der_HomSK_Results.png",sep=""), height=1200, width=1200, pointsize=36 )
 # Original
XLIM <- c(1,ncol(BP.test) ) + c(-.4,.4)
YLIM <- c(0, -log10(Reduce( min, BP.test )) )
plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, main="Breusch-Pagan Test: Transformed",xaxt="n",xlab="Phenotype",ylab="-log10(p)")
abline( h=seq(0,YLIM[2],2), lty=2,col="grey50" )
for ( c in 1:ncol(BP.test) ) {
	points( rep(c,nrow(BP.test)), -log10(BP.test[,c]), col=COLS[c],pch=PCHS,cex=2,lwd=3 )
}
axis( 1, at=1:ncol(BP.test), label=names(DEL), las=2 )
legend( "topleft", pch=PCHS, legend=rownames(BP.test), pt.cex=1.5,pt.lwd=3 )
dev.off()




###########################################################
## END OF DOC #############################################
###########################################################
