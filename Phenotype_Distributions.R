## Plot Distributions of Phenotypes & Residuals ##
## Before and After Transformations ##
## March 4, 2015 ##
## Kristopher Standish ##

library( nlme )
library( gplots )
library( lmtest )

###########################################################
## LOAD DATA ##############################################
###########################################################

## Set Date
DATE <- gsub("-","",Sys.Date())

## Set Paths to Data and to Save
PathToFT <- "/Users/kstandis/Data/Burn/Data/Phenos/Full_Tables/20141229_Full_Table.txt"
PathToRep <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150226_Resp_v_Time.txt"
PathToSave <- paste("/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Plots/",DATE,sep="")
dir.create( PathToSave )

## Load Real Data
FT <- read.table( PathToFT, sep="\t",header=T )
RP <- read.table( PathToRep, sep="\t",header=T )

###########################################################
## ORGANIZE DATA ##########################################
###########################################################

## Get set of Weeks
WKS <- as.numeric( unique( RP$WK ) )

## Set number of patients (before Filtering)
Samps.1 <- as.character( FT$ID_2 )
N.samps.1 <- length( Samps.1 )

## Pull weeks participating in trial for each patient
IN <- FT$IN
BL <- c(0,24,16)[factor(FT$GRP)]
names(IN) <- names(BL) <- FT$ID_2

## Remove Patients w/ little or not Treatment Data
RM.drop <- which( IN < 8 )
RM.drop.names <- names(RM.drop)

## Set number of patients
Samps <- setdiff( Samps.1, RM.drop.names )
N.samps <- length( Samps )

## Calculate metrics at WAG & BL
INIT.0 <- INIT.BL <- array( , c(N.samps,7) )
rownames(INIT.0) <- rownames(INIT.BL) <- Samps
colnames(INIT.0) <- colnames(INIT.BL) <- c("DAS","CRP","SJC","TJC","lCRP","rSJC","rTJC")
WAG.4 <- WAG.12 <- WAG.20 <- WAG.28 <- WK.F <- array( , c(N.samps,7) )
rownames(WAG.4) <- rownames(WAG.12) <- rownames(WAG.20) <- rownames(WAG.28) <- rownames(WK.F) <- Samps
colnames(WAG.4) <- colnames(WAG.12) <- colnames(WAG.20) <- colnames(WAG.28) <- colnames(WK.F) <- c("DAS","CRP","SJC","TJC","lCRP","rSJC","rTJC")
for ( s in 1:N.samps ) {
	samp <- Samps[s]
	BL.samp <- BL[samp]
	TEMP_RP <- RP[ which(RP$IID==samp), ]
	 # Wk 0
	INIT.0[samp,c("DAS","CRP","SJC","TJC")] <- c( TEMP_RP[ which(TEMP_RP$WK==0), c("DAS","CRP","SJC","TJC") ], recursive=T )
	INIT.0[samp,"lCRP"] <- log10( TEMP_RP[ which(TEMP_RP$WK==0), "CRP" ] )
	INIT.0[samp,c("rSJC","rTJC")] <- c( sqrt( TEMP_RP[ which(TEMP_RP$WK==0), c("SJC","TJC") ] ), recursive=T )
	 # BL
	which_BL <- which(TEMP_RP$WK==BL.samp)
	if ( length(which_BL)>0 ) {
		INIT.BL[samp,c("DAS","CRP","SJC","TJC")] <- c( TEMP_RP[ which_BL, c("DAS","CRP","SJC","TJC") ], recursive=T )
		INIT.BL[samp,"lCRP"] <- log10( TEMP_RP[ which_BL, "CRP" ] )
		INIT.BL[samp,c("rSJC","rTJC")] <- c( sqrt( TEMP_RP[ which_BL, c("SJC","TJC") ] ), recursive=T )
	}
	# WAG
	which_WAG4 <- which(TEMP_RP$WK==BL.samp+4)
	if ( length(which_WAG4)>0 ) {
		WAG.4[samp,c("DAS","CRP","SJC","TJC")] <- c( TEMP_RP[ which_WAG4, c("DAS","CRP","SJC","TJC") ], recursive=T )
		WAG.4[samp,"lCRP"] <- log10( TEMP_RP[ which_WAG4, "CRP" ] )
		WAG.4[samp,c("rSJC","rTJC")] <- c( sqrt( TEMP_RP[ which_WAG4, c("SJC","TJC") ] ), recursive=T )
	}
	which_WAG12 <- which(TEMP_RP$WK==BL.samp+12)
	if ( length(which_WAG12)>0 ) {
		WAG.12[samp,c("DAS","CRP","SJC","TJC")] <- c( TEMP_RP[ which_WAG12, c("DAS","CRP","SJC","TJC") ], recursive=T )
		WAG.12[samp,"lCRP"] <- log10( TEMP_RP[ which_WAG12, "CRP" ] )
		WAG.12[samp,c("rSJC","rTJC")] <- c( sqrt( TEMP_RP[ which_WAG12, c("SJC","TJC") ] ), recursive=T )
	}
	which_WAG20 <- which(TEMP_RP$WK==BL.samp+20)
	if ( length(which_WAG20)>0 ) {
		WAG.20[samp,c("DAS","CRP","SJC","TJC")] <- c( TEMP_RP[ which_WAG20, c("DAS","CRP","SJC","TJC") ], recursive=T )
		WAG.20[samp,"lCRP"] <- log10( TEMP_RP[ which_WAG20, "CRP" ] )
		WAG.20[samp,c("rSJC","rTJC")] <- c( sqrt( TEMP_RP[ which_WAG20, c("SJC","TJC") ] ), recursive=T )
	}
	which_WAG28 <- which(TEMP_RP$WK==BL.samp+28)
	if ( length(which_WAG28)>0 ) {
		WAG.28[samp,c("DAS","CRP","SJC","TJC")] <- c( TEMP_RP[ which_WAG28, c("DAS","CRP","SJC","TJC") ], recursive=T )
		WAG.28[samp,"lCRP"] <- log10( TEMP_RP[ which_WAG28, "CRP" ] )
		WAG.28[samp,c("rSJC","rTJC")] <- c( sqrt( TEMP_RP[ which_WAG28, c("SJC","TJC") ] ), recursive=T )
	}
	# WK100
	which_WK_F <- which.max(TEMP_RP$WK)
	if ( length(which_WK_F)>0 ) {
		WK.F[samp,c("DAS","CRP","SJC","TJC")] <- c( TEMP_RP[ which_WK_F, c("DAS","CRP","SJC","TJC") ], recursive=T )
		WK.F[samp,"lCRP"] <- log10( TEMP_RP[ which_WK_F, "CRP" ] )
		WK.F[samp,c("rSJC","rTJC")] <- c( sqrt( TEMP_RP[ which_WK_F, c("SJC","TJC") ] ), recursive=T )
	}
}

###########################################################
## CALCULATE DIFFERENCES/RESIDS ###########################
###########################################################
 # 4, 12, 20, 28 WAG
 # First & Last
WAG_CATS <- c("WAG_4","WAG_12","WAG_20","WAG_28","FL")

###########################################################
## Distributions: Delta

## Create Tables of Delta-Values
DEL <- list()
DEL$DAS <- data.frame( WAG.4[,"DAS"]-INIT.BL[,"DAS"], WAG.12[,"DAS"]-INIT.BL[,"DAS"], WAG.20[,"DAS"]-INIT.BL[,"DAS"], WAG.28[,"DAS"]-INIT.BL[,"DAS"], WK.F[,"DAS"]-INIT.0[,"DAS"])
DEL$CRP <- data.frame( WAG.4[,"CRP"]-INIT.BL[,"CRP"], WAG.12[,"CRP"]-INIT.BL[,"CRP"], WAG.20[,"CRP"]-INIT.BL[,"CRP"], WAG.28[,"CRP"]-INIT.BL[,"CRP"], WK.F[,"CRP"]-INIT.0[,"CRP"])
DEL$SJC <- data.frame( WAG.4[,"SJC"]-INIT.BL[,"SJC"], WAG.12[,"SJC"]-INIT.BL[,"SJC"], WAG.20[,"SJC"]-INIT.BL[,"SJC"], WAG.28[,"SJC"]-INIT.BL[,"SJC"], WK.F[,"SJC"]-INIT.0[,"SJC"])
DEL$TJC <- data.frame( WAG.4[,"TJC"]-INIT.BL[,"TJC"], WAG.12[,"TJC"]-INIT.BL[,"TJC"], WAG.20[,"TJC"]-INIT.BL[,"TJC"], WAG.28[,"TJC"]-INIT.BL[,"TJC"], WK.F[,"TJC"]-INIT.0[,"TJC"])
for ( i in 1:length(DEL) ) { colnames(DEL[[i]]) <- WAG_CATS }

###########################################################
## Distributions: Residual vs Initial Values

## Calculate Residuals vs Initial Values
RES <- DEL
for ( p in 1:length(RES) ) {
	pheno <- names(RES)[p]
	for ( c in 1:4 ) {
		MOD <- lm( DEL[[p]][,c]	~ INIT.BL[,pheno] )
		which_rows <- as.numeric( names(resid(MOD)) )
		RES[[p]][which_rows,c] <- resid( MOD )
	}
	MOD <- lm( DEL[[p]][,5]	~ INIT.0[,pheno] )
	which_rows <- as.numeric( names(resid(MOD)) )
	RES[[p]][which_rows,5] <- resid( MOD )
}

###########################################################
## Distributions: Transformed Delta

## Create Tables of Transformed Delta-Values
DEL.t <- list()
DEL.t$DAS <- data.frame( WAG.4[,"DAS"]-INIT.BL[,"DAS"], WAG.12[,"DAS"]-INIT.BL[,"DAS"], WAG.20[,"DAS"]-INIT.BL[,"DAS"], WAG.28[,"DAS"]-INIT.BL[,"DAS"], WK.F[,"DAS"]-INIT.0[,"DAS"])
DEL.t$lCRP <- data.frame( WAG.4[,"lCRP"]-INIT.BL[,"lCRP"], WAG.12[,"lCRP"]-INIT.BL[,"lCRP"], WAG.20[,"lCRP"]-INIT.BL[,"lCRP"], WAG.28[,"lCRP"]-INIT.BL[,"lCRP"], WK.F[,"lCRP"]-INIT.0[,"lCRP"])
DEL.t$rSJC <- data.frame( WAG.4[,"rSJC"]-INIT.BL[,"rSJC"], WAG.12[,"rSJC"]-INIT.BL[,"rSJC"], WAG.20[,"rSJC"]-INIT.BL[,"rSJC"], WAG.28[,"rSJC"]-INIT.BL[,"rSJC"], WK.F[,"rSJC"]-INIT.0[,"rSJC"])
DEL.t$rTJC <- data.frame( WAG.4[,"rTJC"]-INIT.BL[,"rTJC"], WAG.12[,"rTJC"]-INIT.BL[,"rTJC"], WAG.20[,"rTJC"]-INIT.BL[,"rTJC"], WAG.28[,"rTJC"]-INIT.BL[,"rTJC"], WK.F[,"rTJC"]-INIT.0[,"rTJC"])
for ( i in 1:length(DEL.t) ) { colnames(DEL.t[[i]]) <- WAG_CATS }

###########################################################
## Distributions: Transformed Residual vs Initial Values

## Calculate Residuals vs Initial Values
 # And also test for Homoscedasticity (bptest)
RES.t <- DEL.t
BP.test <- array( ,c(5,4))
rownames(BP.test) <- colnames(RES.t[[1]])
colnames(BP.test) <- names(RES.t)
for ( p in 1:length(RES.t) ) {
	pheno <- names(RES.t)[p]
	for ( c in 1:4 ) {
		MOD <- lm( DEL.t[[p]][,c]	~ INIT.BL[,pheno] )
		which_rows <- as.numeric( names(resid(MOD)) )
		RES.t[[p]][which_rows,c] <- resid( MOD )
		BP <- bptest( MOD )$p.value
		BP.test[c,p] <- BP
	}
	# Calculate Residuals
	MOD <- lm( DEL.t[[p]][,5]	~ INIT.0[,pheno] )
	which_rows <- as.numeric( names(resid(MOD)) )
	RES.t[[p]][which_rows,5] <- resid( MOD )
	BP <- bptest( MOD )$p.value
	BP.test[5,p] <- BP
}

# ###########################################################
# ## Remove Patients w/ little or not Treatment Data
# RM.drop <- which( rownames(DEL[[1]]) %in% FT$ID_2[which( FT$IN<8 )] )
# for ( i in 1:length(DEL) ) {
# 	DEL[[i]] <- DEL[[i]][ -RM.drop, ]
# 	RES[[i]] <- RES[[i]][ -RM.drop, ]
# 	DEL.t[[i]] <- DEL.t[[i]][ -RM.drop, ]
# 	RES.t[[i]] <- RES.t[[i]][ -RM.drop, ]
# }

###########################################################
## FCT: PLOTS #############################################
###########################################################
COLS <- c("firebrick1","gold1","chartreuse1","dodgerblue1")

## Plot Function
SHIZ <- function(LIST, COLNAME, FILENAME, BIN_SIZES, LABEL) {
	P_VALS <- numeric(length(LIST))
	png( paste(PathToSave,"/Dist-",FILENAME,"_",COLNAME,".png",sep=""), height=700, width=2400, pointsize=30 )
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
SHAP.p[["1_Delta"]] <- array( , c(5,4) )
rownames(SHAP.p[["1_Delta"]]) <- WAG_CATS
## Plot it
SHAP.p[["1_Delta"]]["WAG_4",] <- SHIZ( DEL, "WAG_4","1_Delta", c(.25,2,1,1), "Delta" )
SHAP.p[["1_Delta"]]["WAG_12",] <- SHIZ( DEL, "WAG_12","1_Delta", c(.25,2,1,1), "Delta" )
SHAP.p[["1_Delta"]]["WAG_20",] <- SHIZ( DEL, "WAG_20","1_Delta", c(.25,2,1,1), "Delta" )
SHAP.p[["1_Delta"]]["WAG_28",] <- SHIZ( DEL, "WAG_28","1_Delta", c(.25,2,1,1), "Delta" )
SHAP.p[["1_Delta"]]["FL",] <- SHIZ( DEL, "FL","1_Delta", c(.25,2,1,1), "Delta" )

###########################################################
## Distributions: Residuals
SHAP.p[["2_Resid"]] <- array( , c(5,4) )
rownames(SHAP.p[["2_Resid"]]) <- WAG_CATS
## Plot it
SHAP.p[["2_Resid"]]["WAG_4",] <- SHIZ( RES, "WAG_4","2_Resid", c(.25,2,1,1), "Residuals" )
SHAP.p[["2_Resid"]]["WAG_12",] <- SHIZ( RES, "WAG_12","2_Resid", c(.25,2,1,1), "Residuals" )
SHAP.p[["2_Resid"]]["WAG_20",] <- SHIZ( RES, "WAG_20","2_Resid", c(.25,2,1,1), "Residuals" )
SHAP.p[["2_Resid"]]["WAG_28",] <- SHIZ( RES, "WAG_28","2_Resid", c(.25,2,1,1), "Residuals" )
SHAP.p[["2_Resid"]]["FL",] <- SHIZ( RES, "FL","2_Resid", c(.25,2,1,1), "Residuals" )

###########################################################
## Distributions: Transformed Delta
SHAP.p[["3_TrDelta"]] <- array( , c(5,4) )
rownames(SHAP.p[["3_TrDelta"]]) <- WAG_CATS
## Plot it
SHAP.p[["3_TrDelta"]]["WAG_4",] <- SHIZ( DEL.t, "WAG_4","3_TrDelta", c(.25,.25,.5,.5), "Delta" )
SHAP.p[["3_TrDelta"]]["WAG_12",] <- SHIZ( DEL.t, "WAG_12","3_TrDelta", c(.25,.25,.5,.5), "Delta" )
SHAP.p[["3_TrDelta"]]["WAG_20",] <- SHIZ( DEL.t, "WAG_20","3_TrDelta", c(.25,.25,.5,.5), "Delta" )
SHAP.p[["3_TrDelta"]]["WAG_28",] <- SHIZ( DEL.t, "WAG_28","3_TrDelta", c(.25,.25,.5,.5), "Delta" )
SHAP.p[["3_TrDelta"]]["FL",] <- SHIZ( DEL.t, "FL","3_TrDelta", c(.25,.25,.5,.5), "Delta" )

###########################################################
## Shapiro's Test for Normality
SHAP.p[["4_TrResid"]] <- array( , c(5,4) )
rownames(SHAP.p[["4_TrResid"]]) <- WAG_CATS
## Plot it
SHAP.p[["4_TrResid"]]["WAG_4",] <- SHIZ( RES.t, "WAG_4","4_TrResid", c(.25,.25,.5,.5), "Residuals" )
SHAP.p[["4_TrResid"]]["WAG_12",] <- SHIZ( RES.t, "WAG_12","4_TrResid", c(.25,.25,.5,.5), "Residuals" )
SHAP.p[["4_TrResid"]]["WAG_20",] <- SHIZ( RES.t, "WAG_20","4_TrResid", c(.25,.25,.5,.5), "Residuals" )
SHAP.p[["4_TrResid"]]["WAG_28",] <- SHIZ( RES.t, "WAG_28","4_TrResid", c(.25,.25,.5,.5), "Residuals" )
SHAP.p[["4_TrResid"]]["FL",] <- SHIZ( RES.t, "FL","4_TrResid", c(.25,.25,.5,.5), "Residuals" )

## Plot Results
# P-Values for Shapiro Test for Original and Transformed Data
COLS <- c("firebrick1","gold1","chartreuse1","dodgerblue1")
PCHS <- 0:4
png( paste(PathToSave,"/Shapiro_Results.png",sep=""), height=1200, width=2400, pointsize=30 )
par(mfrow=c(1,2))
 # Original
XLIM <- c(1,ncol(SHAP.p[[2]]) )
YLIM <- c(0, -log10(Reduce( min, SHAP.p[c(2,4)] )) )
plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, main="Shapiro Test of Residuals: Original",xaxt="n",xlab="Phenotype",ylab="-log10(p)")
abline( h=seq(0,YLIM[2],5), lty=2,col="grey50" )
for ( c in 1:ncol(SHAP.p[[2]]) ) {
	points( rep(c,nrow(SHAP.p[[2]])), -log10(SHAP.p[[2]][,c]), col=COLS[c],pch=PCHS,cex=2,lwd=3 )
}
axis( 1, at=1:ncol(SHAP.p[[2]]), label=names(DEL), las=2 )
 # Transformed
XLIM <- c(1,ncol(SHAP.p[[4]]) )
YLIM <- c(0, -log10(Reduce( min, SHAP.p[c(2,4)] )) )
plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, main="Shapiro Test of Residuals: Transformed",xaxt="n",xlab="Phenotype",ylab="-log10(p)")
abline( h=seq(0,YLIM[2],5), lty=2,col="grey50" )
for ( c in 1:ncol(SHAP.p[[4]]) ) {
	points( rep(c,nrow(SHAP.p[[4]])), -log10(SHAP.p[[4]][,c]), col=COLS[c],pch=PCHS,cex=2,lwd=3 )
}
axis( 1, at=1:ncol(SHAP.p[[4]]), label=names(DEL.t), las=2 )
legend( "topleft", pch=PCHS, legend=rownames(SHAP.p[[4]]), pt.cex=1.6,pt.lwd=3 )
dev.off()

# # QQ
# LIM <- c( 0, ceiling( -log10(Reduce(min,SHAP.p[c(2,4)])) ) )
# plot( 0,0,type="n", xlim=LIM, ylim=LIM, xlab="Expected",ylab="Observed" )
# abline( h=seq( 0,LIM[2],5 ), lty=2,col="grey50" )
# abline( v=seq( 0,LIM[2],5 ), lty=2,col="grey50" )
# abline( 0,1, lty=1,lwd=2,col="black" )
# points( -log10(1:20/20), -log10(sort(SHAP.p[[2]])), col="mediumpurple2",pch="+" )
# points( -log10(1:20/20), -log10(sort(SHAP.p[[4]])), col="cadetblue2",pch="+" )

###########################################################
## HOMOSCEDASTICITY PLOT ##################################
###########################################################
COLS <- c("firebrick1","gold1","chartreuse1","dodgerblue1")
COLS.4 <- gsub("1","3",COLS)
## Plots for Homoscedasticity
 # Similar to Distribution Plots
for ( c in 1:4 ) {
	# WAG
	c_name <- colnames(RES.t[[1]])[c]
	png( paste(PathToSave,"/HomSK-",c_name,".png",sep=""), height=700, width=2400, pointsize=30 )
	par(mfrow=c(1,4))
	for ( p in 1:4 ) {
		pheno <- names(RES.t)[p]
		MOD <- lm( DEL.t[[p]][,c] ~ INIT.BL[,pheno] )
		YLIM <- range( RES.t[[p]][,c], na.rm=T )
		XLIM <- range( INIT.BL[,pheno], na.rm=T )
		plot( RES.t[[p]][,c] ~ INIT.BL[,pheno], xlim=XLIM,ylim=YLIM, main=paste("Residuals vs Initial",pheno,"-",c_name), xlab=paste("Delta",pheno),ylab=paste("Initial",pheno), col=COLS[p], pch="+" )
		abline( h=seq(-10,10,1), lty=2,lwd=1,col="grey50" )
		abline( h=0, lty=2,lwd=2,col=COLS.4[p] )
		text( XLIM[1],quantile(YLIM,.9), pos=4, label=paste("Breusch-Pagan:",formatC(BP.test[c_name,pheno],format="e",digits=2)) )
	}
	dev.off()
}
# FL
c_name <- colnames(RES.t[[1]])[5]
png( paste(PathToSave,"/HomSK-",c_name,".png",sep=""), height=700, width=2400, pointsize=30 )
par(mfrow=c(1,4))
for ( p in 1:4 ) {
	# Calculate Residuals
	pheno <- names(RES.t)[p]
	MOD <- lm( DEL.t[[p]][,5] ~ INIT.BL[,pheno] )
	YLIM <- range( RES.t[[p]][,5], na.rm=T )
	XLIM <- range( INIT.0[,pheno], na.rm=T )
	plot( RES.t[[p]][,5] ~ INIT.0[,pheno], main=paste("Residuals vs Initial",pheno,"-",c_name), xlab=paste("Delta",pheno),ylab=paste("Initial",pheno), col=COLS[p], pch="+" )
	abline( h=seq(-10,10,1), lty=2,lwd=1,col="grey50" )
	abline( h=0, lty=2,lwd=2,col=COLS.4[p] )
	text( XLIM[1],quantile(YLIM,.9), pos=4, label=paste("Breusch-Pagan:",formatC(BP.test[c_name,pheno],format="e",digits=2)) )
}
dev.off()

# P-Values for BP Test for Original and Transformed Data
COLS <- c("firebrick1","gold1","chartreuse1","dodgerblue1")
PCHS <- 0:4
png( paste(PathToSave,"/HomSK_Results.png",sep=""), height=1200, width=1200, pointsize=30 )
 # Original
XLIM <- c(1,ncol(BP.test) )
YLIM <- c(0, -log10(Reduce( min, BP.test )) )
plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, main="Breusch-Pagan Test: Transformed",xaxt="n",xlab="Phenotype",ylab="-log10(p)")
abline( h=seq(0,YLIM[2],2), lty=2,col="grey50" )
for ( c in 1:ncol(BP.test) ) {
	points( rep(c,nrow(BP.test)), -log10(BP.test[,c]), col=COLS[c],pch=PCHS,cex=2,lwd=3 )
}
axis( 1, at=1:ncol(BP.test), label=names(DEL), las=2 )
legend( "topleft", pch=PCHS, legend=rownames(BP.test), pt.cex=1.6,pt.lwd=3 )
dev.off()

###########################################################
## PLOT CORRELATION b/n TIME POINTS #######################
###########################################################

## Calculate Correlations amongst Phenotypes
 # (after Transformation)
FRAME <- data.frame( DEL.t$DAS, DEL.t$lCRP, DEL.t$rSJC, DEL.t$rTJC )
FRAME.names <- paste( rep(names(DEL.t),rep(5,4)), rep(colnames(DEL.t$DAS),4), sep="_" )
colnames(FRAME) <- FRAME.names
CORR.t <- cor( FRAME, use="pairwise.complete.obs", method="spearman" )

## Heatmap Correlation
COLS.list <- c("black","slateblue3","steelblue2","springgreen2","gold2","chocolate2","firebrick1")
COLS <- colorRampPalette(COLS.list)(100)
BRKS <- seq( 0,1,length.out=101 )
png( paste(PathToSave,"/Corr_tPhenos.png",sep=""), height=1600, width=1600, pointsize=30 )
heatmap.2( CORR.t, col=COLS, trace="none",scale="none", Colv=F,Rowv=F,dendrogram="none", margins=c(8,8), main="Correlation b/n Single Measurements", lhei=c(1,5),lwid=c(1,5) )
dev.off()

###########################################################
## WRITE TABLES OF PHENOTYPES #############################
###########################################################

## Compile Phenotype Measurements into Single Table
INIT_TAB <- data.frame( INIT.BL, INIT.0 )
colnames(INIT_TAB) <- paste( rep(c("Ibl","I0"),rep(ncol(INIT.0),2)), colnames(INIT.BL), sep="_" )
WAG_TAB <- data.frame( WAG.4, WAG.12, WAG.20, WAG.28 )
colnames(WAG_TAB) <- paste( rep(c("WAG4","WAG12","WAG20","WAG28"),rep(ncol(INIT.0),4)), colnames(WAG.4), sep="_" )
DEL_TAB <- data.frame( DEL.t$DAS, DEL.t$lCRP, DEL.t$rSJC, DEL.t$rTJC )
colnames(DEL_TAB) <- paste( "DEL",rep(c("WAG4","WAG12","WAG20","WAG28","FL"),4),rep(names(DEL.t),rep(5,4)), sep="_" )
DEL_TAB.2 <- data.frame( DEL$CRP, DEL$SJC, DEL$TJC )
colnames(DEL_TAB.2) <- paste( "DEL",rep(c("WAG4","WAG12","WAG20","WAG28","FL"),3),rep(names(DEL)[2:4],rep(5,3)), sep="_" )

FULL_TAB <- data.frame( INIT_TAB, WAG_TAB, DEL_TAB, DEL_TAB.2 )
for ( c in 1:ncol(FULL_TAB) ) { print(paste(colnames(FULL_TAB)[c],"-",length(which(is.na( FULL_TAB[,c] ))))) }

## Write Table
FULL_TAB.w <- data.frame( IID=rownames(FULL_TAB), FID=rownames(FULL_TAB), FULL_TAB )
write.table( FULL_TAB.w, gsub("20141229_Full_Table.txt",paste(DATE,"_Single_Pheno_Table.txt",sep=""),PathToFT), sep="\t",row.names=F,col.names=T,quote=F )
 # Write Phenotype List
FULL_TAB.colnames <- colnames(FULL_TAB)
write.table( data.frame(FULL_TAB.colnames), gsub("20141229_Full_Table.txt",paste(DATE,"Pheno_List.txt",sep=""),PathToFT), sep="\t",row.names=F,col.names=F,quote=F )



###########################################################
## END OF DOC #############################################
###########################################################
