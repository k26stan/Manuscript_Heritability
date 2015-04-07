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
DATE <- "20150313"

## Set Paths to Data and to Save
PathToFT <- "/Users/kstandis/Data/Burn/Data/Phenos/Full_Tables/20141229_Full_Table.txt"
PathToRep <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150226_Resp_v_Time.txt"
PathToWAG <- "/Users/kstandis/Data/Burn/Data/Phenos/Full_Tables/20150310_Single_Pheno_Table.txt"
PathToDER <- "/Users/kstandis/Data/Burn/Data/Phenos/Full_Tables/20150313_Derived_Pheno_Table.txt"
PathToSave <- paste("/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Plots/",DATE,sep="")

## Load Real Data
FT <- read.table( PathToFT, sep="\t",header=T )
RP <- read.table( PathToRep, sep="\t",header=T )
WAG <- read.table( PathToWAG, sep="\t",header=T )
DER <- read.table( PathToDER, sep="\t",header=T )

MG <- merge( FT, DER, by.x="ID_2",by.y="FID" )
###########################################################
## ORGANIZE DATA ##########################################
###########################################################

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
DEL_CATS <- c("MNw","MNwo","MNa","MNcd","Bdr","PRC","Bwk","VARdr","VARwk")
DEL <- list()
DEL$DAS <- data.frame( MG[, colnames(DER)[grep("DEL.*DAS",colnames(DER))] ] )
DEL$lCRP <- data.frame( MG[, colnames(DER)[grep("DEL.*CRP",colnames(DER))] ] )
DEL$rSJC <- data.frame( MG[, colnames(DER)[setdiff( grep("DEL.*SJC",colnames(DER)),grep("28",colnames(DER)) ) ] ] )
DEL$rTJC <- data.frame( MG[, colnames(DER)[setdiff( grep("DEL.*TJC",colnames(DER)),grep("28",colnames(DER)) ) ] ] )
for ( i in 1:length(DEL) ) { colnames(DEL[[i]]) <- DEL_CATS }

## Create Tables of Post-Treatment Values
POST_CATS <- c("MNw","MNwo","MNa","Bdr")
POST <- list()
POST$DAS <- data.frame( MG[, colnames(DER)[grep("POST.*DAS",colnames(DER))] ] )
POST$lCRP <- data.frame( MG[, colnames(DER)[grep("POST.*CRP",colnames(DER))] ] )
POST$rSJC <- data.frame( MG[, colnames(DER)[setdiff( grep("POST.*SJC",colnames(DER)),grep("28",colnames(DER)) ) ] ] )
POST$rTJC <- data.frame( MG[, colnames(DER)[setdiff( grep("POST.*TJC",colnames(DER)),grep("28",colnames(DER)) ) ] ] )
for ( i in 1:length(POST) ) { colnames(POST[[i]]) <- POST_CATS }

###########################################################
## Distributions: Residual vs Initial Values

## Create table to write Pheno/Covs to (for GCTA)
WRITE_PHENO_TAB <- array( ,c(0,2) )

## Categories of Pheno using Pre- as Covariate
RES_CATS <- c("MNw","MNwo","MNa","Bdr")

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
	RES[[p]] <- RES[[p]][, RES_CATS] # grep( paste(RES_CATS,collapse="|"), colnames(RES$DAS) ) ]
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

## Plot Results
# P-Values for Shapiro Test for Original and Transformed Data
COLS <- c("firebrick1","gold1","chartreuse1","dodgerblue1")
PCHS <- 0:7
png( paste(PathToSave,"/4b_Der_Shapiro_Results.png",sep=""), height=1200, width=2400, pointsize=30 )
par(mfrow=c(1,2))
 # Original
XLIM <- c(1,ncol(SHAP.p[[1]]) )
YLIM <- c(0, -log10(Reduce( min, SHAP.p[[1]] )) )
plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, main="Shapiro Test: Derived Stats",xaxt="n",xlab="Phenotype",ylab="-log10(p)")
abline( h=seq(0,YLIM[2],5), lty=2,col="grey50" )
for ( c in 1:ncol(SHAP.p[[1]]) ) {
	points( rep(c,nrow(SHAP.p[[1]])), -log10(SHAP.p[[1]][,c]), col=COLS[c],pch=PCHS,cex=2,lwd=3 )
}
axis( 1, at=1:ncol(SHAP.p[[1]]), label=names(DEL), las=2 )
legend( "topleft", pch=PCHS, legend=rownames(SHAP.p[[1]]), pt.cex=1.5,pt.lwd=3 )
 # Transformed
XLIM <- c(1,ncol(SHAP.p[[2]]) )
YLIM <- c(0, -log10(Reduce( min, SHAP.p[[2]] )) )
plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, main="Shapiro Test: Residuals of Derived Stats",xaxt="n",xlab="Phenotype",ylab="-log10(p)")
abline( h=seq(0,YLIM[2],5), lty=2,col="grey50" )
for ( c in 1:ncol(SHAP.p[[2]]) ) {
	points( rep(c,nrow(SHAP.p[[2]])), -log10(SHAP.p[[2]][,c]), col=COLS[c],pch=PCHS,cex=2,lwd=3 )
}
axis( 1, at=1:ncol(SHAP.p[[2]]), label=names(RES), las=2 )
legend( "topleft", pch=PCHS, legend=rownames(SHAP.p[[2]]), pt.cex=1.5,pt.lwd=3 )
dev.off()

# # QQ
# LIM <- c( 0, ceiling( -log10(Reduce(min,SHAP.p)) ) )
# plot( 0,0,type="n", xlim=LIM, ylim=LIM, xlab="Expected",ylab="Observed" )
# abline( h=seq( 0,LIM[2],5 ), lty=2,col="grey50" )
# abline( v=seq( 0,LIM[2],5 ), lty=2,col="grey50" )
# abline( 0,1, lty=1,lwd=2,col="black" )
# points( -log10(1:prod(dim(SHAP.p[[1]]))/prod(dim(SHAP.p[[1]]))), -log10(sort(SHAP.p[[1]])), col="mediumpurple2",pch="+" )
# points( -log10(1:prod(dim(SHAP.p[[2]]))/prod(dim(SHAP.p[[2]]))), -log10(sort(SHAP.p[[2]])), col="cadetblue2",pch="+" )

###########################################################
## HOMOSCEDASTICITY PLOT ##################################
###########################################################
# COLS <- c("firebrick1","gold1","chartreuse1","dodgerblue1")
# COLS.4 <- gsub("1","3",COLS)
# ## Plots for Homoscedasticity
#  # Similar to Distribution Plots
# for ( c in 1:4 ) {
# 	# WAG
# 	c_name <- colnames(RES[[1]])[c]
# 	png( paste(PathToSave,"/4b_Der_HomSK-",c_name,".png",sep=""), height=700, width=2400, pointsize=30 )
# 	par(mfrow=c(1,4))
# 	for ( p in 1:4 ) {
# 		pheno <- names(RES)[p]
# 		MOD <- lm( DEL[[p]][,c] ~ INIT.BL[,pheno] )
# 		YLIM <- range( RES[[p]][,c], na.rm=T )
# 		XLIM <- range( INIT.BL[,pheno], na.rm=T )
# 		plot( RES[[p]][,c] ~ INIT.BL[,pheno], xlim=XLIM,ylim=YLIM, main=paste("Residuals vs Initial",pheno,"-",c_name), xlab=paste("Delta",pheno),ylab=paste("Initial",pheno), col=COLS[p], pch="+" )
# 		abline( h=seq(-10,10,1), lty=2,lwd=1,col="grey50" )
# 		abline( h=0, lty=2,lwd=2,col=COLS.4[p] )
# 		text( XLIM[1],quantile(YLIM,.9), pos=4, label=paste("Breusch-Pagan:",formatC(BPest[c_name,pheno],format="e",digits=2)) )
# 	}
# 	dev.off()
# }
# # FL
# c_name <- colnames(RES[[1]])[5]
# png( paste(PathToSave,"/4b_Der_HomSK-",c_name,".png",sep=""), height=700, width=2400, pointsize=30 )
# par(mfrow=c(1,4))
# for ( p in 1:4 ) {
# 	# Calculate Residuals
# 	pheno <- names(RES)[p]
# 	MOD <- lm( DEL[[p]][,5] ~ INIT.BL[,pheno] )
# 	YLIM <- range( RES[[p]][,5], na.rm=T )
# 	XLIM <- range( INIT.0[,pheno], na.rm=T )
# 	plot( RES[[p]][,5] ~ INIT.0[,pheno], main=paste("Residuals vs Initial",pheno,"-",c_name), xlab=paste("Delta",pheno),ylab=paste("Initial",pheno), col=COLS[p], pch="+" )
# 	abline( h=seq(-10,10,1), lty=2,lwd=1,col="grey50" )
# 	abline( h=0, lty=2,lwd=2,col=COLS.4[p] )
# 	text( XLIM[1],quantile(YLIM,.9), pos=4, label=paste("Breusch-Pagan:",formatC(BPest[c_name,pheno],format="e",digits=2)) )
# }
# dev.off()

# P-Values for BP Test for Original and Transformed Data
COLS <- c("firebrick1","gold1","chartreuse1","dodgerblue1")
PCHS <- 0:7
png( paste(PathToSave,"/4b_Der_HomSK_Results.png",sep=""), height=1200, width=1200, pointsize=30 )
 # Original
XLIM <- c(1,ncol(BP.test) )
YLIM <- c(0, -log10(Reduce( min, BP.test )) )
plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, main="Breusch-Pagan Test: Transformed",xaxt="n",xlab="Phenotype",ylab="-log10(p)")
abline( h=seq(0,YLIM[2],2), lty=2,col="grey50" )
for ( c in 1:ncol(BP.test) ) {
	points( rep(c,nrow(BP.test)), -log10(BP.test[,c]), col=COLS[c],pch=PCHS,cex=2,lwd=3 )
}
axis( 1, at=1:ncol(BP.test), label=names(DEL), las=2 )
legend( "topleft", pch=PCHS, legend=rownames(BP.test), pt.cex=1.5,pt.lwd=3 )
dev.off()

# ###########################################################
# ## PLOT CORRELATION b/n TIME POINTS #######################
# ###########################################################

# ## Calculate Correlations amongst Phenotypes
#  # (after Transformation)
# FRAME <- data.frame( DEL.t$DAS, DEL.t$lCRP, DEL.t$rSJC, DEL.t$rTJC )
# FRAME.names <- paste( rep(names(DEL.t),rep(5,4)), rep(colnames(DEL.t$DAS),4), sep="_" )
# colnames(FRAME) <- FRAME.names
# CORR.t <- cor( FRAME, use="pairwise.complete.obs", method="spearman" )

# ## Heatmap Correlation
# COLS.list <- c("black","slateblue3","steelblue2","springgreen2","gold2","chocolate2","firebrick1")
# COLS <- colorRampPalette(COLS.list)(100)
# BRKS <- seq( 0,1,length.out=101 )
# png( paste(PathToSave,"/Corr_tPhenos.png",sep=""), height=1600, width=1600, pointsize=30 )
# heatmap.2( CORR.t, col=COLS, trace="none",scale="none", Colv=F,Rowv=F,dendrogram="none", margins=c(8,8), main="Correlation b/n Single Measurements", lhei=c(1,5),lwid=c(1,5) )
# dev.off()

# ###########################################################
# ## WRITE TABLES OF PHENOTYPES #############################
# ###########################################################

# ## Compile Phenotype Measurements into Single Table
# INIT_TAB <- data.frame( INIT.BL, INIT.0 )
# colnames(INIT_TAB) <- paste( rep(c("Ibl","I0"),rep(ncol(INIT.0),2)), colnames(INIT.BL), sep="_" )
# WAG_TAB <- data.frame( WAG.4, WAG.12, WAG.20, WAG.28 )
# colnames(WAG_TAB) <- paste( rep(c("WAG4","WAG12","WAG20","WAG28"),rep(ncol(INIT.0),4)), colnames(WAG.4), sep="_" )
# DEL_TAB <- data.frame( DEL.t$DAS, DEL.t$lCRP, DEL.t$rSJC, DEL.t$rTJC )
# colnames(DEL_TAB) <- paste( "DEL",rep(c("WAG4","WAG12","WAG20","WAG28","FL"),4),rep(names(DEL.t),rep(5,4)), sep="_" )

# FULL_TAB <- data.frame( INIT_TAB, WAG_TAB, DEL_TAB )
# for ( c in 1:ncol(FULL_TAB) ) { print(paste(colnames(FULL_TAB)[c],"-",length(which(is.na( FULL_TAB[,c] ))))) }

# ## Write Table
# FULL_TAB.w <- data.frame( IID=rownames(FULL_TAB), FID=rownames(FULL_TAB), FULL_TAB )
# write.table( FULL_TAB.w, gsub("20141229_Full_Table.txt","20150310_Single_Pheno_Table.txt",PathToFT), sep="\t",row.names=F,col.names=T,quote=F )
#  # Write Phenotype List
# FULL_TAB.colnames <- colnames(FULL_TAB)
# write.table( data.frame(FULL_TAB.colnames), gsub("20141229_Full_Table.txt","20150310_Pheno_List.txt",PathToFT), sep="\t",row.names=F,col.names=F,quote=F )



###########################################################
## END OF DOC #############################################
###########################################################

COLS <- c("firebrick1","chocolate1","gold1","chartreuse1","dodgerblue1")
COLS <- c("firebrick1","gold1","dodgerblue1")
COLS.heat <- colorRampPalette(COLS)(100)
## Normalize to Mean (by Measurement)
DER.2 <- DER[,3:ncol(DER)]
rownames(DER.2) <- DER$FID
for ( c in 1:ncol(DER.2) ) {
	TEMP <- ( DER.2[,c] - mean(DER.2[,c]) ) / sd(DER.2[,c])
	DER.2[,c] <- TEMP
}
CORR.2 <- cor( t(DER.2), use="pairwise.complete.obs",method="spearman" )
heatmap.2( CORR.2, scale="none",trace="none", col=COLS.heat, RowSideColors=colorRampPalette(c("black","red1"))(nrow(CORR.2)) )

COLUMNS <- sample( 1:ncol(DER.2), 5 )
pairs( DER.2[,COLUMNS] )












