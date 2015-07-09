## Make Figure 6 (GCTA Results) for Resp_Herit Manuscript ##
## July 8, 2015 ##
## Kristopher Standish ##

## Gameplan
 # Show discordance b/n ranks of patients at different timepoints
 # Show correlations b/n phenotypes
 # Show variance around mean/fit (after treatment?)
 # SIMULATIONS???

library( gplots )

###########################################################
## LOAD DATA ##############################################
###########################################################

## Set Date
DATE <- gsub("-","",Sys.Date())

## Set Paths to Data and to Save
PathToSing <- "/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Data/GCTA_Results/Single/GCTA_Estimates.ALL.Rdata"
PathToDer <- "/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Data/GCTA_Results/Derived/GCTA_Estimates.ALL.Rdata"
PathToSave <- paste("/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Plots/",DATE,"_GCTA",sep="")
dir.create( PathToSave )

## Load Data
load( PathToSing )
SING <- COMPILE
load( PathToDer )
DER <- COMPILE

## Cohort Names
Cohort_Name <- "(MAF>1% - SNP+IND)" # "Full Cohort (MAF>1, SNP+Indel)"
PHENOS <- c("DAS","lCRP","rSJC","rTJC")

###################################################
## FCT: PLOT HERITABILITY ESTIMATES ###############
###################################################

PLOT_GCTA <- function(DATA, tag) {

	## Pull out Information from DATA
	VAR <- DATA$VAR
	SE <- DATA$SE
	MOD <- DATA$MOD
	PHENOS <- rownames(VAR)
	PHENOS.2 <- gsub("DEL_","",PHENOS)

	## Basic Plot for (all) Phenotypes
	COLS.list <- c("firebrick1","gold1","chartreuse1","dodgerblue1") # c("firebrick2","gold2","chartreuse2","deepskyblue2","slateblue2")
	COLS <- rep( COLS.list[1], nrow(VAR) )
	COLS[ grep("SJC",PHENOS) ] <- "gold4"
	COLS[ grep("TJC",PHENOS) ] <- "dodgerblue4"
	COLS[ grep("CRP",PHENOS) ] <- "chartreuse4"
	COLS[ grep("rSJC",PHENOS) ] <- "gold1"
	COLS[ grep("rTJC",PHENOS) ] <- "dodgerblue1"
	COLS[ grep("lCRP",PHENOS) ] <- "chartreuse1"
	LTYS <- 1
	PCHS <- 20
	# LTYS <- rep( 2, nrow(VAR) )
	# LTYS[ grep("DEL",rownames(VAR)) ] <- 1
	# PCHS <- rep( 1, nrow(VAR) )
	# PCHS[ grep("DEL",rownames(VAR)) ] <- 20

	## Set Plot Parameters
	YLIM <- c( min( 0,VAR[,"VgVp"]-SE[,"VgVp"], na.rm=T), max(1,max(VAR[,"VgVp"]+SE[,"VgVp"],na.rm=T)) )
	YLIM <- c( -.3, 1.6 )
	XLIM <- c( 0,nrow(VAR)+1 )
	WHICH_SIG <- which( MOD[,"Pval"] < .05 )
	## Open Plot
	# png( paste(PathToSave,"/GCTA_Estimates.",tag,".png",sep=""), height=1200,width=800+40*length(PHENOS.2), pointsize=36 )
	plot( 0,0,type="n", xlim=XLIM, ylim=YLIM+c(0,.4), main=paste("Heritability Estimate -",Cohort_Name), ylab="% Phenotypic Variance", xlab="Phenotype", yaxt="n", xaxt="n" )
	 # Vertical Grid Lines
	abline( h=seq(-2,XLIM[2]+.4,.1), lty=2, col="grey50", lwd=1 )
	abline( h=c(0,1), lty=1, col="black", lwd=1 )
	 # Horizontal Lines/Data
	arrows( 1:nrow(VAR), 0, 1:nrow(VAR), 1, col="black", lwd=1, length=0 )
	arrows( 1:nrow(VAR), VAR[,"VgVp"]-SE[,"VgVp"], 1:nrow(VAR), VAR[,"VgVp"]+SE[,"VgVp"], col=COLS, lty=LTYS, code=3, angle=90, length=.15, lwd=6 )
	points( 1:nrow(VAR), VAR[,"VgVp"], col=COLS, pch=PCHS, cex=1.4, lwd=3 )
	 # Axis/Labels/Significance
	axis(2, at=seq(0,1,.2), las=2 )
	text( 1:nrow(VAR)-.35, .05+sapply( VAR[,"VgVp"]+SE[,"VgVp"], function(x) max(x,1) ), labels=PHENOS.2, pos=4, cex=.8, col=COLS, srt=90 )
	if ( length(WHICH_SIG)>0 ) { text( WHICH_SIG, rep(-.05,length(WHICH_SIG)), labels="*", col=COLS[WHICH_SIG], cex=1.5 ) }
	# dev.off()

}

###################################################
## PLOT HERITABILITY ESTIMATES ####################
###################################################

## Plot Heritability Estimates
png( paste(PathToSave,"/GCTA_Estimates.Fig6.png",sep=""), height=1200,width=2400, pointsize=36 )
layout( matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(11,4) )
 # Single Measure Data
PLOT_GCTA( SING, "Single" )
 # Derived Means
WHICH <- grep( "MNa", rownames(DER$VAR) )
 # Main
DER.mn <- lapply( DER[2:4], function(x) x[WHICH,] )
PLOT_GCTA( DER.mn, "Derived" )
dev.off()

## Plot Alternate Phenotypes
png( paste(PathToSave,"/GCTA_Estimates.Der.Alt.png",sep=""), height=1200,width=2400, pointsize=36 )
DER.alt <- lapply( DER[2:4], function(x) x[-WHICH,] )
PLOT_GCTA( DER.alt, "Derived.Alt" )
dev.off()









