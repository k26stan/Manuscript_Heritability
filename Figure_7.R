## Make Figure 7 (Permutation Results) for Resp_Herit Manuscript ##
## July 9, 2015 ##
## Kristopher Standish ##

## Gameplan

library( gplots )
library( xtable )
###########################################################
## LOAD DATA ##############################################
###########################################################

## Set Date
DATE <- gsub("-","",Sys.Date())

## Set Paths to Data and to Save
# PathToSing <- "/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Data/GCTA_Results/Single/4-PERM_Compile.Rdata"
# PathToDer <- "/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Data/GCTA_Results/Derived/4-PERM_Compile.Rdata"
# PathToResTabs <- "/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Plots/20160216_GCTA/TAB"
 # PC4 Model
PathToSing <- "/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Data/GCTA_Results/PC4_Sing/4-PERM_Compile.Rdata"
PathToDer <- "/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Data/GCTA_Results/PC4_Der/4-PERM_Compile.Rdata"
PathToResTabs <- "/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Plots/20160622_GCTA_PC4/TAB"
PathToSave <- paste("/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Plots/",DATE,"_GCTAperm_PC4",sep="")
 # No PC Model
PathToSing <- "/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Data/GCTA_Results/PC0_Sing/4-PERM_Compile.Rdata"
PathToDer <- "/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Data/GCTA_Results/PC0_Der/4-PERM_Compile.Rdata"
PathToResTabs <- "/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Plots/20160623_GCTA_PC0/TAB"
PathToSave <- paste("/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Plots/",DATE,"_GCTAperm_PC0",sep="")

dir.create( PathToSave )

## Load Data
load( PathToSing )
SING <- COMPILE.full
load( PathToDer )
DER <- COMPILE.full

## Load Table of GCTA Results
RES.sing <- read.table( paste(PathToResTabs,".sing.txt",sep=""), sep="\t",header=T )
RES.der <- read.table( paste(PathToResTabs,".der.txt",sep=""), sep="\t",header=T )

## Cohort Names
Cohort_Name <- "(MAF>1% - SNP+IND)" # "Full Cohort (MAF>1, SNP+Indel)"

###################################################
## FCT: MAKE PERMUTATION PLOTS ####################
###################################################

PLOT_PERM <- function(COMPILE, tag) {

	## Pull out Information from DATA
	VAR <- COMPILE$VAR
	SE <- COMPILE$SE
	MOD <- COMPILE$MOD
	PHENOS <- names(VAR)
	PHENOS.2 <- gsub("DEL_","",PHENOS)
	Num_Perms <- nrow(MOD[[1]])-1

	## Plotting Parameters for different Phenotypes
	 # colors
	COLS.list <- c("firebrick1","gold1","chartreuse1","dodgerblue1") # c("firebrick2","gold2","chartreuse2","deepskyblue2","slateblue2")
	COLS <- rep("firebrick1",length(VAR)) # COLS.list # rep( COLS.list[1], nrow(VAR) )
	COLS[ grep("SJC",names(VAR)) ] <- "gold4"
	COLS[ grep("TJC",names(VAR)) ] <- "dodgerblue4"
	COLS[ grep("CRP",names(VAR)) ] <- "chartreuse4"
	COLS[ grep("rSJC",names(VAR)) ] <- "gold1"
	COLS[ grep("rTJC",names(VAR)) ] <- "dodgerblue1"
	COLS[ grep("lCRP",names(VAR)) ] <- "chartreuse1"
	 # pch
	PCHS <- rep(1,length(VAR))
	if ( grepl("sing",tag, ignore.case=T) ) {
		PCHS <- numeric(length(VAR)) # rep( "o", length(VAR) )
		PCHS[ grep("WAG4",names(VAR)) ] <- 0
		PCHS[ grep("WAG12",names(VAR)) ] <- 1
		PCHS[ grep("WAG20",names(VAR)) ] <- 2
		PCHS[ grep("WAG28",names(VAR)) ] <- 3
		PCHS[ grep("FL",names(VAR)) ] <- 4
	}else{
		PCHS <- numeric(length(VAR)) 
		PCHS[ grep("MNa",names(VAR)) ] <- 0
		PCHS[ grep("MNcd",names(VAR)) ] <- 1
		PCHS[ grep("PRC",names(VAR)) ] <- 2
		PCHS[ grep("Bwk",names(VAR)) ] <- 3
		PCHS[ grep("VARwk",names(VAR)) ] <- 4

	}
	if ( grepl("alt",tag, ignore.case=T) ) {
		PCHS <- numeric(length(VAR)) 
		PCHS[ grep("MNa",names(VAR)) ] <- 0
		PCHS[ grep("MNcd",names(VAR)) ] <- 1
		PCHS[ grep("PRC",names(VAR)) ] <- 2
		PCHS[ grep("Bwk",names(VAR)) ] <- 3
		PCHS[ grep("VARwk",names(VAR)) ] <- 4
	}

	## Calculate Permuted P-Values
	P.perm.comp <- unlist(lapply( MOD, function(x) (1+length(which( abs(x[-nrow(x),"LRT"]) > abs(x[nrow(x),"LRT"])) )) / (nrow(x)) ))
	P.dat.comp <- unlist(lapply( MOD, function(x) x["True","Pval"] ))

	## Plot Permuted vs Actual P-Values
	P.dat.comp <- unlist(lapply( MOD, function(x) x["True","Pval"] ))
	PCHS.leg <- data.frame( LABS=unlist(lapply(strsplit(names(P.dat.comp),"_"),function(x) paste(x[1],x[2],sep="_"))), PCHS )
	PCHS.leg <- PCHS.leg[ which(!duplicated(PCHS.leg$PCHS)), ]
	COLS.leg <- data.frame( LABS=unlist(lapply(strsplit(names(P.dat.comp),"_"),function(x) x[3])), COLS )
	COLS.leg <- COLS.leg[ which(!duplicated(COLS.leg$COLS)), ]
	LIM <- c(0, max(4,max(-log10( c(P.dat.comp,P.perm.comp) ))) )
	 # Create Plot
	plot( 0,0, type="n", xlim=LIM,ylim=LIM, xlab="GCTA: -log10(p)",ylab="Permuted: -log10(p)", main="Permuted vs GCTA P-Values" )
	abline( h=seq(0,LIM[2]+1,1),lty=2,col="grey50") ; abline( v=seq(0,LIM[2]+1,1),lty=2,col="grey50")
	abline( 0,1, lty=1,col="black",lwd=2 )
	abline( h=-log10(.05), lty=2,col="magenta2",lwd=3 )
	abline( v=-log10(.05), lty=2,col="magenta2",lwd=3 )
	abline( h=-log10(1/(1+Num_Perms)), lty=2,col="chocolate2",lwd=3 )
	 # Populate w/ Data
	points( -log10(P.dat.comp), -log10(P.perm.comp), col=COLS,pch=PCHS,cex=2, lwd=3 )
	 # Create Legend
	legend( "topleft", legend=PCHS.leg$LABS[order(PCHS.leg$PCHS)], pch=PCHS.leg$PCHS[order(PCHS.leg$PCHS)], col="black", cex=1.2, pt.lwd=3,pt.cex=2 )
	legend( "bottomright", legend=COLS.leg$LABS[order(COLS.leg$COLS)], col=as.character(COLS.leg$COLS[order(COLS.leg$COLS)]), pch=20, cex=1.2, ncol=2 )

	## Return Data
	COMP_OUT <- list( DAT=P.dat.comp, PERM=P.perm.comp )
	return(COMP_OUT)
}


##########################################
## MAKE PERMUTATION PLOTS ################

## Save Figures Together in 1 Plot
png( paste(PathToSave,"/Perm_Pvals.Fig7.png",sep=""), width=2000,height=1000, pointsize=30 )
par(mfrow=c(1,2))

## Plot Permutations for Single Delta Stats
PHENOS_ALL <- names(SING$MOD)
WHICH_PHENOS <- Reduce( intersect, list(grep("DEL",PHENOS_ALL), grep("VARdr",PHENOS_ALL,invert=T), grep("JC28",PHENOS_ALL,invert=T)) )
SING.2 <- lapply( SING, function(x) x[WHICH_PHENOS] )
OUT.sing <- PLOT_PERM( SING.2, "sing" )

## Plot Permutations for Derived Delta Stats
PHENOS_ALL <- names(DER$MOD)
WHICH_PHENOS <- Reduce( intersect, list(grep("DEL",PHENOS_ALL), grep("VARdr",PHENOS_ALL,invert=T), grep("JC28",PHENOS_ALL,invert=T)) )
DER.2 <- lapply( DER, function(x) x[WHICH_PHENOS] )
OUT.der <- PLOT_PERM( DER.2, "der" )

dev.off()


###################################################
## COMPILE RESULTS TABLES #########################
###################################################

## Create Output Tables for Manuscript
 # Single Time Point
TAB.sing <- data.frame( RES.sing, P.perm=OUT.sing$PERM )
 # Derived Phenotypes
TAB.mn <- data.frame( RES.der, P.perm=OUT.der$PERM[grep("MNa",names(OUT.der$PERM))] )
 # Single + Mean
TAB.manu <- rbind( TAB.sing, TAB.mn )
TAB.manu.meta <- t(sapply(strsplit( rownames(TAB.manu),"_" ),"[",2:3))
TAB.manu <- data.frame( TAB.manu.meta[,2:1], TAB.manu )
colnames(TAB.manu)[1:2] <- c("Phenotype","Time")
PHENOS.uniq <- paste("_",as.character(unique(TAB.manu$Phenotype)),sep="")
TAB.manu <- Reduce( rbind, lapply( PHENOS.uniq, function(x)rbind(TAB.sing[grep(x,rownames(TAB.sing)),],TAB.mn[grep(x,rownames(TAB.mn)),]) ) )
TAB.manu.meta <- t(sapply(strsplit( rownames(TAB.manu),"_" ),"[",2:3))
TAB.manu <- data.frame( TAB.manu.meta[,2:1], TAB.manu )
colnames(TAB.manu)[1:2] <- c("Phenotype","Time")

## Write Table
TAB.manu.2 <- TAB.manu
rownames(TAB.manu.2) <- NULL
TAB.manu.2[,"Phenotype"] <- as.character(TAB.manu.2[,"Phenotype"])
TAB.manu.2[which(duplicated(TAB.manu.2$Phenotype)),"Phenotype"] <- ""
TAB.manu.2 <- cbind( TAB.manu.2, c("","*")[factor(TAB.manu.2$P<.05 & TAB.manu.2$P.perm<.05)] )
colnames(TAB.manu.2)[ncol(TAB.manu.2)] <- ""
# for ( col in 3:4 ) { TAB.manu.2[,col] <- round(TAB.manu.2[,col],2) }
# for ( col in 5:6 ) { TAB.manu.2[,col] <- round(TAB.manu.2[,col],4) }
write.table(TAB.manu.2, paste(PathToSave,"/TAB_PermResults.csv",sep=""), sep=",",row.names=F,col.names=T,quote=F )
writeLines( print(xtable(TAB.manu.2,digits=c(0,1,1,2,2,4,4,1)), include.rownames=F ), paste(PathToSave,"/TAB_PermResults.LaTeX.txt",sep="") )


###################################################
## END OF DOC #####################################
###################################################







	# FRAME <- data.frame( GCTA=P.dat.comp, PERM=P.perm.comp, PCHS, COLS )
	# # Plot Confidence Intervals of Permuted and Actual Data
	# XLIM <- c( 0,length(COMPILE$MOD)+1 )
	# YLIM <- c( -.5,1.5 )
	# # COLS <- c("firebrick1","gold1","chartreuse1","dodgerblue1","chartreuse1","dodgerblue1") 
	# plot( 0,0,type="n", ylim=YLIM,xlim=XLIM, main="Heritability Estimates for True/Permuted Data", ylab="Heritability Estimate (%)", xlab="Phenotype",xaxt="n" )
	# axis( 1, at=1:length(COMPILE$MOD), label=names(COMPILE$MOD),las=2 )
	# abline( h=seq(-2,2,.2),lty=2,col="grey50" )
	# abline( h=c(0,1),lty=1,col="black" )
	# for ( p in 1:length(COMPILE$MOD) ) {
	# 	arrows( p, COMPILE$VAR[[p]][,"VgVp"]+COMPILE$SE[[p]][,"VgVp"], p, COMPILE$VAR[[p]][,"VgVp"]-COMPILE$SE[[p]][,"VgVp"], code=3,angle=90, col=COLS[p],lwd=1 )
	# }
	# for ( p in 1:length(COMPILE$MOD) ) {
	# 	arrows( p, COMPILE$VAR[[p]][Num_Perms+1,"VgVp"]+COMPILE$SE[[p]][Num_Perms+1,"VgVp"], p, COMPILE$VAR[[p]][Num_Perms+1,"VgVp"]-COMPILE$SE[[p]][Num_Perms+1,"VgVp"], code=3,angle=90, col=gsub("1","4",COLS[p]),lwd=3 )
	# 	points( p, COMPILE$VAR[[p]][Num_Perms+1,"VgVp"], col=gsub("1","4",COLS[p]),pch=20 )
	# }




