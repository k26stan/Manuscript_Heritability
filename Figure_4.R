## Make Figure 4 (Single-Measurement Limitations) for Resp_Herit Manuscript ##
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
# PathToFT <- "/Users/kstandis/Data/Burn/Data/Phenos/Full_Tables/20141229_Full_Table.txt"
PathToFT <- "/Users/kstandis/Data/Janssen/Data/Pheno/Derived/20151015_Full_Table.txt"
# PathToRep <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150530_Resp_v_Time.txt"
PathToRep <- "/Users/kstandis/Data/Janssen/Data/Pheno/Derived/20160112_Resp_v_Time.txt"
PathToSing <- "/Users/kstandis/Data/Janssen/Data/Pheno/Derived/20150520_Single_Pheno_Table.txt"
PathToDer <- "/Users/kstandis/Data/Janssen/Data/Pheno/Derived/20150619_Derived_Pheno_Table.txt"
PathToSave <- paste("/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Plots/",DATE,"_SingleMeas/",sep="")
dir.create( PathToSave )

## Load Real Data
TAB <- read.table( PathToSing, sep="\t",header=T )
FT <- read.table( PathToFT, sep="\t",header=T )
RP <- read.table( PathToRep, sep="\t",header=T )
SING <- read.table( PathToSing, sep="\t",header=T )
DER <- read.table( PathToDer, sep="\t",header=T )

## Specify Phenotypes of Interest
WHICH_PHENOS <- c("DAS","lCRP","rSJC","rTJC")
COLS <- c("firebrick1","gold1","chartreuse1","dodgerblue1")
COLS.4 <- gsub("1","4",COLS)
COLS <- c("firebrick1","gold2","chartreuse1","dodgerblue1")

###########################################################
## PATIENT RANKS ##########################################
###########################################################

# ## Get Ranks for each Phenotype at each Timepoint
# RANKS <- list()
# for ( p in 1:length(WHICH_PHENOS) ) {
# 	pheno <- WHICH_PHENOS[p]
# 	RANKS[[pheno]] <- apply( SING[,intersect(grep("DEL",colnames(SING)),grep(pheno,colnames(SING)))], 2, rank )
# 	rownames(RANKS[[pheno]]) <- SING$IID
# } # Rank 1 is MOST improved

# ## Sort by 4WAG Measurement
# N <- nrow(RANKS[[1]])
# RANKS.sort.20 <- lapply( RANKS, function(x) x[order(x[,grep("WAG20",colnames(x))],decreasing=T),] )
# RANKS.sort.4 <- lapply( RANKS, function(x) x[order(x[,grep("WAG4",colnames(x))],decreasing=T),] )

# ## Combine to Single Data Frame
# RANKS.df.20 <- cbind( RANKS.sort.20[["DAS"]], RANKS.sort.20[["lCRP"]]+1*N, RANKS.sort.20[["rSJC"]]+2*N, RANKS.sort.20[["rTJC"]]+3*N )
# RANKS.df.4 <- cbind( RANKS.sort.4[["DAS"]], RANKS.sort.4[["lCRP"]]+1*N, RANKS.sort.4[["rSJC"]]+2*N, RANKS.sort.4[["rTJC"]]+3*N )

# ## Create Heatmap of Patient Improvement Ranks
#  # Colors
# COLS.DAS <- colorRampPalette(c("white",COLS[1],"black"))(N)
# COLS.lCRP <- colorRampPalette(c("white",COLS[2],"black"))(N)
# COLS.rSJC <- colorRampPalette(c("white",COLS[3],"black"))(N)
# COLS.rTJC <- colorRampPalette(c("white",COLS[4],"black"))(N)
# COLS.heat <- c( COLS.DAS, COLS.lCRP, COLS.rSJC, COLS.rTJC )
# COLS.brks <- seq( 0.5,max(RANKS.df)+.5,1 )
#  # Create Plot (20WAG)
# png( paste(PathToSave,"4_A_20WAG.png",sep="/"), height=1000,width=1600, pointsize=26 )
# heatmap.2( t(RANKS.df.20), col=COLS.heat, breaks=COLS.brks, main="Patient Rank in Improvement\nAfter Treatment",xlab="Pateint (By Rank 20WAG)", scale="none",trace="none",dendrogram="none",Rowv=F,Colv=F, labCol="",labRow=gsub("DEL_","",colnames(RANKS.df)), key=F,lhei=c(1,6),lwid=c(1,99),rowsep=seq(0,50,5),margins=c(2.5,10) )
# dev.off()
#  # Create Plot (4WAG)
# png( paste(PathToSave,"4_A_4WAG.png",sep="/"), height=1000,width=1600, pointsize=26 )
# heatmap.2( t(RANKS.df.4), col=COLS.heat, breaks=COLS.brks, main="Patient Rank in Improvement\nAfter Treatment",xlab="Pateint (By Rank 20WAG)", scale="none",trace="none",dendrogram="none",Rowv=F,Colv=F, labCol="",labRow=gsub("DEL_","",colnames(RANKS.df)), key=F,lhei=c(1,6),lwid=c(1,99),rowsep=seq(0,50,5),margins=c(2.5,10) )
# dev.off()

# # w <- 1
# # plot( RANKS.df.4[,1], RANKS.df.4[,1+w], col=COLS[1], pch="+" )
# # points( RANKS.df.4[,6]-N, RANKS.df.4[,6+w]-N, col=COLS[2], pch="+" )
# # points( RANKS.df.4[,11]-N*2, RANKS.df.4[,11+w]-N*2, col=COLS[3], pch="+" )
# # points( RANKS.df.4[,16]-N*3, RANKS.df.4[,16+w]-N*3, col=COLS[4], pch="+" )

# # abline( lm(RANKS.df.4[,1+w]~RANKS.df.4[,1]),col=COLS[1],lwd=2,lty=2 )
# # abline( lm(I(RANKS.df.4[,6+w]-N*1)~I(RANKS.df.4[,6]-N*1)),col=COLS[2],lwd=2,lty=2 )
# # abline( lm(I(RANKS.df.4[,11+w]-N*2)~I(RANKS.df.4[,11]-N*2)),col=COLS[3],lwd=2,lty=2 )
# # abline( lm(I(RANKS.df.4[,16+w]-N*3)~I(RANKS.df.4[,16]-N*3)),col=COLS[4],lwd=2,lty=2 )

###########################################################
## VARIANCE AROUND FIT ####################################
###########################################################

## Linear Regression w/ PLAC+DRUG as fixed effect Grouped by Patient
MODS <- list()
for ( p in 1:length(WHICH_PHENOS) ) {
	pheno <- WHICH_PHENOS[p]
	TEMP_DAT <- RP[, c("IID","DRUG","PLAC","WK",pheno) ] ; colnames(TEMP_DAT)[5] <- "Pheno"
	TEMP_DAT <- TEMP_DAT[ which(!is.na(TEMP_DAT[,"Pheno"])), ]
	# MOD <- lmList( Pheno ~ DRUG+PLAC | IID, data=TEMP_DAT )
	MOD <- lmList( Pheno ~ DRUG+PLAC+WK | IID, data=TEMP_DAT )
	# Compile
	MODS[[pheno]] <- MOD
}

## Calculate Variance Around Fit (including PLAC & DRUG)
SDS <- matrix( unlist(lapply( MODS, function(x) aggregate( resid(x), by=list(ID=names(resid(x))), sd, na.rm=T )[,2] )), ncol=length(MODS),byrow=F )

## Create Plot of Variance Distributions
 # Binsizes for each Phenotype
BINS <- c( .1, .05, .1, .1 )
# png( paste(PathToSave,"4_A.png",sep="/"), height=500,width=2000, pointsize=36 )
# par(mfrow=c(1,4))
png( paste(PathToSave,"4_A.png",sep="/"), height=1600,width=1600, pointsize=38 )
par(mfrow=c(2,2))
par(mar=c(4,4,3,1))
for ( i in 1:4 ) {
	XLIM <- range( SDS[,i], na.rm=T ) * c(1,1.05)
	BRKS <- seq( floor(XLIM[1]), XLIM[2]+BINS[i], BINS[i] )
	TEMP <- hist( SDS[,i], breaks=BRKS, plot=F )
	YLIM <- c(0,1.2*max(TEMP$counts))
	TEMP <- hist( SDS[,i], col=COLS[i], xlab="St.Dev of Residuals",ylab="# Patients",main=paste("Within Patient Variation:",WHICH_PHENOS[i]), breaks=BRKS,xlim=XLIM,ylim=YLIM )
	# abline( v=mean(SDS[,i],na.rm=T), lty=2,col=COLS.4[i],lwd=4 )
	# abline( v=mean(SDS[,i],na.rm=T), lty=3,col=COLS.4[i],lwd=8 )
	abline( v=median(SDS[,i],na.rm=T), lty=3,col="black",lwd=8 )
	text( median(SDS[,i],na.rm=T),1.15*max(TEMP$counts),paste("Range =",paste(round(range(SDS[,i],na.rm=T),2),collapse=" - ")), pos=4 )
	text( median(SDS[,i],na.rm=T),1.05*max(TEMP$counts),paste("Median =",round(median(SDS[,i],na.rm=T),2)), pos=4 )
}
dev.off()

###########################################################
## PHENOTYPE PAIR PLOTS ###################################
###########################################################

SING.corr.arr <- SING[,grep("DEL",colnames(SING))]
SING.corr.arr <- SING.corr.arr[ , which( sapply(strsplit(colnames(SING.corr.arr),"_"),"[",3) %in% WHICH_PHENOS ) ]
PHENOS <- c("DAS","lCRP","rSJC","rTJC")
for ( p in 1:length(PHENOS) ) { pheno <- PHENOS[p]
	DAT <- SING.corr.arr[,grep(pheno,colnames(SING.corr.arr))]
	png( paste(PathToSave,"/4_B.",pheno,".png",sep=""), height=2400,width=2400, pointsize=40 )
	par(mfrow=c(ncol(DAT),ncol(DAT)))
	par(mar=c(4,4,3,1))
	for ( c in 1:ncol(DAT) ) {
		col_c <- gsub("DEL_","",colnames(DAT)[c])
		for ( r in 1:ncol(DAT) ) {
			col_r <- gsub("DEL_","",colnames(DAT)[r])
			col_r <- gsub("_DAT","",col_r)
			if ( r==c ) { hist(DAT[,c], main=paste("Histogram:",col_c),xlab=col_c,col=COLS[p]) }
			if ( r<c ) { plot( DAT[,c], DAT[,r], pch="+",main=paste(col_r,"vs",col_c),xlab=col_c,ylab=col_r,col=COLS[p] ) ; abline(0,1,lty=2,lwd=2) }
			if ( r>c ) { plot( 0,0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",xaxt="n",yaxt="n" ) ; text(.5,.5, label=round(cor(DAT[,r],DAT[,c],use="pairwise.complete.obs",method="pearson"),2),col=COLS[p],cex=2.2 ) }
		}	
	}
	dev.off()
}

###########################################################
## PHENOTYPE CORRELATIONS #################################
###########################################################

## Calculate Correlations b/n Phenotypes
SING.corr.arr <- SING[,grep("DEL",colnames(SING))]
SING.corr.arr <- SING.corr.arr[ , which( sapply(strsplit(colnames(SING.corr.arr),"_"),"[",3) %in% WHICH_PHENOS ) ]
CORR <- cor( SING.corr.arr, use="pairwise.complete.obs",method="pearson" )
colnames(CORR) <- rownames(CORR) <- gsub("DEL_","",colnames(CORR))

## Create Heatmap of Patient Improvement Ranks
 # Colors
# COLS.list <- c("gold1","chocolate2","firebrick3","black","slateblue3","steelblue2","springgreen1") ; COLS.brks <- seq( -1,1,,length.out=201 )
COLS.list <- c("black","slateblue3","steelblue2","springgreen1","gold1","chocolate2","firebrick1") ; COLS.brks <- c( -1, seq( 0,1,length.out=200) )
# COLS.list <- c("black","chocolate1") ; COLS.brks <- c( -1, seq( 0,1,length.out=200) )
COLS.heat <- colorRampPalette(COLS.list)(200)
COLS.rows <- rep( COLS, rep(5,4) )
 # Create Plot
png( paste(PathToSave,"4_C.png",sep="/"), height=1600,width=1600, pointsize=36 )
heatmap.2( CORR, col=COLS.heat, breaks=COLS.brks, main="Correlation between\nResponse Phenotypes", 
	scale="none",trace="none",dendrogram="none",Rowv=F,Colv=F,lhei=c(1,5),lwid=c(1,5),margins=c(5,5),
	RowSideColors=COLS.rows,ColSideColors=COLS.rows,density.info="none",srtRow=-45,srtCol=45,
	cellnote=round(CORR,2),notecol="white",notecex=.8,
	colsep=seq(0,20,5),rowsep=seq(0,20,5) )
dev.off()



###########################################################
## ADDITIONAL FIGURES???? #################################
###########################################################

## Another Figure?

## Create Plot of ...
# png( paste(PathToSave,"4_D.png",sep="/"), height=700,width=2800, pointsize=36 )


# dev.off()



###########################################################
## MAKE FIGURE 4 ##########################################
###########################################################














###########################################################
## END OF DOC #############################################
###########################################################