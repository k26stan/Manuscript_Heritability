## Create Plot showing GO-FURTHER trial design ##
## February 12, 2016 ##
## Kristopher Standish ##

library(nlme)
library(gplots)

##############################################################
## LOAD DATA #################################################
##############################################################

## Set Date
DATE <- gsub("-","",Sys.Date())

## New Mac Paths
PathToRawFiles <- "/Users/kstandis/Data/Janssen/Data/Pheno/Raw_Files/"
PathToData <- "/Users/kstandis/Data/Janssen/Data/Pheno/Derived/20160112_Resp_v_Time.txt"
PathToFT <- "/Users/kstandis/Data/Janssen/Data/Pheno/Derived/20151015_Full_Table.txt"
PathToDer <- "/Users/kstandis/Data/Janssen/Data/Pheno/Derived/20150619_Derived_Pheno_Table.txt"
# PathToPlot <- paste("/Users/kstandis/Data/Janssen/Plots_Mac/",DATE,"_Phenotype_Plots/",sep="" )
PathToPlot <- paste("/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Plots/",DATE,"_TrialDesign/",sep="" )
dir.create( PathToPlot )

## Previously Compiled Data
FT <- read.table( PathToFT,sep="\t",header=T )
DER <- read.table( PathToDer,sep="\t",header=T )
DAT.l <- read.table( PathToData,sep="\t",header=T )

## Get Weeks of Measurement
WKS <- as.numeric(unique( DAT.l$WK ))
N.WKS <- length(WKS)

##############################################################
## FILTER DATA ###############################################
##############################################################

## Take out Patients who left within 4 weeks of getting DRUG
RM.exit.id <- as.character( FT$ID_2[which(FT$IN<=4)] )
RM.exit <- which( DAT.l$IID %in% RM.exit.id )
DAT <- DAT.l[-RM.exit,c(1:15,17)]

## Remove NA Values for DAS
DAT <- DAT[ which(!is.na(DAT$DAS)), ]
dim(DAT)

##############################################################
## PLOT TRIAL DESIGN (2) #####################################
##############################################################

####################################
## Parameters ######################

## Weeks of Injections
WKS.inj <- c(0,4,12,16,20,24,28+(0:9)*8)
INJ <- list()
INJ$GOL <- c(1,2)[factor(WKS.inj %in% c(16,24))]
INJ$EE <- c(1,2)[factor(WKS.inj %in% c(0:15,24))]
INJ$PBO <- c(1,2)[factor(WKS.inj %in% c(0:20))]
ARMS <- names(INJ)
LWD.inj <- 8

## Arm Labels
LABS <- c("A","B2","B1")
LABS <- c("GOL","PBO-EE","PBO-NE")
names(LABS) <- ARMS

## Point Shape/Size
PCHS <- c(1,4,16)
PCHS <- c(16,4)
CEX.pt <- 3
LWD.pt <- 12
LWD.ln <- 10

## WAG Labels
WAGS <- list()
WAGS$GOL <- c(0,4,12,20,28)
WAGS$PBO <- 24+c(0,4,12,20,28)
WAGS$EE <- 16+c(0,4,12,20,28)
CEX.wag <- 1

## Timeline Y values
GAP <- .22
yvals <- c(.5,2)
YVALS <- list()
YVALS$GOL <- rep(yvals[2],N.WKS)
YVALS$PBO <- rep( c(yvals[1],yvals[2]-2*GAP),c(8,8))
YVALS$EE <- rep( c(yvals[1]+GAP,yvals[2]-GAP),c(6,10))
for ( y in names(YVALS) ) { names(YVALS[[y]]) <- WKS }

## Plot Colors
COLS <- c("seagreen3","chocolate1","mediumpurple2") ; names(COLS) <- ARMS
COLS <- colorRampPalette(c("seagreen3","black"))(4)[1:3] ; names(COLS) <- ARMS
COLS <- c("seagreen3","brown3","goldenrod2") ; names(COLS) <- ARMS
COLS <- c("seagreen3","deeppink2","royalblue2") ; names(COLS) <- ARMS
COLS.ann <- c("chocolate2","mediumpurple2","gold2","steelblue2","cadetblue3","firebrick1")
XLIM <- c(-6,100)
YLIM <- c(0,3.25)

####################################
## Create Plot #####################

## Open File
png( paste(PathToPlot,"1-TrialDesign.png",sep=""),height=1200,width=2000,pointsize=36 )
par(mar=c(4,1,4,1))

## Create Plot
plot( 0,0,type="n",xlim=XLIM,ylim=YLIM,xaxt="n",yaxt="n",xlab="Week",ylab="",main="GO-FURTHER Trial Design")
axis( 1,at=WKS,las=2 )
abline( v=WKS, lty=3,col="black",lwd=1 )
## Label Time Point
 # Baseline
abline( v=c(0,16,24),col=COLS.ann[2],lty=2,lwd=LWD.inj )
text( 3+c(0,16-1,24),2.75, label=c("Baseline\nMeasurements","Early Escape","Non-Early\nCrossover"),srt=90,col=COLS.ann[2] )
## Plot ARMs of Trial Design
 # Trial Arms
SCRAP <- lapply( ARMS, function(x)points(WKS.inj,YVALS[[x]][as.character(WKS.inj)],col=COLS[[x]],pch=PCHS[INJ[[x]]],lwd=LWD.pt,cex=CEX.pt ))
SCRAP <- lapply( ARMS, function(x)points(WKS,YVALS[[x]],col=COLS[[x]],type="l",lwd=LWD.ln ))
## Annotate Plot
 # Initial Randomization
if ( XLIM[1]<0 ) {
	TEMP_YVALS <- c( YVALS$GOL[1], mean(c(YVALS$PBO[1],YVALS$EE[1])) )
	arrows( rep(XLIM[1]-1,2), rep(mean(TEMP_YVALS),2), rep(-2,2), TEMP_YVALS+c(-GAP,GAP), lwd=c(2,1)*7,length=.3 )
	text( XLIM[1]+2,mean(TEMP_YVALS), label="2:1\nRandomization",pos=4,cex=1.2 )
	text( XLIM[1]-1, TEMP_YVALS, label=c("GOL","PBO"),cex=1.2 )
	points( XLIM[1]-1, mean(TEMP_YVALS), pch=20,cex=2 )
}
 # Indicate Observations at Bottom of Plot
arrows( WKS,rep(.15,N.WKS),WKS,rep(-.1,N.WKS),col=COLS.ann[1],lwd=LWD.inj,angle=30,length=.3 )
text( XLIM[2],.25,label="Clinical Observations",col=COLS.ann[1],pos=2)
 # Legends
legend( 66,YLIM[2],legend=LABS[c(1,3,2)],lty=1,lwd=LWD.ln,col=COLS[ARMS][c(1,3,2)],title="Group")
legend( 88,YLIM[2],legend=c("GOL","PBO"),pch=PCHS,title="Injection",pt.lwd=LWD.pt,pt.cex=.8*CEX.pt)
 # WAG
SCRAP <- lapply( ARMS, function(x)text( WAGS[[x]], max(YVALS[[x]]), label=WAGS$GOL,cex=CEX.wag) )
SCRAP <- lapply( ARMS, function(x)points( WAGS[[x]], YVALS[[x]][as.character(WAGS[[x]])], pch=1,col="black", cex=CEX.pt,lwd=4 ) )
text( tail(WAGS$PBO,1), tail(YVALS$PBO,1)-GAP, label="(WAG)", cex=CEX.wag )
dev.off()




##############################################################
## END OF DOC ################################################
##############################################################















# ##############################################################
# ## PLOT TRIAL DESIGN (1) #####################################
# ##############################################################
# ARM <- list()
# ARM$GOL <- rep( c(1,3),c(1,15) )
# ARM$PBO <- rep( 1:3,c(1,8,7))
# ARM$EE <- rep( 1:3,c(1,6,9))
# # ARM$GOL <- rep( 1,16 )
# # ARM$PBO <- rep( 2:1,c(8,8))
# # ARM$EE <- rep( 2:1,c(6,10))
# ARMS <- names(ARM)
# PCHS <- c(1,4,16)

# YVALS <- list()
# YVALS$GOL <- rep(3,N.WKS)
# YVALS$PBO <- rep(2,N.WKS)
# YVALS$EE <- rep(2:1,c(7,9))

# COLS <- c("seagreen3","chocolate1","mediumpurple2") ; names(COLS) <- ARMS
# COLS.ann <- c("steelblue2","cadetblue3","gold2","firebrick1")
# XLIM <- c(0,110)
# YLIM <- c(0,4.25)
# ## Open Plot
# png( paste(PathToPlot,"2-TrialDesign.png",sep=""),height=1400,width=2000,pointsize=36 )
# plot( 0,0,type="n",xlim=XLIM,ylim=YLIM,xaxt="n",yaxt="n",xlab="Week",ylab="",main="GO-FURTHER Trial Design")
# axis( 1,at=WKS )
# abline( v=WKS, lty=3,col="black",lwd=1 )
# ## Plot ARMs of Trial Design
#  # Initial Randomization
# if ( XLIM[1]<0 ) {
# 	points( c(0,XLIM[1],0),c(2,2.5,3),col="black",lwd=6,type="l" )
# 	text( XLIM[1],2.5, label="2:1 Randomization",pos=4 )
# }
#  # Trial Arms
# lapply(ARMS,function(x)points(WKS,YVALS[[x]],col=COLS[[x]],pch=PCHS[ARM[[x]]],type="o",lwd=6,cex=2 ))
# text( 100, 3:1, label=ARMS, col=COLS, pos=4 )
# ## Annotate Plot
#  # Baseline
# abline( v=0,col=COLS.ann[1],lty=2,lwd=4 )
# text( 3,3.75, label="Baseline\nMeasurements",srt=90,col=COLS.ann[1] )
#  # Observations
# arrows( WKS,rep(.15,N.WKS),WKS,rep(-.1,N.WKS),col=COLS.ann[2],lwd=6,angle=30,length=.3 )
# text( XLIM[2],.25,label="Clinical Observations",col=COLS.ann[2],pos=2)
#  # PBO/EE Crossovers
# arrows( WKS[c(7,9)],c(2.3,1.7),WKS[c(7,9)],c(2.1,1.9),col=COLS.ann[4],lwd=4,angle=30,length=.3 )
# text( WKS[c(7,9)]+4,c(2.4,1.6), label=c("EE Crossover","PBO Crossover"),col=COLS.ann[4] )
#  # Legend
# legend("topright",legend=c("Baseline","On PBO","On GOL"),pch=PCHS,title="Clinical Observation:",pt.lwd=6,pt.cex=2)
# # legend("topright",legend=c("GOL","PBO"),pch=PCHS,title="Injection",pt.lwd=6,pt.cex=2)
# dev.off()








# ##############################################################
# ## PLOT DAS vs TIME ##########################################
# ##############################################################

# ## Specify Inputs
# N_SAMP.1 <- 20
# N_SAMP.2 <- 3
# TAG <- "Test"

# ## Create Function to Plot DAS (for some patients) vs Time
# DAS_v_TIME <- function( DAT, N_SAMP.1, N_SAMP.2, FT, TAG ) {
# 	Samps.all <- as.character(unique( FT$ID_2 ))
# 	Samps.which <- sample( Samps.all, N_SAMP.1, replace=F )
# 	if (N_SAMP.2==3) {
# 		Samps.bold <- unlist(lapply(unique(FT$GRP),function(x)sample(as.character(unique(FT$ID_2[FT$GRP==x])),1) ))
# 	}else{
# 		Samps.bold <- sample( Samps.which, N_SAMP.2, replace=F )	
# 	}
	
# 	## Plotting Parameters
# 	XLIM <- c(0,100)
# 	YLIM <- c(0,9)
# 	 # Plot Colors
# 	# COLS <- c("deepskyblue1","chartreuse1","firebrick1")
# 	COLS <- c("seagreen3","chocolate1","mediumpurple2")
# 	names(COLS) <- c("G","P","PE")
# 	PCHS <- c(1,4,16)
# 	names(PCHS) <- c("00","01","10")
# 	## Save to File Location
# 	png( paste(PathToPlot,"1b_DAS_v_Time.",TAG,".png",sep=""), height=1000,width=1600,pointsize=28 )
# 	## Open Plot
# 	plot( 0,0,type="n", xlim=XLIM,xaxt="n",ylim=YLIM,xlab="Week",ylab="DAS", main="Disease Activity Score vs Time\nSampling of Patients" )
# 	axis( 1, at=WKS )
# 	abline( h=seq(0,10,1), lty=3,col="grey50",lwd=1 )
# 	abline( v=WKS, lty=3,col="grey50",lwd=1)
# 	## Sample Data
# 	for ( s in 1:N_SAMP.1 ) {
# 		samp <- Samps.which[s]
# 		COLOR <- COLS[ FT[which(FT$ID_2==samp),"GRP"] ]
# 		COLOR.4 <- gsub("1","4",COLOR)
# 		COLOR.4 <- colorRampPalette(c(COLOR,"black"))(5)[3] 
# 		WHICH <- which(DAT$IID==samp)
# 		XVALS <- DAT$WK[WHICH]
# 		YVALS <- DAT$DAS[WHICH]
# 		pchs <- PCHS[paste(DAT$DRUG[WHICH],DAT$PLAC[WHICH],sep="")]# factor(DAT$DRUG[WHICH])]
# 		points( XVALS, YVALS, pch=pchs,type="o", col="grey85", lwd=2 )
# 	}
# 	abline( v=c(0,24,16), lty=3,col=COLS,lwd=3)
# 	for ( s in 1:N_SAMP.2 ) {
# 		samp <- Samps.bold[s]
# 		COLOR <- COLS[ FT[which(FT$ID_2==samp),"GRP"] ]
# 		COLOR.4 <- colorRampPalette(c(COLOR,"black"))(5)[3] # gsub("1","4",COLOR)
# 		WHICH <- which(DAT$IID==samp)
# 		XVALS <- DAT$WK[WHICH]
# 		YVALS <- DAT$DAS[WHICH]
# 		pchs <- PCHS[paste(DAT$DRUG[WHICH],DAT$PLAC[WHICH],sep="")]
# 		points( XVALS, YVALS, pch=pchs,type="o", col=COLOR, lwd=4,lty=1,cex=2 )
# 	}
# 	dev.off()
# }
# # DAS_v_TIME( DAT, N_SAMP.1, N_SAMP.2, FT, "Test" )

# DAS_v_TIME( DAT, 50, 3, FT, "50" )

# ##############################################################
# ## COVARIATE PLOT ############################################
# ##############################################################

# ## Show reason for Covariate Correction in GWAS
#  # Save to File
# png( paste(PathToPlot,"2_DEL_v_Init.png",sep=""), height=1000,width=2000,pointsize=32 )
#  # 2 Plots in 1 Window
# par(mfrow=c(1,2))
#  # Plotting Colors
# COLS.list <- c("firebrick1","chocolate1","gold2","springgreen2","steelblue2","slateblue3")
# COLS <- colorRampPalette(COLS.list)(nrow(FT))
# # COLS <- colorRampPalette( rep(COLS.list,2) )(nrow(FT))
# COLS.init <- COLS[ rank(as.numeric(FT$DEL_MNe_MN)) ]

# ## Plot Delta-DAS vs Initial DAS
#  # Plot Limits
# XLIM <- range( FT$DAS_BL_MN )
# YLIM <- range( FT$DEL_MNe_MN, na.rm=T )
#  # Open First Plot
# plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, xlab="Pre-Treatment DAS",ylab="Change in DAS",main="Improvement vs Initial Disease State" )
#  # Lines
# abline( h=seq(-10,10,1),lty=3,col="grey50",lwd=1 )
# abline( v=seq(0,10,1),lty=3,col="grey50",lwd=1 )
#  # Plot Data & Line of Fit
# points( FT$DAS_BL_MN, FT$DEL_MNe_MN, col=COLS.init, pch="+" )
# MOD <- lm(FT$DEL_MNe_MN~FT$DAS_BL_MN)
# abline( MOD, col="black",lwd=4,lty=2 )
# ## Plot Residuals vs Initial DAS
# RESIDS <- resid(MOD)
#  # Plot Limits
# YLIM <- range( RESIDS, na.rm=T )
#  # Open Plot
# plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, xlab="Pre-Treatment DAS",ylab="Change in DAS",main="Residuals vs Initial Disease State" )
#  # Lines
# abline( h=seq(-10,10,1),lty=3,col="grey50",lwd=1 )
# abline( v=seq(0,10,1),lty=3,col="grey50",lwd=1 )
#  # Plot Data/ Line of Fit (y=0)
# points( RESIDS ~ FT$DAS_BL_MN[as.numeric(names(RESIDS))], col=COLS.init[as.numeric(names(RESIDS))], pch="+" )
# abline( h=0, col="black",lwd=4,lty=2 )

# dev.off()

# ##############################################################
# ## LONGITUDINAL ANALYSIS #####################################
# ##############################################################

# ##################################################
# ## Population Level Regression ###################

# ## Model Population Effects
# MOD.LM <- lm( DAS ~ DRUG*WK + PLAC, data=DAT )
#  # Simulate Patient from each Arm
# G.arm <- data.frame( WK=WKS,DRUG=c(0,rep(1,15)),PLAC=rep(0,16) )
# P.arm <- data.frame( WK=WKS,DRUG=c(rep(0,9),rep(1,7)),PLAC=c(0,rep(1,8),rep(0,7)) )
# PE.arm <- data.frame( WK=WKS,DRUG=c(rep(0,7),rep(1,9)),PLAC=c(0,rep(1,6),rep(0,9)) )
#  # Predict Disease State of Patient from Each Arm
# PRED.G <- predict( MOD.LM, newdata=G.arm )
# PRED.P <- predict( MOD.LM, newdata=P.arm )
# PRED.PE <- predict( MOD.LM, newdata=PE.arm )

# ## Plot Individuals over Time w/ Population Estimates
#  # All Samples
# Samps.all <- as.character(unique( FT$ID_2 ))
#  # Limits
# XLIM <- c(0,100)
# YLIM <- c(0,9)
#  # Plot Colors
# COLS <- c("deepskyblue1","chartreuse1","firebrick1")
# names(COLS) <- c("G","P","PE")
# COLS.4 <- gsub("1","4",COLS)
#  # Save to File Location
# png( paste(PathToPlot,"Long1_DAS_v_Time.Pop.png",sep=""), height=1000,width=1600,pointsize=28 )
#  # Open Plot
# plot( 0,0,type="n", xlim=XLIM,ylim=YLIM,xlab="Week",ylab="DAS", main="Disease Activity Score vs Time\nPopulation Estimates" )
# abline( h=seq(0,10,1), lty=3,col="grey50",lwd=1 )
# abline( v=WKS, lty=3,col="grey50",lwd=1)
# abline( v=c(0,24,16), lty=2,col=COLS,lwd=3)
#  # Sample Data
# for ( s in 1:length(Samps.all) ) {
# 	samp <- Samps.all[s]
# 	COLOR <- COLS[ FT[which(FT$ID_2==samp),"GRP"] ]
# 	COLOR.4 <- gsub("1","4",COLOR)
# 	WHICH <- which(DAT$IID==samp)
# 	XVALS <- DAT$WK[WHICH]
# 	YVALS <- DAT$DAS[WHICH]
# 	points( XVALS, YVALS, type="l", col="grey70", lwd=2 )
# }
#  # Plot Population Estimates
# points( XVALS, PRED.G, type="l", col=COLS.4["G"], lwd=4 )
# points( XVALS, PRED.P, type="l", col=COLS.4["P"], lwd=4 )
# points( XVALS, PRED.PE, type="l", col=COLS.4["PE"], lwd=4 )
# dev.off()

# ##################################################
# ## Patient Level Regression ######################

# ## Model Individual Effects
# DAT.LIS <- DAT[,c("IID","DAS","DRUG","WK","IID")]
# # DAT.LIS <- DAT.LIS[ -which(is.na(DAT.LIS),arr.ind=T)[,1], ]
# MOD.LIS <- lmList( DAS ~ DRUG*WK | IID, data=DAT.LIS )
# PRED.LIS <- predict( MOD.LIS )

# ## Plot Individuals over Time w/ a Few Individual Estimates
#  # All Samples
# Samps.all <- as.character(unique( FT$ID_2 ))
# N_FIT <- 5
# Samps.fit <- sample( Samps.all, N_FIT )
#  # Limits
# XLIM <- c(0,100)
# YLIM <- c(0,9)
#  # Plot Colors
# COLS <- c("deepskyblue1","chartreuse1","firebrick1")
# names(COLS) <- c("G","P","PE")
#  # Save to File Location
# png( paste(PathToPlot,"Long2_DAS_v_Time.Ind.png",sep=""), height=1000,width=1600,pointsize=28 )
#  # Open Plot
# plot( 0,0,type="n", xlim=XLIM,ylim=YLIM,xlab="Week",ylab="DAS", main="Disease Activity Score vs Time\nSampling of Patients" )
# abline( h=seq(0,10,1), lty=3,col="grey50",lwd=1 )
# abline( v=WKS, lty=3,col="grey50",lwd=1)
# abline( v=c(0,24,16), lty=2,col=COLS,lwd=3)
#  # Sample Data
# for ( s in 1:length(Samps.all) ) {
# 	samp <- Samps.all[s]
# 	COLOR <- COLS[ FT[which(FT$ID_2==samp),"GRP"] ]
# 	COLOR.4 <- gsub("1","4",COLOR)
# 	WHICH <- which(DAT.LIS$IID==samp)
# 	XVALS <- DAT.LIS$WK[WHICH]
# 	YVALS <- DAT.LIS$DAS[WHICH]
# 	points( XVALS, YVALS, type="l", col="grey70", lwd=2 )
# }
#  # Fit Data
# for ( s in 1:N_FIT ) {
# 	samp <- Samps.fit[s]
# 	COLOR <- COLS[ FT[which(FT$ID_2==samp),"GRP"] ]
# 	COLOR.4 <- gsub("1","4",COLOR)
# 	 # Raw Data
# 	WHICH <- which(DAT.LIS$IID==samp)
# 	XVALS <- DAT.LIS$WK[WHICH]
# 	YVALS <- DAT.LIS$DAS[WHICH]
# 	points( XVALS, YVALS, type="o", col=COLOR, lwd=2,lty=3, pch=16,cex=.8 )
# 	 # Fit Data
# 	YVALS <- PRED.LIS[ which(names(PRED.LIS)==samp) ]
# 	XVALS <- DAT.LIS$WK[ which(names(PRED.LIS)==samp) ]
# 	points( XVALS, YVALS, type="l", col=COLOR.4, lwd=4 )
# }
# dev.off()

# ## Plot Individual Estimates
# PAIRWISE_HIST <- function( COEFS, tag, color ) {
# 	COEF.N <- ncol(COEFS)
# 	COEF.name <- colnames(COEFS)
# 	COLOR <- color # "slateblue3" # "tomato2"
# 	COLOR2 <- "cadetblue2" # "seagreen3"
# 	 # Save to File Location
# 	png( paste(PathToPlot,"Long2_Coef_Hist.",tag,".png",sep=""), height=800+300*COEF.N,width=800+300*COEF.N,pointsize=36 )
# 	par(mfrow=c(COEF.N,COEF.N))
# 	 # 
# 	for ( r in 1:COEF.N ) {
# 		for ( c in 1:COEF.N ) {
# 			if ( r==c ) {
# 				hist( COEFS[,r], xlab=colnames(COEFS)[r], main=paste("Distribution -",colnames(COEFS)[r]), col=COLOR )
# 				abline( v=mean(COEFS[,r],na.rm=T), col=COLOR2,lty=2,lwd=5 )
# 			}
# 			if ( c > r ) {
# 				plot( 0,0, type="n", xlim=c(0,1),ylim=c(0,1), xlab="",ylab="", xaxt="n",yaxt="n" )
# 				text( 0.5,0.5, label=round(cor(COEFS[,c(r,c)],use="pairwise.complete.obs")[1,2],2), col=COLOR, cex=2 )
# 			}
# 			if ( r > c ) {
# 				plot( COEFS[,c], COEFS[,r], xlab=colnames(COEFS)[c], ylab=colnames(COEFS)[r], main=paste(colnames(COEFS)[c],"vs",colnames(COEFS)[r]), pch="+",col=COLOR )
# 				abline( v=mean(COEFS[,c],na.rm=T), col=COLOR2,lty=2,lwd=5 )
# 				abline( h=mean(COEFS[,r],na.rm=T), col=COLOR2,lty=2,lwd=5 )
# 			}
# 		}
# 	}
# 	dev.off()
# }
# PAIRWISE_HIST( coef(MOD.LIS), "LIS4", "slateblue3" )
# PAIRWISE_HIST( coef(MOD.LIS)[,1:3], "LIS3", "slateblue3" )

# ##################################################
# ## Mixed-Effects Modeling ########################

# DAT.LME <- DAT

# ## Mixed Effects Model
# MOD.LME <- lme( fixed = DAS ~ DRUG + PLAC, random = ~ DRUG | IID, data=DAT.LME )
# MOD.LME.0 <- lme( fixed = DAS ~ DRUG*WK + PLAC, random = ~ DRUG+WK | IID, data=DAT.LME )
# MOD.LME.30 <- lme( fixed = DAS ~ DRUG*I(WK-30) + PLAC, random = ~ DRUG+I(WK-30) | IID, data=DAT.LME )
# MOD.LME.50 <- lme( fixed = DAS ~ DRUG*I(WK-50) + PLAC, random = ~ DRUG+I(WK-50) | IID, data=DAT.LME )
# MOD.LME.cor <- lme( fixed = DAS ~ DRUG + PLAC, random = ~ DRUG | IID, data=DAT.LME, correlation=corCAR1(value = .5, form = ~1 | IID) )
# MOD.LME.cor.30 <- lme( fixed = DAS ~ DRUG*I(WK-30) + PLAC, random = ~ DRUG+I(WK-30) | IID, data=DAT.LME, correlation=corCAR1(value = .5, form = ~1 | IID) )
# MOD.LME <- MOD.LME.cor.30
# pairs( data.frame(coef(MOD.LME.0)[,1],coef(MOD.LME.30)[,1],coef(MOD.LME.50)[,1]) )
# MOD.LME <- MOD.LME.50

# ## Predict DAS Values for each Patient
# PRED.LME <- predict( MOD.LME )

# ## Plot Individuals over Time w/ a Few Individual Estimates
#  # All Samples
# Samps.all <- as.character(unique( FT$ID_2 ))
#  # Limits
# XLIM <- c(0,100)
# YLIM <- c(0,9)
#  # Plot Colors
# COLS <- c("deepskyblue1","chartreuse1","firebrick1")
# names(COLS) <- c("G","P","PE")
#  # Save to File Location
# png( paste(PathToPlot,"Long2_DAS_v_Time.LME.png",sep=""), height=1000,width=1600,pointsize=28 )
#  # Open Plot
# plot( 0,0,type="n", xlim=XLIM,ylim=YLIM,xlab="Week",ylab="DAS", main="Disease Activity Score vs Time\nSampling of Patients" )
# abline( h=seq(0,10,1), lty=3,col="grey50",lwd=1 )
# abline( v=WKS, lty=3,col="grey50",lwd=1)
# abline( v=c(0,24,16), lty=2,col=COLS,lwd=3)
#  # Sample Data
# for ( s in 1:length(Samps.all) ) {
# 	samp <- Samps.all[s]
# 	COLOR <- COLS[ FT[which(FT$ID_2==samp),"GRP"] ]
# 	COLOR.4 <- gsub("1","4",COLOR)
# 	WHICH <- which(DAT.LME$IID==samp)
# 	XVALS <- DAT.LME$WK[WHICH]
# 	YVALS <- DAT.LME$DAS[WHICH]
# 	points( XVALS, YVALS, type="l", col="grey70", lwd=2 )
# }
#  # Fit Data
# for ( s in 1:N_FIT ) {
# 	samp <- Samps.fit[s]
# 	COLOR <- COLS[ FT[which(FT$ID_2==samp),"GRP"] ]
# 	COLOR.4 <- gsub("1","4",COLOR)
# 	 # Raw Data
# 	WHICH <- which(DAT.LME$IID==samp)
# 	XVALS <- DAT.LME$WK[WHICH]
# 	YVALS <- DAT.LME$DAS[WHICH]
# 	points( XVALS, YVALS, type="o", col=COLOR, lwd=2,lty=3, pch=16,cex=.8 )
# 	 # Fit LIS Data
# 	YVALS <- PRED.LIS[ which(names(PRED.LIS)==samp) ]
# 	XVALS <- DAT.LIS$WK[ which(names(PRED.LIS)==samp) ]
# 	points( XVALS, YVALS, type="l", col=COLOR.4, lwd=3,lty=2 )
# 	 # Fit LME Data
# 	YVALS <- PRED.LME[ which(names(PRED.LME)==samp) ]
# 	XVALS <- DAT.LME$WK[ which(names(PRED.LME)==samp) ]
# 	points( XVALS, YVALS, type="l", col=COLOR.4, lwd=4 )
# }
# dev.off()

# ## Plot Mixed Effect Estimates
# PAIRWISE_HIST( coef(MOD.LME)[,1:3], "LME3", "chocolate1" )


# ##############################################################
# ## END OF DOC ################################################
# ##############################################################

# # ## Plot Residuals vs Initial DAS, PC1, PC2, & RF_ACPA
# # MOD <- lm( DEL_MNe_MN ~ DAS_BL_MN+PC1+PC2+RF_ACPA, data=FT )
# # RESIDS <- resid(MOD)
# #  # Plot Limits
# # YLIM <- range( RESIDS, na.rm=T )
# #  # Open Plot
# # plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, xlab="Pre-Treatment DAS",ylab="Change in DAS",main="Residuals vs Initial Disease State" )
# #  # Lines
# # abline( h=seq(-10,10,1),lty=3,col="grey50",lwd=1 )
# # abline( v=seq(0,10,1),lty=3,col="grey50",lwd=1 )
# #  # Plot Data/ Line of Fit (y=0)
# # points( RESIDS ~ FT$DAS_BL_MN[as.numeric(names(RESIDS))], col=COLS.init[as.numeric(names(RESIDS))], pch="+" )
# # abline( h=0, col="black",lwd=4,lty=2 )






