## Simulate Data for Repeated Measures ##
## January 5, 2015 ##
## Kristopher Standish ##

DATE <- "20150316"

PathToSave <- paste("/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Plots/",DATE,sep="")

###############################################
## 1 - Show Distributions for Groups ##########
###############################################

## Start by showing only 2 groups ##
MEAN <- 0
ST_DEV <- .5
GRP_DIFF <- .1
DEL_DAS_RANGE <- seq( MEAN-4*ST_DEV, MEAN+4*ST_DEV, .1 )
GRP_RANGE <- list()
GRP_RANGE$A <- dnorm( DEL_DAS_RANGE, MEAN, ST_DEV )
GRP_RANGE$B <- dnorm( DEL_DAS_RANGE, MEAN+GRP_DIFF, ST_DEV )

## Plot 2 distributions on same axes
SETS <- names(GRP_RANGE)
XLIM <- range(DEL_DAS_RANGE)
YLIM <- c(0,1)
COLS <- c( "black", "firebrick1","chocolate1","gold1","springgreen1","cadetblue1","steelblue2","slateblue3")
png( paste(PathToSave,"/3_SimRP-1_2_Groups.png",sep=""), height=1000,width=1400,pointsize=28)
plot( 0,0,type="n", xlim=XLIM, ylim=YLIM, main="Distribution of Mean Values", xlab="Delta-DAS", ylab="Frequency" )
abline( h=seq( 0,1,.1), lty=2, col="grey50" )
abline( v=seq( floor(XLIM[1]),ceiling(XLIM[2]),1), lty=2, col="grey50" )
for ( s in 1:length(SETS) ) {
	set <- SETS[s]
	lines( DEL_DAS_RANGE, GRP_RANGE[[set]], col=COLS[s], lwd=3, lty=1 )
	arrows( DEL_DAS_RANGE[which.max(GRP_RANGE[[set]])], 0, DEL_DAS_RANGE[which.max(GRP_RANGE[[set]])], max(GRP_RANGE[[set]])+.05, length=0, col=COLS[s], lwd=2 )
	text( DEL_DAS_RANGE[which.max(GRP_RANGE[[set]])], max(GRP_RANGE[[set]])+.075, labels=set, col=COLS[s])
}
dev.off()

## Plot >2 distributions on same axes
GRP_RANGE$C <- dnorm( DEL_DAS_RANGE, MEAN+2*GRP_DIFF, ST_DEV )
GRP_RANGE$D <- dnorm( DEL_DAS_RANGE, MEAN+3*GRP_DIFF, ST_DEV )
GRP_RANGE$E <- dnorm( DEL_DAS_RANGE, MEAN+4*GRP_DIFF, ST_DEV )
GRP_RANGE$F <- dnorm( DEL_DAS_RANGE, MEAN+5*GRP_DIFF, ST_DEV )
GRP_RANGE$G <- dnorm( DEL_DAS_RANGE, MEAN+6*GRP_DIFF, ST_DEV )
GRP_RANGE$H <- dnorm( DEL_DAS_RANGE, MEAN+7*GRP_DIFF, ST_DEV )
GRP_RANGE$I <- dnorm( DEL_DAS_RANGE, MEAN+8*GRP_DIFF, ST_DEV )
GRP_RANGE$J <- dnorm( DEL_DAS_RANGE, MEAN+9*GRP_DIFF, ST_DEV )
GRP_RANGE$K <- dnorm( DEL_DAS_RANGE, MEAN+10*GRP_DIFF, ST_DEV )
SETS <- names(GRP_RANGE)
XLIM <- range(DEL_DAS_RANGE)
YLIM <- c(0,1)
COLS.list <- c( "black", "firebrick1","chocolate1","gold1","springgreen1","cadetblue1","steelblue2","slateblue3")
COLS <- colorRampPalette(COLS.list)(length(SETS))
png( paste(PathToSave,"/3_SimRP-1_10_Groups.png",sep=""), height=1000,width=1400,pointsize=28)
plot( 0,0,type="n", xlim=XLIM, ylim=YLIM, main="Distribution of Mean Values", xlab="Delta-DAS", ylab="Frequency" )
abline( h=seq( 0,1,.1), lty=2, col="grey50" )
abline( v=seq( floor(XLIM[1]),ceiling(XLIM[2]),1), lty=2, col="grey50" )
for ( s in 1:length(SETS) ) {
	set <- SETS[s]
	lines( DEL_DAS_RANGE, GRP_RANGE[[set]], col=COLS[s], lwd=3, lty=1 )
	arrows( DEL_DAS_RANGE[which.max(GRP_RANGE[[set]])], 0, DEL_DAS_RANGE[which.max(GRP_RANGE[[set]])], max(GRP_RANGE[[set]])+.05, length=0, col=COLS[s], lwd=2 )
	text( DEL_DAS_RANGE[which.max(GRP_RANGE[[set]])], max(GRP_RANGE[[set]])+.075, labels=set, col=COLS[s])
}
dev.off()

###############################################
## 2 - Simulate Mean Values ###################
###############################################

## Simulate mean Values for 2 groups
SAMP_SIZE <- 200
MEAN <- 0
ST_DEV <- .5
GRP_DIFF <- .1
SIM_DAT <- list()
SIM_DAT$A <- rnorm( SAMP_SIZE, MEAN, ST_DEV )
SIM_DAT$B <- rnorm( SAMP_SIZE, MEAN+GRP_DIFF, ST_DEV )
GRP_RANGE <- list()
GRP_RANGE$A <- dnorm( DEL_DAS_RANGE, MEAN, ST_DEV )
GRP_RANGE$B <- dnorm( DEL_DAS_RANGE, MEAN+GRP_DIFF, ST_DEV )

## Plot it
SETS <- names(GRP_RANGE)
XLIM <- range(SIM_DAT,DEL_DAS_RANGE)
YLIM <- c(0, SAMP_SIZE/4)
COLS <- c( "grey30", "firebrick1","chocolate1","gold1","springgreen1","cadetblue1","steelblue1","slateblue1")
COLS.4 <- c( "black","firebrick3")
BRKS <- seq( floor(XLIM[1]),ceiling(XLIM[2]),.25 )
png( paste(PathToSave,"/3_SimRP-2_2grp_Samp.png",sep=""), height=1000,width=1400,pointsize=28)
 # Grp A
hist(SIM_DAT$A, breaks=BRKS, col=COLS.4[1],density=20,angle=45, xlim=XLIM, ylim=YLIM, main="Distribution of Mean Values", xlab="Delta-DAS", ylab="Frequency" )
abline( h=seq( 0,YLIM[2],10), lty=2, col="grey50" )
abline( v=seq( floor(XLIM[1]),ceiling(XLIM[2]),1), lty=2, col="grey50" )
hist(SIM_DAT$A, breaks=BRKS, col=COLS.4[1],density=20,angle=45, xlim=XLIM, ylim=YLIM, main="Distribution of Mean Values", xlab="Delta-DAS", ylab="Frequency", add=T )
hist(SIM_DAT$B, breaks=BRKS, col=COLS.4[2],density=20,angle=-45, xlim=XLIM, ylim=YLIM, main="Distribution of Mean Values", xlab="Delta-DAS", ylab="Frequency", add=T )
for ( s in 1:length(SETS) ) {
	set <- SETS[s]
	lines( DEL_DAS_RANGE, 50*GRP_RANGE[[set]], col=COLS[s], lwd=3, lty=1 )
	arrows( DEL_DAS_RANGE[which.max(GRP_RANGE[[set]])], 0, DEL_DAS_RANGE[which.max(GRP_RANGE[[set]])], 50*max(GRP_RANGE[[set]])+.05, length=0, col=COLS[s], lwd=2 )
	# text( DEL_DAS_RANGE[which.max(GRP_RANGE[[set]])], 50*max(GRP_RANGE[[set]])+.075, labels=set, col=COLS[s])
}
dev.off()

#####################################################
## 2- SIMULATE RESPONSE DATA ########################
#####################################################

## Assumptions ##
# There is a mean change from the initial DAS
# That mean change depends on which group somebody is in (A or B)
# There is variance associated with the mean change
# That variance is IID at each time point

## Specify Number of Individuals
N_IND <- 200

## Specify Difference b/n Groups
DIFF.GRP <- .1

## Specify Number of Time Points
N_TIME <- 10

## Create Function to Simulate Data
SIM.DAT <- function( N_IND, DIFF.GRP, N_TIME, TO_PLOT ) {
	## Specify Group Sizes for A & B
	N_IND.A <- N_IND.B <- N_IND

	## Specify Sample Distributions for Change in DAS
	DEL.DAS.MN.A <- 0 # mean( FT$DEL_MNe_MN+DIFF.GRP, na.rm=T )
	DEL.DAS.MN.B <- 0+DIFF.GRP # mean( FT$DEL_MNe_MN, na.rm=T )
	DEL.DAS.SD <- 0.5 # sd( FT$DEL_MNe_MN, na.rm=T )

	## Sample Delta-DAS Scores
	SIM.DEL.DAS.A <- rnorm( N_IND.A, DEL.DAS.MN.A, DEL.DAS.SD )
	SIM.DEL.DAS.B <- rnorm( N_IND.B, DEL.DAS.MN.B, DEL.DAS.SD )

	## Specify Sample Distributions for Variance
	# VAR.DAS.DIS <- sqrt(VAR.DAS)
	SIM.VAR.DEL.DAS.A <- 10*rbeta( N_IND.A, 1, 10 ) # sample( VAR.DAS.DIS, N_IND.A, replace=T )
	SIM.VAR.DEL.DAS.B <- 10*rbeta( N_IND.B, 1, 10 ) # sample( VAR.DAS.DIS, N_IND.B, replace=T )
	SIM.SD.DEL.DAS.A <- sqrt( SIM.VAR.DEL.DAS.A )
	SIM.SD.DEL.DAS.B <- sqrt( SIM.VAR.DEL.DAS.B )

	## Sample Multiple Timepoints for each Individual
	SIM.TIME.DAS.A <- array( , c(N_IND.A,N_TIME) )
	SIM.TIME.DAS.B <- array( , c(N_IND.B,N_TIME) )
	for ( i in 1:N_IND.A ) {
		SIM.TIME.DAS.A[i,] <- rnorm( N_TIME, SIM.DEL.DAS.A[i], SIM.SD.DEL.DAS.A[i] )
	}
	for ( i in 1:N_IND.B ) {
		SIM.TIME.DAS.B[i,] <- rnorm( N_TIME, SIM.DEL.DAS.B[i], SIM.SD.DEL.DAS.B[i] )
	}
	X <- 1:10
	SIM.TIME.DAS.A[X,X]
	SIM.TIME.DAS.B[X,X]

	## Plot a few Patients over Time
	if ( TO_PLOT==1 ) {
		N_PLOT <- 2
		COLS.LIST.A <- "black" # c("firebrick2","chocolate2","gold2")
		COLS.LIST.B <- "slateblue1" # c("springgreen2","steelblue2","slateblue3")
		COLS.A <- colorRampPalette(COLS.LIST.A)(N_PLOT)
		COLS.B <- colorRampPalette(COLS.LIST.B)(N_PLOT)
		WHICH_PLOT.A <- sample( 1:N_IND.A, N_PLOT )
		WHICH_PLOT.B <- sample( 1:N_IND.B, N_PLOT )
		XLIM <- c( 1,N_TIME )
		YLIM <- c(-3,3) # range( SIM.TIME.DAS.A, SIM.TIME.DAS.B )
		plot( 0,0,type="n", xlim=XLIM, ylim=YLIM, xlab="Time of Measurement",ylab="Delta-DAS", main="Simulated Change in DAS over Time")
		abline( h=seq(-10,10,1), lty=2, col="grey50" )
		for ( i in 1:N_PLOT ) {
			ind.A <- WHICH_PLOT.A[i]
			points( 1:N_TIME, SIM.TIME.DAS.A[ind.A,], col=COLS.A[i], pch=20, type="o", lty=2, lwd=1)
			ind.B <- WHICH_PLOT.B[i]
			points( 1:N_TIME, SIM.TIME.DAS.B[ind.B,], col=COLS.B[i], pch=20, type="o", lty=2, lwd=1)
		}
		abline( h=SIM.DEL.DAS.A[WHICH_PLOT.A], col=COLS.A, lty=1, lwd=2 )
		abline( h=SIM.DEL.DAS.B[WHICH_PLOT.B], col=COLS.B, lty=1, lwd=2 )
	}

	## Return Simulated Data
	COMPILE.PARAMS <- list( SIM.DEL.DAS.A, SIM.DEL.DAS.B, SIM.VAR.DEL.DAS.A, SIM.VAR.DEL.DAS.B )
	names(COMPILE.PARAMS) <- c("MN_A","MN_B","SD_A","SD_B")
	COMPILE.TIME <- list( SIM.TIME.DAS.A, SIM.TIME.DAS.B )
	names(COMPILE.TIME) <- c("Time_A","Time_B")
	COMPILE <- list( COMPILE.PARAMS, COMPILE.TIME )
	names(COMPILE) <- c("Par","Time")
	return(COMPILE)
} # Close "SIM.DAT" Function
TEST <- SIM.DAT( N_IND, DIFF.GRP, N_TIME, TO_PLOT=1 )

## Plot 2 for the slides
# png( paste(PathToSave,"/3_SimRP-4_DWAI.png",sep=""), height=1000,width=1000,pointsize=28 )
# TEST <- SIM.DAT(200, .1, 20, 1)
# dev.off()
# png( paste(PathToSave,"/3_SimRP-4_DWAI2.png",sep=""), height=1000,width=1000,pointsize=28 )
# TEST <- SIM.DAT(200, 1, 20, 1)
# dev.off()

#####################################################
## P-VAL CALCULATIONS ###############################
#####################################################

## Specify Number of Individuals
N_IND <- 200

## Specify Difference b/n Groups
DIFF.GRP <- .1

## Specify Number of Time Points
N_TIME <- 10

## Simulate Data
DAT.1 <- SIM.DAT( N_IND, DIFF.GRP, N_TIME, 1 )

## Calculate Test Statistics/Significance at each Time Point and w/ Mean
P.VALS <- numeric( N_TIME+2 )
P.VALS[1] <- t.test( DAT.1$Par$MN_A, DAT.1$Par$MN_B, paired=F, alternative="two.sided" )$p.value
P.VALS[2] <- t.test( rowMeans(DAT.1$Time$Time_A), rowMeans(DAT.1$Time$Time_B), paired=F, alternative="two.sided" )$p.value
for ( i in 1:N_TIME ) {
	P.VALS[i+2] <- t.test( DAT.1$Time$Time_A[,i], DAT.1$Time$Time_B[,i], paired=F, alternative="two.sided" )$p.value
}

## Plot P-Values for a pair of Groups
png( paste(PathToSave,"/3_SimRP-5_Test_AB.png",sep=""), height=1000,width=1400,pointsize=34 )
plot( 1:length(P.VALS), -log10(P.VALS), pch="+", col="firebrick2", main="Identifying Difference b/n Groups", xlab="Repeated Measures", ylab="-log10(p)", xaxt="n", type="b" )
abline( v=1:length(P.VALS), lty=2, col="grey50" )
abline( h=seq(0,10,.5), lty=2, col="grey50" )
abline( h=-log10(.05), lty=2, lwd=3, col="chartreuse" )
axis( 1, at=1:length(P.VALS), label=c("Par_Mean","Obs_Mean",3:length(P.VALS)-1), las=2 )
dev.off()

## Specify Difference b/n Groups
DIFF.GRP <- .2

## Simulate Data
DAT.1 <- SIM.DAT( N_IND, DIFF.GRP, N_TIME, 1 )

## Calculate Test Statistics/Significance at each Time Point and w/ Mean
P.VALS <- numeric( N_TIME+2 )
P.VALS[1] <- t.test( DAT.1$Par$MN_A, DAT.1$Par$MN_B, paired=F, alternative="two.sided" )$p.value
P.VALS[2] <- t.test( rowMeans(DAT.1$Time$Time_A), rowMeans(DAT.1$Time$Time_B), paired=F, alternative="two.sided" )$p.value
for ( i in 1:N_TIME ) {
	P.VALS[i+2] <- t.test( DAT.1$Time$Time_A[,i], DAT.1$Time$Time_B[,i], paired=F, alternative="two.sided" )$p.value
}

## Plot P-Values for a pair of Groups
png( paste(PathToSave,"/3_SimRP-5_Test_AC.png",sep=""), height=1000,width=1400,pointsize=34 )
plot( 1:length(P.VALS), -log10(P.VALS), pch="+", col="firebrick2", main="Identifying Difference b/n Groups", xlab="Repeated Measures", ylab="-log10(p)", xaxt="n", type="b" )
abline( v=1:length(P.VALS), lty=2, col="grey50" )
abline( h=seq(0,10,.5), lty=2, col="grey50" )
abline( h=-log10(.05), lty=2, lwd=3, col="chartreuse" )
axis( 1, at=1:length(P.VALS), label=c("Par_Mean","Obs_Mean",3:length(P.VALS)-1), las=2 )
dev.off()

#####################################################
## POWER CALCULATIONS ###############################
#####################################################

## Specify Number of Simulations
N_ITER <- 1000

## Specify Range of Sample Sizes
SAMP_SIZE <- seq( 20, 200, 20 )

## Specify Range of Differences b/n Groups
EFF_SIZE <- seq( .1, 1, .1 )

## Specify Number of Time Points
N_TIME <- 10

## Calculate Test Statistics/Significance at each Time Point and w/ Mean

P.VALS <- list()
## Loop Through Sample Sizes
start <- proc.time()
for ( s in 1:length(SAMP_SIZE) ) {
	size <- SAMP_SIZE[s]
	print(paste( "#### Running SAMPLE size:",size,"-",round( proc.time()-start, 2)[3] ))
	P.VALS[[s]] <- list()
	## Loop Thru Effect Sizes
	for ( e in 1:length(EFF_SIZE) ) {
		eff <- EFF_SIZE[e]
		P.VALS[[s]][[e]] <- array(, c(N_ITER,N_TIME+2) )
		## Sample a bunch of times
		for ( i in 1:N_ITER ) {
			## Simulate Data
			DAT <- SIM.DAT( size, eff, N_TIME, 0 )
			## Calculate P-Values
			P.VALS[[s]][[e]][i,1] <- t.test( DAT$Par$MN_A, DAT$Par$MN_B, paired=F, alternative="two.sided" )$p.value
			P.VALS[[s]][[e]][i,2] <- t.test( rowMeans(DAT$Time$Time_A), rowMeans(DAT$Time$Time_B), paired=F, alternative="two.sided" )$p.value
			for ( t in 1:N_TIME ) {
				P.VALS[[s]][[e]][i,t+2] <- t.test( DAT$Time$Time_A[,t], DAT$Time$Time_B[,t], paired=F, alternative="two.sided" )$p.value
			} # Close Iter Loop
		} # Close Effect Size Loop
		colnames(P.VALS[[s]][[e]]) <- c( "MN_Par","MN_Obs", paste("t",1:N_TIME,sep="") )
		rownames(P.VALS[[s]][[e]]) <- paste("i",1:N_ITER,sep="")
		print(paste( "Done with effect size:",eff,"-",round( proc.time()-start, 2)[3] ))
	} # Close Sample Size Loop
	names(P.VALS[[s]]) <- paste("e",EFF_SIZE,sep="")
}
names(P.VALS) <- paste("s",SAMP_SIZE,sep="")

## Go through Sample/Effect Sizes and Calculate Power
## Loop Through Sample Sizes
POWER.1 <- POWER.MN.1 <- POWER.MN.2 <- array( , c(length(SAMP_SIZE),length(EFF_SIZE)) )
rownames(POWER.1) <- rownames(POWER.MN.1) <- rownames(POWER.MN.2) <- names(P.VALS)
colnames(POWER.1) <- colnames(POWER.MN.1) <- colnames(POWER.MN.2) <- names(P.VALS[[1]])
start <- proc.time()
for ( s in 1:length(SAMP_SIZE) ) {
	size <- SAMP_SIZE[s]
	## Loop Thru Effect Sizes
	for ( e in 1:length(EFF_SIZE) ) {
		eff <- EFF_SIZE[e]
		POWER.MN.1[s,e] <- length(which( P.VALS[[s]][[e]][,1]<.05 )) / N_ITER
		POWER.MN.2[s,e] <- length(which( P.VALS[[s]][[e]][,2]<.05 )) / N_ITER
		POWER.1[s,e] <- length(which( P.VALS[[s]][[e]][,1:N_TIME+2]<.05 )) / N_ITER / N_TIME
	}
}

## Plot this Shiz
 # Set up Parameters
COLS.LIST <- c("maroon1","firebrick2","chocolate2","gold2","springgreen2","steelblue2","slateblue3","black")
COLS <- colorRampPalette(COLS.LIST)(length(EFF_SIZE))
COLS <- sample( colorRampPalette(COLS.LIST)(length(EFF_SIZE)) )
XLIM <- range( SAMP_SIZE ) + c(0,40)
YLIM <- c(0,1)
 # Plot 1 effect size
png( paste(PathToSave,"/3_SimRP-6_Power_2.png",sep=""), height=1000,width=1400,pointsize=34 )
plot( 0,0,type="n", xlim=XLIM, ylim=YLIM, xlab="Sample Size", ylab="Power", main="Power vs Sample/Effect Size" )
abline( h=seq( 0,1,.1), lty=2, col="grey50" )
for ( e in 1 ) {
	points( SAMP_SIZE, POWER.MN.1[,e], col=COLS[e], lty=1, type="o", lwd=3, pch=20 )
	points( SAMP_SIZE, POWER.MN.2[,e], col=COLS[e], lty=3, type="o", lwd=3, pch=20 )
	points( SAMP_SIZE, POWER.1[,e], col=COLS[e], lty=2, type="o", lwd=3, pch=20 )
}
legend( "bottomright", legend=c("Mean_Par","Mean_Obs","Single"), lty=c(1,3,2), lwd=3, cex=.8 )
dev.off()
 # Plot all effect sizes
png( paste(PathToSave,"/3_SimRP-6_Power_ALL.png",sep=""), height=1000,width=1400,pointsize=34 )
plot( 0,0,type="n", xlim=XLIM, ylim=YLIM, xlab="Sample Size", ylab="Power", main="Power vs Sample/Effect Size" )
abline( h=seq( 0,1,.1), lty=2, col="grey50" )
for ( e in 1:length(EFF_SIZE) ) {
	points( SAMP_SIZE, POWER.MN.1[,e], col=COLS[e], lty=1, type="o", lwd=3, pch=20 )
	points( SAMP_SIZE, POWER.MN.2[,e], col=COLS[e], lty=3, type="o", lwd=3, pch=20 )
	points( SAMP_SIZE, POWER.1[,e], col=COLS[e], lty=2, type="o", lwd=3, pch=20 )
}
legend( "bottomright", legend=c("Mean_Par","Mean_Obs","Single"), lty=c(1,3,2), lwd=3, cex=.8 )
legend( "topright", legend=EFF_SIZE, lty=1,col=COLS,title="Effect Sizes", lwd=3, cex=.8 )
dev.off()



# library(gplots)
# heatmap.2( POWER.MN.1, Rowv=F, Colv=F, dendrogram="none",scale="none",trace="none" )










#####################################################
## END OF DOC #######################################
#####################################################






############################################################################
## PART 2 ##################################################################
############################################################################

## Test out Variance Differences in Phenotypes ##
## Argue for using Mean rather than Single Timepoint ##
## January 2, 2015 ##
## Kristopher Standish ##

#####################################################
## LOAD DATA ########################################
#####################################################

## Set Paths
PathToData <- "/Users/kstandis/Data/Burn/Data/Phenos/Full_Tables/20141229_Full_Table.txt"
PathToSave <- paste("/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Plots/",DATE,sep="")
# PathToSave <- "/Users/kstandis/Dropbox/Schork/JNJ11/Slides/20150105_Q3_Deliver/"

## Load Data
FT <- read.table( PathToData, sep="\t", header=T )

#####################################################
## CALCULATE VARIANCE AFTER RESPONSE ################
#####################################################

## Specify Measurement Weeks
WKS <- c(0,2,4,8,12,14,16,20,24,28,36,44,52,68,84,100)

## Specify Arms/Groups: Gol, Plac, PlacEE
G <- which( FT$GRP=="G" )
P <- which( FT$GRP=="P" )
PE <- which( FT$GRP=="PE" )

## Specify Baseline Week for each Group
BL <- rep(0,nrow(FT))
BL[P] <- 24
BL[PE] <- 16

## Calculate Post-Treatment Variance for each Patient
VAR.DAS <- VAR.lCRP <- VAR.rSJC <- VAR.rTJC <- numeric( nrow(FT) )
for ( r in 1:nrow(FT) ) {
	start_wk <- BL[r]
	which_wks <- WKS[ (1+which(WKS==start_wk)):length(WKS) ]
	# DAS
	col_names <- paste( "DAS_",which_wks,"wk", sep="" )
	VAR.DAS[r] <- var( c(FT[r,col_names],recursive=T), na.rm=T )
	# lCRP
	col_names <- paste( "CRP_",which_wks,"wk", sep="" )
	VAR.lCRP[r] <- var( log10(c(FT[r,col_names],recursive=T)), na.rm=T )
	# rSJC
	col_names <- paste( "SJC_",which_wks,"wk", sep="" )
	VAR.rSJC[r] <- var( sqrt(c(FT[r,col_names],recursive=T)), na.rm=T )
	# rTJC
	col_names <- paste( "TJC_",which_wks,"wk", sep="" )
	VAR.rTJC[r] <- var( sqrt(c(FT[r,col_names],recursive=T)), na.rm=T )
}

## Plot Distribution of VAR.DAS
png( paste(PathToSave,"/3_SimRP-3_Hist_VAR_DAS.png",sep=""), height=1000,width=1000,pointsize=28 )
hist( VAR.DAS, breaks=seq(0,4,.25), col="purple", main="Distribution of Variance around Mean (Post-Gol)", xlab="Post-Treatment Variance" )
dev.off()

## Plot Correlations b/n Variance and Mean Post-Treatment DAS
par(mfrow=c(2,2))
 # DAS
plot( VAR.DAS ~ FT$DAS_PG_MNe )
MOD.DAS <- lm( VAR.DAS ~ FT$DAS_PG_MNe )
if ( anova(MOD.DAS)[1,5]<.05 ) { COL="red" }else{ COL="black" }
abline( MOD.DAS, col=COL )
 # lCRP
plot( VAR.lCRP ~ FT$lCRP_PG_MNe )
MOD.lCRP <- lm( VAR.lCRP ~ FT$lCRP_PG_MNe )
if ( anova(MOD.lCRP)[1,5]<.05 ) { COL="red" }else{ COL="black" }
abline( MOD.lCRP, col=COL )
 # rSJC
plot( VAR.rSJC ~ FT$rSJC_PG_MNe )
MOD.rSJC <- lm( VAR.rSJC ~ FT$rSJC_PG_MNe )
if ( anova(MOD.rSJC)[1,5]<.05 ) { COL="red" }else{ COL="black" }
abline( MOD.rSJC, col=COL )
 # rTJC
plot( VAR.rTJC ~ FT$rTJC_PG_MNe )
MOD.rTJC <- lm( VAR.rTJC ~ FT$rTJC_PG_MNe )
if ( anova(MOD.rTJC)[1,5]<.05 ) { COL="red" }else{ COL="black" }
abline( MOD.rTJC, col=COL )

## Plot Correlations b/n Variance and Mean Changed
par(mfrow=c(2,2))
 # DAS
plot( VAR.DAS ~ FT$DEL_MNe_MN )
MOD.DAS <- lm( VAR.DAS ~ FT$DEL_MNe_MN )
if ( anova(MOD.DAS)[1,5]<.05 ) { COL="red" }else{ COL="black" }
abline( MOD.DAS, col=COL )
 # lCRP
plot( VAR.lCRP ~ FT$DEL_lCRP_MNe_MN )
MOD.lCRP <- lm( VAR.lCRP ~ FT$DEL_lCRP_MNe_MN )
if ( anova(MOD.lCRP)[1,5]<.05 ) { COL="red" }else{ COL="black" }
abline( MOD.lCRP, col=COL )
 # rSJC
plot( VAR.rSJC ~ FT$DEL_rSJC_MNe_MN )
MOD.rSJC <- lm( VAR.rSJC ~ FT$DEL_rSJC_MNe_MN )
if ( anova(MOD.rSJC)[1,5]<.05 ) { COL="red" }else{ COL="black" }
abline( MOD.rSJC, col=COL )
 # rTJC
plot( VAR.rTJC ~ FT$DEL_rTJC_MNe_MN )
MOD.rTJC <- lm( VAR.rTJC ~ FT$DEL_rTJC_MNe_MN )
if ( anova(MOD.rTJC)[1,5]<.05 ) { COL="red" }else{ COL="black" }
abline( MOD.rTJC, col=COL )

## Plot Correlations b/n Variance and Residuals of Mean ~ Init
par(mfrow=c(2,2))
 # DAS
MOD.DAS.1 <- resid( lm( FT$DEL_MNe_MN ~ FT$DAS_BL_MN ) )
plot( VAR.DAS[as.numeric(names(MOD.DAS.1))] ~ MOD.DAS.1 )
MOD.DAS <- lm( VAR.DAS[as.numeric(names(MOD.DAS.1))] ~ MOD.DAS.1 )
if ( anova(MOD.DAS)[1,5]<.05 ) { COL="red" }else{ COL="black" }
abline( MOD.DAS, col=COL )
 # lCRP
MOD.lCRP.1 <- resid( lm( FT$DEL_lCRP_MNe_MN ~ FT$lCRP_BL_MN ) )
plot( VAR.lCRP[as.numeric(names(MOD.lCRP.1))] ~ MOD.lCRP.1 )
MOD.lCRP <- lm( VAR.lCRP[as.numeric(names(MOD.lCRP.1))] ~ MOD.lCRP.1 )
if ( anova(MOD.lCRP)[1,5]<.05 ) { COL="red" }else{ COL="black" }
abline( MOD.lCRP, col=COL )
 # rSJC
MOD.rSJC.1 <- resid( lm( FT$DEL_rSJC_MNe_MN ~ FT$rSJC_BL_MN ) )
plot( VAR.rSJC[as.numeric(names(MOD.rSJC.1))] ~ MOD.rSJC.1 )
MOD.rSJC <- lm( VAR.rSJC[as.numeric(names(MOD.rSJC.1))] ~ MOD.rSJC.1 )
if ( anova(MOD.rSJC)[1,5]<.05 ) { COL="red" }else{ COL="black" }
abline( MOD.rSJC, col=COL )
 # rTJC
MOD.rTJC.1 <- resid( lm( FT$DEL_rTJC_MNe_MN ~ FT$rTJC_BL_MN ) )
plot( VAR.rTJC[as.numeric(names(MOD.rTJC.1))] ~ MOD.rTJC.1 )
MOD.rTJC <- lm( VAR.rTJC[as.numeric(names(MOD.rTJC.1))] ~ MOD.rTJC.1 )
if ( anova(MOD.rTJC)[1,5]<.05 ) { COL="red" }else{ COL="black" }
abline( MOD.rTJC, col=COL )

