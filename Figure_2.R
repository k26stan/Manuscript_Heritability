## Make Figure 2 (Phenotype Distributions) for Resp_Herit Manuscript ##
## Before and After Transformations ##
## June 23, 2015 ##
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
PathToSing <- "Data/Burn/Data/Phenos/Full_Tables/20150520_Single_Pheno_Table.txt"
PathToSave <- paste("/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Plots/",DATE,sep="")
dir.create( PathToSave )

## Load Real Data
TAB <- read.table( PathToSing, sep="\t",header=T )
FT <- read.table( PathToFT, sep="\t",header=T )
RP <- read.table( PathToRep, sep="\t",header=T )

###########################################################
## MAKE FIGURE 2 ##########################################
###########################################################
## 4 rows, 5 columns
 # Rows
   # 1 - Distributions of Delta Values (raw)
   # 2 - Residuals (raw)
   # 3 - Residuals (transformed)
   # 4 - Residuals vs Initial Values (transformed)
 # Columns
   # 1 - DAS (red)
   # 2 - CRP (yellow)
   # 3 - SJC (green)
   # 4 - TJC (blue)
   # 5 - Tests

## Open File for Plot
png( paste(PathToSave,"/Figure_2.png",sep=""), height=2400,width=3000,pointsize=40 )
par(mfrow=c(4,5))

## Figure 2A ##############################################
## Distribution of Delta-Values (before Transformation) (week 20)

## Plotting Parameters
CATS <- c("DAS","CRP","SJC","TJC")
COLS <- c("firebrick1","gold2","chartreuse1","dodgerblue1")

## Compile P-Values
TIMES <- c("WAG4","WAG12","WAG20","WAG28","FL")
SHAP.p.a <- array(,c(5,4)) ; colnames(SHAP.p.a) <- CATS ; rownames(SHAP.p.a) <- TIMES
for ( time in TIMES ) {
	for ( cat in CATS ) {
		VALS <- TAB[,paste("DEL",time,cat,sep="_")]
		SHAP.p.a[time,cat] <- shapiro.test(VALS)$p.value
	}
}


 # DAS
cat <- CATS[1]
VALS <- TAB$DEL_WAG20_DAS
COLOR <- COLS[1]
BINSIZE <- .5
BRKS <- seq( floor(min(VALS,na.rm=T))-BINSIZE, ceiling(max(VALS,na.rm=T))+BINSIZE, BINSIZE )
TICKS.x <- BRKS[ which( round(BRKS,0)==BRKS ) ]
HIST <- hist( VALS, breaks=BRKS,col=COLOR,xaxt="n",main=paste("Change in",cat,"- 20WAG"),xlab=paste("Change in",cat) )
axis( 1, at=TICKS.x )
text( quantile(BRKS,.05),.9*max(HIST$counts), label=paste("Shapiro P\n",formatC(SHAP.p.a["WAG20",cat],digits=2,format="e")),col="black",pos=4 )
 # CRP
cat <- CATS[2]
VALS <- TAB$DEL_WAG20_CRP
COLOR <- COLS[2]
BINSIZE <- 10
BRKS <- seq( round(min(VALS,na.rm=T)-BINSIZE,-1), round(max(VALS,na.rm=T)+BINSIZE,-1), BINSIZE )
TICKS.x <- seq(-200,200,100)
HIST <- hist( VALS, breaks=BRKS,col=COLOR,xaxt="n",main=paste("Change in",cat,"- 20WAG"),xlab=paste("Change in",cat) )
axis( 1, at=TICKS.x )
text( quantile(BRKS,.05),.9*max(HIST$counts), label=paste("Shapiro P\n",formatC(SHAP.p.a["WAG20",cat],digits=2,format="e")),col="black",pos=4 )
 # SJC
cat <- CATS[3]
VALS <- TAB$DEL_WAG20_SJC
COLOR <- COLS[3]
BINSIZE <- 5
BRKS <- seq( round(min(VALS,na.rm=T)-BINSIZE,-1), round(max(VALS,na.rm=T)+BINSIZE,-1), BINSIZE )
TICKS.x <- seq(-100,100,20)
HIST <- hist( VALS, breaks=BRKS,col=COLOR,xaxt="n",main=paste("Change in",cat,"- 20WAG"),xlab=paste("Change in",cat) )
axis( 1, at=TICKS.x )
text( quantile(BRKS,.05),.9*max(HIST$counts), label=paste("Shapiro P\n",formatC(SHAP.p.a["WAG20",cat],digits=2,format="e")),col="black",pos=4 )
 # TJC
cat <- CATS[4]
VALS <- TAB$DEL_WAG20_TJC
COLOR <- COLS[4]
BINSIZE <- 5
BRKS <- seq( round(min(VALS,na.rm=T)-BINSIZE,-1), round(max(VALS,na.rm=T)+BINSIZE,-1), BINSIZE )
TICKS.x <- seq(-100,100,20)
HIST <- hist( VALS, breaks=BRKS,col=COLOR,xaxt="n",main=paste("Change in",cat,"- 20WAG"),xlab=paste("Change in",cat) )
axis( 1, at=TICKS.x )
text( quantile(BRKS,.05),.9*max(HIST$counts), label=paste("Shapiro P\n",formatC(SHAP.p.a["WAG20",cat],digits=2,format="e")),col="black",pos=4 )
 # Plot P-Values for each Timepoint
XLIM <- c( 1, length(CATS) )
YLIM <- c( 0, -log10(min(SHAP.p.a)) )
PCHS <- 1:nrow(SHAP.p.a) - 1
plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, xlab="Phenotype",ylab="-log10(p)",main="Shapiro Test",xaxt="n" )
axis( 1, at=1:length(CATS), label=CATS )
abline( h=seq(0,YLIM[2],5),lty=3,col="grey50",lwd=1 )
# legend( "topright", pch=PCHS,legend=TIMES,ncol=2, cex=.8 )
legend( "topleft", pch=PCHS,legend=TIMES,ncol=1, cex=.8 )
for ( c in 1:length(CATS) ) {
	cat <- CATS[c]
	points( rep(c,length(TIMES)), -log10(SHAP.p.a[,c]), col=COLS[c], pch=PCHS )
}

## Figure 2B ##############################################
## Distribution of Residuals (before Transformation) (week 20)

## Plotting Parameters
CATS <- c("DAS","CRP","SJC","TJC")
COLS <- c("firebrick1","gold2","chartreuse1","dodgerblue1")

## Calculate Residuals & Compile P-Values
TIMES <- c("WAG4","WAG12","WAG20","WAG28","FL")
SHAP.p.b <- array(,c(5,4)) ; colnames(SHAP.p.b) <- CATS ; rownames(SHAP.p.b) <- TIMES
for ( time in TIMES ) {
	for ( cat in CATS ) {
		if ( grepl("WAG",time) ) {
			VALS <- resid( lm(TAB[,paste("DEL",time,cat,sep="_")] ~ TAB[,paste("Ibl",cat,sep="_")] ))
		}else{
			VALS <- resid( lm(TAB[,paste("DEL",time,cat,sep="_")] ~ TAB[,paste("I0",cat,sep="_")] ))
		}
		SHAP.p.b[time,cat] <- shapiro.test(VALS)$p.value
	}
}

 # DAS
cat <- CATS[1]
VALS <- resid( lm(TAB[,paste("DEL",time,cat,sep="_")] ~ TAB[,paste("Ibl",cat,sep="_")] ))
COLOR <- COLS[1]
BINSIZE <- .5
BRKS <- seq( floor(min(VALS,na.rm=T))-BINSIZE, ceiling(max(VALS,na.rm=T))+BINSIZE, BINSIZE )
TICKS.x <- BRKS[ which( round(BRKS,0)==BRKS ) ]
HIST <- hist( VALS, breaks=BRKS,col=COLOR,xaxt="n",main=paste("Residual:",cat,"- 20WAG"),xlab=paste("Residual:",cat) )
axis( 1, at=TICKS.x )
text( quantile(BRKS,.05),.9*max(HIST$counts), label=paste("Shapiro P\n",formatC(SHAP.p.b["WAG20",cat],digits=2,format="e")),col="black",pos=4 )
 # CRP
cat <- CATS[2]
VALS <- resid( lm(TAB[,paste("DEL",time,cat,sep="_")] ~ TAB[,paste("Ibl",cat,sep="_")] ))
COLOR <- COLS[2]
BINSIZE <- 10
BRKS <- seq( round(min(VALS,na.rm=T)-BINSIZE,-1), round(max(VALS,na.rm=T)+BINSIZE,-1), BINSIZE )
TICKS.x <- seq(-200,200,100)
HIST <- hist( VALS, breaks=BRKS,col=COLOR,xaxt="n",main=paste("Residual:",cat,"- 20WAG"),xlab=paste("Residual:",cat) )
axis( 1, at=TICKS.x )
text( quantile(BRKS,.05),.9*max(HIST$counts), label=paste("Shapiro P\n",formatC(SHAP.p.b["WAG20",cat],digits=2,format="e")),col="black",pos=4 )
 # SJC
cat <- CATS[3]
VALS <- resid( lm(TAB[,paste("DEL",time,cat,sep="_")] ~ TAB[,paste("Ibl",cat,sep="_")] ))
COLOR <- COLS[3]
BINSIZE <- 5
BRKS <- seq( round(min(VALS,na.rm=T)-BINSIZE,-1), round(max(VALS,na.rm=T)+BINSIZE,-1), BINSIZE )
TICKS.x <- seq(-100,100,20)
HIST <- hist( VALS, breaks=BRKS,col=COLOR,xaxt="n",main=paste("Residual:",cat,"- 20WAG"),xlab=paste("Residual:",cat) )
axis( 1, at=TICKS.x )
text( quantile(BRKS,.05),.9*max(HIST$counts), label=paste("Shapiro P\n",formatC(SHAP.p.b["WAG20",cat],digits=2,format="e")),col="black",pos=4 )
 # TJC
cat <- CATS[4]
VALS <- resid( lm(TAB[,paste("DEL",time,cat,sep="_")] ~ TAB[,paste("Ibl",cat,sep="_")] ))
COLOR <- COLS[4]
BINSIZE <- 5
BRKS <- seq( round(min(VALS,na.rm=T)-BINSIZE,-1), round(max(VALS,na.rm=T)+BINSIZE,-1), BINSIZE )
TICKS.x <- seq(-100,100,20)
HIST <- hist( VALS, breaks=BRKS,col=COLOR,xaxt="n",main=paste("Residual:",cat,"- 20WAG"),xlab=paste("Residual:",cat) )
axis( 1, at=TICKS.x )
text( quantile(BRKS,.05),.9*max(HIST$counts), label=paste("Shapiro P\n",formatC(SHAP.p.b["WAG20",cat],digits=2,format="e")),col="black",pos=4 )
 # Plot P-Values for each Timepoint
XLIM <- c( 1, length(CATS) )
YLIM <- c( 0, -log10(min(SHAP.p.b)) )
PCHS <- 1:nrow(SHAP.p.b) - 1
plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, xlab="Phenotype",ylab="-log10(p)",main="Shapiro Test",xaxt="n" )
axis( 1, at=1:length(CATS), label=CATS )
abline( h=seq(0,YLIM[2],5),lty=3,col="grey50",lwd=1 )
# legend( "topright", pch=PCHS,legend=TIMES,ncol=2, cex=.8 )
for ( c in 1:length(CATS) ) {
	cat <- CATS[c]
	points( rep(c,length(TIMES)), -log10(SHAP.p.b[,c]), col=COLS[c], pch=PCHS )
}

## Figure 2C ##############################################
## Distribution of Residuals (after Transformation) (week 20)

## Plotting Parameters
CATS <- c("DAS","lCRP","rSJC","rTJC")
COLS <- c("firebrick1","gold2","chartreuse1","dodgerblue1")

## Calculate Residuals & Compile P-Values
TIMES <- c("WAG4","WAG12","WAG20","WAG28","FL")
SHAP.p.c <- array(,c(5,4)) ; colnames(SHAP.p.c) <- CATS ; rownames(SHAP.p.c) <- TIMES
for ( time in TIMES ) {
	for ( cat in CATS ) {
		if ( grepl("WAG",time) ) {
			VALS <- resid( lm(TAB[,paste("DEL",time,cat,sep="_")] ~ TAB[,paste("Ibl",cat,sep="_")] ))
		}else{
			VALS <- resid( lm(TAB[,paste("DEL",time,cat,sep="_")] ~ TAB[,paste("I0",cat,sep="_")] ))
		}
		SHAP.p.c[time,cat] <- shapiro.test(VALS)$p.value
	}
}

 # DAS
cat <- CATS[1]
VALS <- resid( lm(TAB[,paste("DEL",time,cat,sep="_")] ~ TAB[,paste("Ibl",cat,sep="_")] ))
COLOR <- COLS[1]
BINSIZE <- .5
BRKS <- seq( floor(min(VALS,na.rm=T))-BINSIZE, ceiling(max(VALS,na.rm=T))+BINSIZE, BINSIZE )
TICKS.x <- BRKS[ which( round(BRKS,0)==BRKS ) ]
HIST <- hist( VALS, breaks=BRKS,col=COLOR,xaxt="n",main=paste("Residual:",cat,"- 20WAG"),xlab=paste("Residual:",cat) )
axis( 1, at=TICKS.x )
text( quantile(BRKS,.05),.9*max(HIST$counts), label=paste("Shapiro P\n",formatC(SHAP.p.c["WAG20",cat],digits=2,format="e")),col="black",pos=4 )
 # CRP
cat <- CATS[2]
VALS <- resid( lm(TAB[,paste("DEL",time,cat,sep="_")] ~ TAB[,paste("Ibl",cat,sep="_")] ))
COLOR <- COLS[2]
BINSIZE <- .25
BRKS <- seq( round(min(VALS,na.rm=T)-BINSIZE,0), round(max(VALS,na.rm=T)+BINSIZE,0), BINSIZE )
TICKS.x <- seq(-10,10,1)
HIST <- hist( VALS, breaks=BRKS,col=COLOR,xaxt="n",main=paste("Residual:",cat,"- 20WAG"),xlab=paste("Residual:",cat) )
axis( 1, at=TICKS.x )
text( quantile(BRKS,.05),.9*max(HIST$counts), label=paste("Shapiro P\n",formatC(SHAP.p.c["WAG20",cat],digits=2,format="e")),col="black",pos=4 )
 # SJC
cat <- CATS[3]
VALS <- resid( lm(TAB[,paste("DEL",time,cat,sep="_")] ~ TAB[,paste("Ibl",cat,sep="_")] ))
COLOR <- COLS[3]
BINSIZE <- 1
BRKS <- seq( round(min(VALS,na.rm=T)-BINSIZE,0), round(max(VALS,na.rm=T)+BINSIZE,0), BINSIZE )
TICKS.x <- seq(-20,20,5)
HIST <- hist( VALS, breaks=BRKS,col=COLOR,xaxt="n",main=paste("Residual:",cat,"- 20WAG"),xlab=paste("Residual:",cat) )
axis( 1, at=TICKS.x )
text( quantile(BRKS,.05),.9*max(HIST$counts), label=paste("Shapiro P\n",formatC(SHAP.p.c["WAG20",cat],digits=2,format="e")),col="black",pos=4 )
 # TJC
cat <- CATS[4]
VALS <- resid( lm(TAB[,paste("DEL",time,cat,sep="_")] ~ TAB[,paste("Ibl",cat,sep="_")] ))
COLOR <- COLS[4]
BINSIZE <- 1
BRKS <- seq( round(min(VALS,na.rm=T)-BINSIZE,-1), round(max(VALS,na.rm=T)+BINSIZE,-1), BINSIZE )
TICKS.x <- seq(-20,20,5)
HIST <- hist( VALS, breaks=BRKS,col=COLOR,xaxt="n",main=paste("Residual:",cat,"- 20WAG"),xlab=paste("Residual:",cat) )
axis( 1, at=TICKS.x )
text( quantile(BRKS,.05),.9*max(HIST$counts), label=paste("Shapiro P\n",formatC(SHAP.p.c["WAG20",cat],digits=2,format="e")),col="black",pos=4 )
 # Plot P-Values for each Timepoint
XLIM <- c( 1, length(CATS) )
YLIM <- c( 0, -log10(min(SHAP.p.c)) )
PCHS <- 1:nrow(SHAP.p.c) - 1
plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, xlab="Phenotype",ylab="-log10(p)",main="Shapiro Test",xaxt="n" )
axis( 1, at=1:length(CATS), label=CATS )
abline( h=seq(0,YLIM[2],5),lty=3,col="grey50",lwd=1 )
# legend( "topright", pch=PCHS,legend=TIMES,ncol=2, cex=.8 )
for ( c in 1:length(CATS) ) {
	cat <- CATS[c]
	points( rep(c,length(TIMES)), -log10(SHAP.p.c[,c]), col=COLS[c], pch=PCHS )
}


## Figure 2D ##############################################
## Residuals (after Transformation) vs Initial Values (week 20)

## Plotting Parameters
CATS <- c("DAS","lCRP","rSJC","rTJC")
COLS <- c("firebrick1","gold2","chartreuse1","dodgerblue1")

## Calculate Residuals & Compile P-Values
TIMES <- c("WAG4","WAG12","WAG20","WAG28","FL")
BP.p.d <- array(,c(5,4)) ; colnames(BP.p.d) <- CATS ; rownames(BP.p.d) <- TIMES
for ( time in TIMES ) {
	for ( cat in CATS ) {
		if ( grepl("WAG",time) ) {
			MOD <- lm(TAB[,paste("DEL",time,cat,sep="_")] ~ TAB[,paste("Ibl",cat,sep="_")] )
		}else{
			MOD <- lm(TAB[,paste("DEL",time,cat,sep="_")] ~ TAB[,paste("I0",cat,sep="_")] )
		}
		BP.p.d[time,cat] <- bptest(MOD)$p.value
	}
}

 # DAS
cat <- CATS[1]
COLOR <- COLS[1]
Y.VALS <- resid( lm(TAB[,paste("DEL",time,cat,sep="_")] ~ TAB[,paste("Ibl",cat,sep="_")] ))
X.VALS <- TAB[,paste("Ibl",cat,sep="_")][as.numeric(names(Y.VALS))]
XLIM <- range(X.VALS) ; YLIM <- range(Y.VALS)
plot( X.VALS,Y.VALS, col=COLOR, main=paste("Residuals vs Initial",cat,"- 20WAG"),xlab=paste("Initial",cat),ylab="Residuals" )
abline( h=seq(-10,10,1),lty=3,col="grey50",lwd=1 ) ; abline( h=0,lty=1,col="black",lwd=1 )
text( quantile(XLIM,.05),quantile(YLIM,.9), label=paste("BP Test P\n",formatC(BP.p.d["WAG20",cat],digits=2,format="e")),col="black",pos=4 )
 # CRP
cat <- CATS[2]
COLOR <- COLS[2]
Y.VALS <- resid( lm(TAB[,paste("DEL",time,cat,sep="_")] ~ TAB[,paste("Ibl",cat,sep="_")] ))
X.VALS <- TAB[,paste("Ibl",cat,sep="_")][as.numeric(names(Y.VALS))]
XLIM <- range(X.VALS) ; YLIM <- range(Y.VALS)
plot( X.VALS,Y.VALS, col=COLOR, main=paste("Residuals vs Initial",cat,"- 20WAG"),xlab=paste("Initial",cat),ylab="Residuals" )
abline( h=seq(-10,10,1),lty=3,col="grey50",lwd=1 ) ; abline( h=0,lty=1,col="black",lwd=1 )
text( quantile(XLIM,.05),quantile(YLIM,.9), label=paste("BP Test P\n",formatC(BP.p.d["WAG20",cat],digits=2,format="e")),col="black",pos=4 )
 # SJC
cat <- CATS[3]
COLOR <- COLS[3]
Y.VALS <- resid( lm(TAB[,paste("DEL",time,cat,sep="_")] ~ TAB[,paste("Ibl",cat,sep="_")] ))
X.VALS <- TAB[,paste("Ibl",cat,sep="_")][as.numeric(names(Y.VALS))]
XLIM <- range(X.VALS) ; YLIM <- range(Y.VALS)
plot( X.VALS,Y.VALS, col=COLOR, main=paste("Residuals vs Initial",cat,"- 20WAG"),xlab=paste("Initial",cat),ylab="Residuals" )
abline( h=seq(-10,10,1),lty=3,col="grey50",lwd=1 ) ; abline( h=0,lty=1,col="black",lwd=1 )
text( quantile(XLIM,.05),quantile(YLIM,.9), label=paste("BP Test P\n",formatC(BP.p.d["WAG20",cat],digits=2,format="e")),col="black",pos=4 )
 # TJC
cat <- CATS[4]
COLOR <- COLS[4]
Y.VALS <- resid( lm(TAB[,paste("DEL",time,cat,sep="_")] ~ TAB[,paste("Ibl",cat,sep="_")] ))
X.VALS <- TAB[,paste("Ibl",cat,sep="_")][as.numeric(names(Y.VALS))]
XLIM <- range(X.VALS) ; YLIM <- range(Y.VALS)
plot( X.VALS,Y.VALS, col=COLOR, main=paste("Residuals vs Initial",cat,"- 20WAG"),xlab=paste("Initial",cat),ylab="Residuals" )
abline( h=seq(-10,10,1),lty=3,col="grey50",lwd=1 ) ; abline( h=0,lty=1,col="black",lwd=1 )
text( quantile(XLIM,.05),quantile(YLIM,.9), label=paste("BP Test P\n",formatC(BP.p.d["WAG20",cat],digits=2,format="e")),col="black",pos=4 )
 # Plot P-Values for each Timepoint
XLIM <- c( 1, length(CATS) )
YLIM <- c( 0, -log10(min(BP.p.d)) )
PCHS <- 1:nrow(BP.p.d) - 1
plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, xlab="Phenotype",ylab="-log10(p)",main="Breusch-Pagan Test",xaxt="n" )
axis( 1, at=1:length(CATS), label=CATS )
abline( h=seq(0,YLIM[2],5),lty=3,col="grey50",lwd=1 )
# legend( "topright", pch=PCHS,legend=TIMES,ncol=2, cex=.8 )
for ( c in 1:length(CATS) ) {
	cat <- CATS[c]
	points( rep(c,length(TIMES)), -log10(BP.p.d[,c]), col=COLS[c], pch=PCHS )
}

dev.off()


























###########################################################
## END OF DOC #############################################
###########################################################
