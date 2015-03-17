## Calculate Various Derived Statistics that may Represent Response ##
## March 2, 2015 ##
## Kristopher Standish ##

library( nlme )
library( gplots )

###########################################################
## LOAD DATA ##############################################
###########################################################

## Set Date
DATE <- "20150313"

## Set Paths to Data and to Save
PathToFT <- "/Users/kstandis/Data/Burn/Data/Phenos/Full_Tables/20141229_Full_Table.txt"
PathToRep <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150226_Resp_v_Time.txt"
PathToWAG <- "/Users/kstandis/Data/Burn/Data/Phenos/Full_Tables/20150310_Single_Pheno_Table.txt"
PathToSave <- paste("/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Plots/",DATE,sep="")

## Load Real Data
FT <- read.table( PathToFT, sep="\t",header=T )
RP <- read.table( PathToRep, sep="\t",header=T )
WAG <- read.table( PathToWAG, sep="\t",header=T )

## Game Plan
#

###########################################################
## ORGANIZE DATA ##########################################
###########################################################

## Set number of patients
Samps <- as.character( unique( RP$IID ) )
N.samps <- length( Samps )
WKS <- as.numeric( unique( RP$WK ) )

## Calculate Transformations of Phenotypes
 # DAS/lCRP/rSJC/rTJC/rSJC28/rTJC28
PH <- data.frame( lCRP=log10(RP$CRP), rSJC=sqrt(RP$SJC), rTJC=sqrt(RP$TJC), rSJC28=sqrt(RP$SJC28), rTJC28=sqrt(RP$TJC28) )
DAT.1 <- data.frame( RP, PH )
PH_COLS <- c(16,21:25)
names(PH_COLS) <- colnames(DAT.1)[PH_COLS]

## Remove Samples in Study for less than 8 weeks (as before)
RM.samps <- as.character( FT$ID_2[ which(FT$IN<8) ] )
which.RM.samps <- which( DAT.1$IID %in% RM.samps )
DAT.2 <- DAT.1[ -which.RM.samps, ]

## Remove subjects w/o any DRUG==1 timepoints
RM.samps <- c()
for ( s in 1:N.samps ) {
	samp <- Samps[s]
	TEMP_DAT <- DAT.2[ which(DAT.2$IID==samp), ]
	if ( length(which( TEMP_DAT$DRUG==1 ))==0 ) {
		RM.samps <- c( RM.samps, samp )
	}
}
which.RM.samps <- which( DAT.2$IID %in% RM.samps )
DAT.3 <- DAT.2[ -which.RM.samps, ]
## Remove subjects w/o all missing values for DRUG==1 timepoints
RM.samps <- c()
for ( s in 1:N.samps ) {
	samp <- Samps[s]
	TEMP_DAT <- DAT.3[ which(DAT.3$IID==samp), ]
	TEMP_DAT.1 <- TEMP_DAT[ which(TEMP_DAT$DRUG==1), PH_COLS[-1] ]
	if ( all(is.na(TEMP_DAT.1)) ) {
		RM.samps <- c( RM.samps, samp )
	}
}
which.RM.samps <- which( DAT.3$IID %in% RM.samps )
if ( length(which.RM.samps)>0 ) {
	DAT.4 <- DAT.3[ -which.RM.samps, ]
}else{ DAT.4 <- DAT.3 }

TAB <- DAT.4
TAB$IID <- as.character( TAB$IID )
###########################################################
## DERIVED STATS ##########################################
###########################################################
 # Mean (w/ LOCF)
 # Mean (w/o LOCF)
 # Mean (Area)
 # Cohen's D Stat
 # Beta Value
 # % Improvement
 # Trajectory
 # Variance after treatment

## Set number of patients
Samps <- as.character( unique( TAB$IID ) )
N.samps <- length( Samps )
WKS <- as.numeric( unique( TAB$WK ) )

###########################################################
## 1 - Mean (w/ LOCF)
D.MN.1.pre <- D.MN.1.post <- array( , c(N.samps,6) )
colnames(D.MN.1.pre) <- colnames(D.MN.1.post) <- names(PH_COLS)
for ( s in 1:N.samps ) {
	samp <- Samps[s]
	TEMP_TAB <- TAB[ which(TAB$IID==samp), ]
	if ( nrow(TEMP_TAB)<length(WKS) ) {
		wks_missing <- setdiff( WKS, TEMP_TAB$WK )
		TEMP_TAB.2 <- merge( data.frame(WK=WKS), TEMP_TAB, by="WK", all=T )
		for ( r in 2:length(WKS) ) {
			which_na <- which(is.na( TEMP_TAB.2[r,] ))
			if ( length( which_na )>0 ) {
				TEMP_TAB.2[ r, which_na ] <- TEMP_TAB.2[ r-1, which_na ]
			}
		}
		TEMP_TAB <- TEMP_TAB.2
	}
	D.MN.1.pre[s,] <- colMeans( TEMP_TAB[which(TEMP_TAB$DRUG==0),PH_COLS], na.rm=T )
	D.MN.1.post[s,] <- colMeans( TEMP_TAB[which(TEMP_TAB$DRUG==1),PH_COLS], na.rm=T )
}
## Check-up
D.MN.1.diff <- D.MN.1.post - D.MN.1.pre
length(which(is.na(D.MN.1.diff)))

###########################################################
## 2 - Mean (w/o LOCF)
D.MN.2.pre <- D.MN.2.post <- array( , c(N.samps,6) )
colnames(D.MN.2.pre) <- colnames(D.MN.2.post) <- names(PH_COLS)
for ( s in 1:N.samps ) {
	samp <- Samps[s]
	TEMP_TAB <- TAB[ which(TAB$IID==samp), ]
	D.MN.2.pre[s,] <- colMeans( TEMP_TAB[which(TEMP_TAB$DRUG==0),PH_COLS], na.rm=T )
	D.MN.2.post[s,] <- colMeans( TEMP_TAB[which(TEMP_TAB$DRUG==1),PH_COLS], na.rm=T )
}
## Check-up
D.MN.2.diff <- D.MN.2.post - D.MN.2.pre
length(which(is.na(D.MN.2.diff)))

###########################################################
## 3 - Mean (Area)
 # Trapezoid Rule: A = .5*h*(b1+b2)
   # h = WK[i+1] - WK[i]
   # b1 = Pheno[WK[i+1]]
   # b2 = Pheno[WK[i]]
D.MN.3.pre <- D.MN.3.post <- array( , c(N.samps,6) )
colnames(D.MN.3.pre) <- colnames(D.MN.3.post) <- names(PH_COLS)
for ( s in 1:N.samps ) {
	samp <- Samps[s]
	TEMP_TAB <- TAB[ which(TAB$IID==samp), ]
	which_0 <- which(TEMP_TAB$DRUG==0) ; n_0 <- length(which_0)
	which_1 <- which(TEMP_TAB$DRUG==1) ; n_1 <- length(which_1)
	if (n_0>1 ) {
		h_0 <- TEMP_TAB$WK[ which_0[2:n_0] ] - TEMP_TAB$WK[ which_0[2:n_0-1] ]
		b_0 <- TEMP_TAB[ which_0[2:n_0], PH_COLS ] + TEMP_TAB[ which_0[2:n_0-1], PH_COLS ]
		A_0 <- colSums( .5*h_0*b_0, na.rm=T )
		Sc_0 <- A_0 / ( TEMP_TAB$WK[which_0[n_0]] - TEMP_TAB$WK[which_0[1]])
	}else{
		b_0 <- TEMP_TAB[ which_0, PH_COLS ]
		A_0 <- c( b_0, recursive=T )
		Sc_0 <- A_0
	}
	if (n_1>1 ) {
		h_1 <- TEMP_TAB$WK[ which_1[2:n_1] ] - TEMP_TAB$WK[ which_1[2:n_1-1] ]
		b_1 <- TEMP_TAB[ which_1[2:n_1], PH_COLS ] + TEMP_TAB[ which_1[2:n_1-1], PH_COLS ]
		A_1 <- colSums( .5*h_1*b_1, na.rm=T )
		Sc_1 <- A_1 / ( TEMP_TAB$WK[which_1[n_1]] - TEMP_TAB$WK[which_1[1]])
	}else{
		b_1 <- TEMP_TAB[ which_1, PH_COLS ]
		A_1 <- c( b_1, recursive=T )
		Sc_1 <- A_1
	}
	D.MN.3.pre[s,] <- Sc_0
	D.MN.3.post[s,] <- Sc_1
}
## Check-up
D.MN.3.diff <- D.MN.3.post - D.MN.3.pre
length(which(is.na(D.MN.3.diff)))

###########################################################
## 4 - Mean: Cohen's D Statistic
 # ( mean(x1)-mean(x2) ) / s
 # s = sqrt( ( (n1-1)s1^2 + (n2-1)s2^2 ) / (n1+n2-2) )
D.MN.4.stat <- array( , c(N.samps,6) )
colnames(D.MN.4.stat) <- names(PH_COLS)
for ( s in 1:N.samps ) {
	samp <- Samps[s]
	TEMP_TAB <- TAB[ which(TAB$IID==samp), ]
	which_0 <- which(TEMP_TAB$DRUG==0) ; n_0 <- length(which_0)
	which_1 <- which(TEMP_TAB$DRUG==1) ; n_1 <- length(which_1)
	Mn_0 <- colMeans( TEMP_TAB[ which_0, PH_COLS ], na.rm=T )
	Mn_1 <- colMeans( TEMP_TAB[ which_1, PH_COLS ], na.rm=T )
	if ( n_0>1 ) {
		Sd_0 <- apply( TEMP_TAB[ which_0, PH_COLS ], 2, sd, na.rm=T )
	}else{ Sd_0 <- rep(0,length(Mn_0)) }
	if ( n_1>1 ) {
		Sd_1 <- apply( TEMP_TAB[ which_1, PH_COLS ], 2, sd, na.rm=T )
	}else{ Sd_1 <- rep(0,length(Mn_1)) }
	# Sd_1 <- apply( TEMP_TAB[ which_1, PH_COLS ], 2, sd )
	S <- sqrt( ( (n_1-1)*Sd_1^2 + (n_0-1)*Sd_0^2 ) / (n_0+n_1-2) )
	D <- ( Mn_1 - Mn_0 ) / S
	D[ which(S==0) ] <- ( Mn_1 - Mn_0 )[ which(S==0) ]
	D.MN.4.stat[s,] <- D
}
## Check-up
length(which(is.na(D.MN.4.stat)))
length(which(D.MN.4.stat=="Inf"))
# plot( D.MN.1.diff, D.MN.4.stat )

###########################################################
## 5/9 - Beta Value for Linear Model
 # Use lmList for individual Beta Values
D.MN.5.B <- D.MN.5.t <- D.MN.5.B.int <- D.MN.5.t.int <- array( , c(N.samps,6) )
colnames(D.MN.5.B) <- colnames(D.MN.5.t) <- colnames(D.MN.5.B.int) <- colnames(D.MN.5.t.int) <- names(PH_COLS)
D.MN.9.res <- array( , c(N.samps,6) )
colnames(D.MN.9.res) <- names(PH_COLS)
for ( p in 1:length(PH_COLS) ) {
	col <- PH_COLS[p]
	pheno <- names(PH_COLS)[p]
	TEMP_TAB <- TAB[ , c("IID","DRUG",pheno) ]
	colnames(TEMP_TAB)[3] <- "Pheno"
	TEMP_TAB <- TEMP_TAB[ which(!is.na(TEMP_TAB[,"Pheno"])) ,]
	MOD <- lmList( Pheno ~ DRUG | IID, data=TEMP_TAB )
	D.MN.5.B.int[,p] <- summary(MOD)$coefficients[,"Estimate",1]
	D.MN.5.t.int[,p] <- summary(MOD)$coefficients[,"t value",1]
	D.MN.5.B[,p] <- summary(MOD)$coefficients[,"Estimate",2]
	D.MN.5.t[,p] <- summary(MOD)$coefficients[,"t value",2]
	D.MN.9.res[,p] <- sapply( lapply( MOD, function(x) sd(resid(x)) ), "[",1)
}
length(which(is.na(D.MN.5.B)))

###########################################################
## 6 - % Improvement
 # Calculate Mean Difference (D.MN.2.diff), take % change from baseline
D.MN.6.perc <- D.MN.2.diff / D.MN.2.pre
D.MN.6.perc[,"lCRP"] <- ( 10^(D.MN.2.post[,"lCRP"])-10^(D.MN.2.pre[,"lCRP"]) ) / 10^(D.MN.2.pre[,"lCRP"])

## Check-up
length(which(is.na(D.MN.6.perc)))
# hist( D.MN.6.perc[,"lCRP"], breaks=seq(-1,40,.25) )
# hist( 10^(D.MN.2.diff[,"lCRP"]) )
# hist( 10^(D.MN.2.pre[,"lCRP"]) )
# plot( 10^(D.MN.2.pre[,"lCRP"]), 10^(D.MN.2.diff[,"lCRP"]) )
# TEMP <- which( D.MN.2.pre[,"lCRP"] < D.MN.2.diff[,"lCRP"] )
# TAB.6.lCRP <- data.frame(D.MN.2.pre[,2],D.MN.2.post[,2],D.MN.2.diff[,2],D.MN.6.perc[,2])
# pairs( TAB.6.lCRP )
# head( TAB.6.lCRP )
# TEMP <- which(D.MN.6.perc[,2]>0)
# TEMP <- which( D.MN.2.diff[,2] > D.MN.2.pre[,2] )
# TAB.6.lCRP[TEMP,]

###########################################################
## 7/8 - Trajectory: Beta Value for WK in Linear Model
 # Use lmList and table including only DRUG==1, and calculate trajectory
 # Also calculate residuals around WK fit
TAB.1 <- TAB[ which(TAB$DRUG==1), ]
D.MN.7.B.wk <- D.MN.7.B.int <- array( , c(N.samps,6) )
colnames(D.MN.7.B.wk) <- colnames(D.MN.7.B.int) <- names(PH_COLS)
D.MN.8.res <- array( , c(N.samps,6) )
colnames(D.MN.8.res) <- names(PH_COLS)
for ( p in 1:length(PH_COLS) ) {
	col <- PH_COLS[p]
	pheno <- names(PH_COLS)[p]
	TEMP_TAB.1 <- TAB.1[ , c("IID","WK",pheno) ]
	colnames(TEMP_TAB.1)[3] <- "Pheno"
	TEMP_TAB.1 <- TEMP_TAB.1[ which(!is.na(TEMP_TAB.1[,"Pheno"])) ,]
	MOD <- lmList( Pheno ~ WK | IID, data=TEMP_TAB.1 )
	D.MN.7.B.int[,p] <- sapply( lapply( MOD, function(x) coef(x)[1] ), "[",1)# summary(MOD)$coefficients[,"Estimate",1]
	D.MN.7.B.wk[,p] <- sapply( lapply( MOD, function(x) coef(x)[2] ), "[",1)
	D.MN.8.res[,p] <- sapply( lapply( MOD, function(x) sd(resid(x)) ), "[",1)
}
length(which(is.na(D.MN.7.B.int)))
length(which(is.na(D.MN.7.B.wk)))
D.MN.7.B.wk[ which(is.na(D.MN.7.B.wk)) ] <- 0
D.MN.8.res[ which(is.na(D.MN.8.res)) ] <- 0

###########################################################
## MAKE SOME PLOTS ########################################
###########################################################

## Sanity Check
# DATA <- D.MN.6.perc
# DATA <- D.MN.7.B.wk
# DATA <- D.MN.8.res
# par(mfrow=c(2,3))
# for ( col in 1:ncol(DATA) ) {
# 	BIN <- 0.1*diff(range(DATA[,col]))
# 	BRKS <- seq( min(DATA[,col]), max(DATA[,col])+BIN, BIN )
# 	XLIM <- range(BRKS)
# 	hist( DATA[,col], breaks=BRKS, xlim=XLIM, main=colnames(DATA)[col] )
# }

## Combine Differences into Single Table
DIFF_FRAME <- data.frame( D.MN.1.diff, D.MN.2.diff, D.MN.3.diff, D.MN.4.stat, D.MN.5.B, D.MN.6.perc, D.MN.7.B.wk, D.MN.8.res, D.MN.9.res )
DERIV_NAMES <- paste("DEL",c("MNw","MNwo","MNa","MNcd","Bdr","PRC","Bwk","VARdr","VARwk"),sep="_") ; N.Deriv <- length(DERIV_NAMES)
PHENO_NAMES <- names(PH_COLS) ; N.Pheno <- length(PHENO_NAMES)
FRAME_NAMES <- paste( rep(DERIV_NAMES,rep(N.Pheno,N.Deriv)), rep(PHENO_NAMES,N.Deriv), sep="_" )
colnames(DIFF_FRAME) <- FRAME_NAMES
rownames(DIFF_FRAME) <- Samps

## Combine Post-Treatment Scores into Single Table
D.MN.5.post <- D.MN.5.B.int+D.MN.5.B
POST_FRAME <- data.frame( D.MN.1.post, D.MN.2.post, D.MN.3.post, D.MN.5.post )
DERIV_NAMES <- paste("POST",c("MNw","MNwo","MNa","Bdr"),sep="_") ; N.Deriv <- length(DERIV_NAMES)
PHENO_NAMES <- names(PH_COLS) ; N.Pheno <- length(PHENO_NAMES)
FRAME_NAMES <- paste( rep(DERIV_NAMES,rep(N.Pheno,N.Deriv)), rep(PHENO_NAMES,N.Deriv), sep="_" )
colnames(POST_FRAME) <- FRAME_NAMES
rownames(POST_FRAME) <- Samps

## Combine Pre-Treatment Scores into Single Table
PRE_FRAME <- data.frame( D.MN.1.pre, D.MN.2.pre, D.MN.3.pre, D.MN.5.B.int )
DERIV_NAMES <- paste("PRE",c("MNw","MNwo","MNa","Bdr"),sep="_") ; N.Deriv <- length(DERIV_NAMES)
PHENO_NAMES <- names(PH_COLS) ; N.Pheno <- length(PHENO_NAMES)
FRAME_NAMES <- paste( rep(DERIV_NAMES,rep(N.Pheno,N.Deriv)), rep(PHENO_NAMES,N.Deriv), sep="_" )
colnames(PRE_FRAME) <- FRAME_NAMES
rownames(PRE_FRAME) <- Samps

## Heatmap of Correlation Amongst Phenotypes
COLS.list <- c("black","slateblue3","steelblue2","springgreen2","gold2","chocolate2","firebrick1")
COLS.list <- c("gold2","chocolate2","firebrick1","black","slateblue3","steelblue2","springgreen2")
COLS <- colorRampPalette(COLS.list)(100)
 # Delta-Stats
CORR.diff <- cor( DIFF_FRAME, use="pairwise.complete.obs", method="spearman" )
BRKS <- seq( -1,1,length.out=101 )
png( paste(PathToSave,"/4_Derived_Diff_Heatmap.png",sep=""), height=1600,width=1600,pointsize=32 )
heatmap.2( CORR.diff, main="Correlation: Derived Delta Stats", col=COLS, trace="none", breaks=BRKS, scale="none", margins=c(7,7) )
dev.off()
 # Pre/Post-Treatment Stats
CORR.pp <- cor( data.frame(PRE_FRAME,POST_FRAME), use="pairwise.complete.obs", method="spearman" )
BRKS <- seq( -1,1,length.out=101 )
png( paste(PathToSave,"/4_Derived_PrePost_Heatmap.png",sep=""), height=1600,width=1600,pointsize=32 )
heatmap.2( CORR.pp, main="Correlation: Derived Pre/Post Stats", col=COLS, trace="none", breaks=BRKS, scale="none", margins=c(7,7) )
dev.off()

###########################################################
## WRITE TABLES ###########################################
###########################################################

## Combine Derived Stats into one Table
FULL_TABLE <- data.frame( FID=rownames(DIFF_FRAME),IID=rownames(DIFF_FRAME), DIFF_FRAME, PRE_FRAME, POST_FRAME )

## Write Table
PathToWrite <- gsub("20150310_Single","20150313_Derived",PathToWAG)
write.table( FULL_TABLE, PathToWrite, sep="\t",col.names=T,row.names=F,quote=F)

## Write Phenotype/Covariate List
PHENO_NAMES.write <- colnames(FULL_TABLE)[3:ncol(FULL_TABLE)]
PHENO_NAMES.pre <- PHENO_NAMES.write[ grep("PRE",PHENO_NAMES.write) ]
PHENO_NAMES.write <- PHENO_NAMES.write[ grep("PRE",PHENO_NAMES.write,invert=T) ]
COV_NAMES.write <- PHENO_NAMES.write
COV_NAMES.write <- gsub("POST","PRE",COV_NAMES.write)
COV_NAMES.write[ grep("MNw_",COV_NAMES.write) ] <- COV_NAMES.write[ grep("PRE_MNw_",COV_NAMES.write) ]
COV_NAMES.write[ grep("MNwo_",COV_NAMES.write) ] <- COV_NAMES.write[ grep("PRE_MNwo_",COV_NAMES.write) ]
COV_NAMES.write[ grep("MNa_",COV_NAMES.write) ] <- COV_NAMES.write[ grep("PRE_MNa_",COV_NAMES.write) ]
COV_NAMES.write[ grep("Bdr_",COV_NAMES.write) ] <- COV_NAMES.write[ grep("PRE_Bdr_",COV_NAMES.write) ]
COV_NAMES.write[ grep("DEL",COV_NAMES.write) ] <- ""
 # Combine Lists into Single Table
PHENO_COV_LIST <- data.frame( PHENO_NAMES.write, COV_NAMES.write )

## Write Pheno/Cov Table
PathToWrite <- gsub("20150310_Single_Pheno","20150313_PhenoCov",PathToWAG)
write.table( PHENO_COV_LIST, PathToWrite, sep="\t",row.names=F,col.names=F,quote=F )

# ## Test Correlation b/n Initial Values & Delta-Values
# for ( r in 1:nrow(PHENO_COV_LIST) ) {
# 	resp <- as.character( PHENO_COV_LIST[r,1] )
# 	pred <- as.character( PHENO_COV_LIST[r,2] )
# 	if (pred!="") {
# 		TEMP_TAB <- FULL_TABLE[ , c(resp,pred) ]
# 		colnames(TEMP_TAB) <- c("resp","pred")
# 		MOD <- lm( resp ~ pred, data=TEMP_TAB )
# 		print(paste( "######",resp,"~",pred ))
# 		print( summary(MOD)$coefficients["pred","Pr(>|t|)"] )	
# 	}
# }







###########################################################
## END OF DOC #############################################
###########################################################



# for ( col in 1:ncol(D.MN.1.diff) ) {
# 	par(ask=T)
# 	pairs( data.frame( D.MN.1.diff[,col],D.MN.2.diff[,col],D.MN.3.diff[,col],D.MN.4.stat[,col],D.MN.5.B[,col]) )
# }


# pairs( D.MN.4.stat )





















