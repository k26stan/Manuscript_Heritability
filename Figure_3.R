## Quantify Placebo Effect using Time-Series Data ##
## July 8, 2015 ##
## Kristopher Standish ##

library( nlme )

###########################################################
## LOAD DATA ##############################################
###########################################################

## Set Date
DATE <- gsub("-","",Sys.Date())

## Set Paths to Data and to Save
# PathToFT <- "/Users/kstandis/Data/Burn/Data/Phenos/Full_Tables/20141229_Full_Table.txt"
PathToFT <- "/Users/kstandis/Data/Janssen/Data/Pheno/Derived/20151015_Full_Table.txt"
# PathToRep <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150226_Resp_v_Time.txt"
# PathToRep <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150530_Resp_v_Time.txt"
PathToRep <- "/Users/kstandis/Data/Janssen/Data/Pheno/Derived/20160112_Resp_v_Time.txt"
PathToSave <- paste("/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Plots/",DATE,"_Plac/",sep="")
dir.create( PathToSave )

## Load Real Data
FT <- read.table( PathToFT, sep="\t",header=T )
RP <- read.table( PathToRep, sep="\t",header=T )

###########################################################
## ORGANIZE DATA ##########################################
###########################################################

## Calculate Transformations of Phenotypes
 # DAS/lCRP/rSJC/rTJC/rSJC28/rTJC28
PH <- data.frame( lCRP=log10(RP$CRP), rSJC=sqrt(RP$SJC), rTJC=sqrt(RP$TJC), rSJC28=sqrt(RP$SJC28), rTJC28=sqrt(RP$TJC28) )
DAT.1 <- data.frame( RP, PH )
WHICH_PHENOS <- c("DAS","lCRP","rSJC","rTJC")
PH_COLS <- which( colnames(DAT.1) %in% WHICH_PHENOS )# c(16,21:25)
names(PH_COLS) <- colnames(DAT.1)[PH_COLS]

## Remove timepoints w/ DRUG==1
DAT.2 <- DAT.1[ which(DAT.1$DRUG==0), ]

###########################################################
## MODEL DRUG and PLACEBO EFFECTS #########################
###########################################################

## Gameplan
 # Show Multiple Linear Regression Coefficients
   # Barplot of Betas for Intercept, Placebo, Drug
   # Print P-Values & **
 # Show Placebo as proportion of Drug Effect & Intercept

###########################################################
## Simple Linear Regression w/ PLAC+DRUG as fixed effect
MODS <- list()
for ( c in 1:length(PH_COLS) ) {
	col <- PH_COLS[c]
	pheno <- names(PH_COLS)[c]
	FORMULA <- as.formula( paste( pheno, "~ DRUG+PLAC" ) )
	MOD <- lm( FORMULA, data=DAT.1 )
	# Compile
	MODS[[pheno]] <- MOD
}

###########################################################
## COMPILE MODEL RESULTS ##################################
###########################################################

## Pull out Coefficients, Standard Errors, and P-Values
COEFS <- matrix( unlist(lapply( MODS, function(x) summary(x)$coefficients[,"Estimate"] )), byrow=T,ncol=3 )
STERS <- matrix( unlist(lapply( MODS, function(x) summary(x)$coefficients[,"Std. Error"] )), byrow=T,ncol=3 )
PS <- matrix( unlist(lapply( MODS, function(x) summary(x)$coefficients[,"Pr(>|t|)"] )), byrow=T,ncol=3 )
# colnames(COEFS) <- colnames(STERS) <- colnames(PS) <- c("INT","DRUG","PLAC")
colnames(COEFS) <- colnames(STERS) <- colnames(PS) <- c("BL","GOL","PBO")
rownames(COEFS) <- rownames(STERS) <- rownames(PS) <- names(MODS)

## Get Ratio of Placebo Effect to Intercept & Drug EFfect
RAT.pi <- abs( COEFS[,"PBO"] / COEFS[,"BL"] )
RAT.pd <- abs( COEFS[,"PBO"] / COEFS[,"GOL"] )
RAT <- rbind( RAT.pi, RAT.pd )

## Significance
Stars <- PS
Stars[ which(PS<.05) ] <- "*"
Stars[ which(PS<.01) ] <- "**"
Stars[ which(PS<.001) ] <- "***"

###########################################################
## PLOT FIGURE 3 ##########################################
###########################################################

## Open File to Write Figure to
png( paste(PathToSave,"3_Full.png",sep="/"), height=1000,width=2000, pointsize=32 )

## Set Plotting Parameters
COLS <- c("chocolate2","slateblue3","turquoise3")
par(mfrow=c(1,2)) # layout( matrix(c(1,2,3,4), 2, 2, byrow = TRUE), widths=c(3,2), heights=c(1,1) )
par(mar=c(3,4,4,1))
## Coefs for Intercept, Drug, and Placebo
YLIM <- extendrange( COEFS, f=.15 )
X_BARS <- barplot( t(COEFS), col=COLS, beside=T, ylim=YLIM, ylab="Population Estimate",main="Population Coefficient Estimates" )
abline( h=seq(floor(YLIM[1]),YLIM[2]+1,1), lty=3,col="grey50",lwd=1 )
barplot( t(COEFS), col=COLS, beside=T, add=T )
 # Error Bars
X.arrow <- c(t( X_BARS )) # rep( c(1.5,2.5), rep(6,2) ) + rep( seq(0,15,3),2 )
arrows( X.arrow, c(COEFS)+c(STERS), X.arrow, c(COEFS)-c(STERS), code=3,angle=90,length=0.2,lwd=2 )
 # Text
text( X.arrow[9:12], (c(COEFS)-c(STERS)-.2)[9:12], label=c(Stars)[9:12] )
text( X.arrow[9:12], 1.2, label=paste("p=",formatC(PS[9:12],format="e",digits=2),sep=""), srt=90 )
 # Legend
legend( "topright", fill=COLS, legend=colnames(COEFS), bg="white",ncol=3 )

## Placebo Relative to Intercept & Drug
YLIM <- c( 0, max(RAT) ) * c(1,1.2)
X_BARS <- barplot( RAT, col=COLS[1:2], beside=T, ylim=YLIM, ylab="Placebo Effect Size Ratio",main="Proportional Placebo Effect Size" )
abline( h=seq(0,1,.1), lty=3,col="grey50",lwd=1 )
barplot( RAT, col=COLS[1:2], beside=T, add=T )
 # Put actual % Values on Plot
X.arrow <- c( X_BARS ) 
text( X.arrow, c(RAT)+.04, label=round(c(RAT),3), srt=90 )
 # Legend
legend( "topleft", fill=COLS[1:2], legend=paste("PBO/",colnames(COEFS)[1:2],sep=""),ncol=2 )

dev.off()


###########################################################
## END OF DOC #############################################
###########################################################
