## Simulate Association Data varying Linear Regression Assumptions ##
## February 25, 2015 ##
## Kristopher Standish ##

###########################################################
## Get Set Up #############################################
###########################################################

## Set Date
DATE <- "20150226"

## Set Paths to Data and to Save
PathToData <- "/Users/kstandis/Data/Burn/20141229_Full_Table.txt"
PathToSave <- paste("/Users/kstandis/Dropbox/Schork/JNJ11/Writing/Resp_Herit/Plots/",DATE,sep="")

## Load Real Data
TAB <- read.table( PathToData, sep="\t",header=T )

## Game Plan
 # Pick # Patients
 # Simulate Phenotype Values
   # (Under Different Assumptions)
     # Start w/ DAS simulation ( rnorm(5,1,N.pats) )
 # 

## Number of Patients
N.pats <- 200

## Number of Simulations
N.sims <- 5000

###########################################################
## FCT: Simulate Data #####################################
###########################################################

## Create Function to Simulate/Test Data
SIM.DAT <- function( N.pats, N.sims, MULT ) {

	## Create Compile Variables
	P.vals <- T.vals <- array( , c(N.sims,8) )
	colnames(P.vals) <- colnames(T.vals) <- c("TTT.1","TTT.2","TTT.12.1","TTT.12.2","TTF.1","TTF.2","TTF.12.1","TTF.12.2")

	## Loop Through Simulations
	RAND_UNIF <- runif( N.sims, 0, 1 )
	start_time <- proc.time()
	for ( i in 1:N.sims ) {

		##############################################
		## SIMULATE ##
		
		## Simulate Var.1 Values
		X1 <- rnorm( N.pats, 0, 1 )

		## Simulate Var.2 Values
		X2 <- sample( c(0,1), N.pats, replace=T )

		## Simulate Phenotype Values
		 # Null, Normal, Homoscedastic == TTT
		Y.TTT <- rnorm( N.pats, 0, 1+MULT*mean(X1-min(X1)) )
		Y.TTF <- rnorm( N.pats, 0, 1+MULT*(X1-min(X1)) )
		# Y.TFT <- 
		# Y.TFF <- 
		# Y.FTT <- BETA*X1 + rnorm( N.pats, 0, 1+mean(X1-min(X1)) )
		# Y.FTF <- BETA*X1 + rnorm( N.pats, 0, 1+(X1-min(X1)) )
		# Y.FTT <- rnorm( N.pats, X1, 1+mean(X1-min(X1)) )
		# Y.FTF <- rnorm( N.pats, X1, 1+(X1-min(X1)) )
		# Y.FFT <- 
		# Y.FFF <- 

		##############################################
		## MODEL ##

		## Model Y Vars vs X Vars
		 # TTT
		LM.1.TTT <- lm( Y.TTT ~ X1 )
		LM.2.TTT <- lm( Y.TTT ~ X2 )
		LM.12.TTT <- lm( Y.TTT ~ X1 + X2 )
		 # TTF
		LM.1.TTF <- lm( Y.TTF ~ X1 )
		LM.2.TTF <- lm( Y.TTF ~ X2 )
		LM.12.TTF <- lm( Y.TTF ~ X1 + X2 )

		## Pull P-Values
		 # TTT
		P.1.TTT <- summary( LM.1.TTT )$coefficients["X1","Pr(>|t|)"]
		P.2.TTT <- summary( LM.2.TTT )$coefficients["X2","Pr(>|t|)"]
		P.12.TTT.1 <- summary( LM.12.TTT )$coefficients["X1","Pr(>|t|)"]
		P.12.TTT.2 <- summary( LM.12.TTT )$coefficients["X2","Pr(>|t|)"]
		 # TTF
		P.1.TTF <- summary( LM.1.TTF )$coefficients["X1","Pr(>|t|)"]
		P.2.TTF <- summary( LM.2.TTF )$coefficients["X2","Pr(>|t|)"]
		P.12.TTF.1 <- summary( LM.12.TTF )$coefficients["X1","Pr(>|t|)"]
		P.12.TTF.2 <- summary( LM.12.TTF )$coefficients["X2","Pr(>|t|)"]

		## Pull t-Statistics
		 # TTT
		T.1.TTT <- summary( LM.1.TTT )$coefficients["X1","t value"]
		T.2.TTT <- summary( LM.2.TTT )$coefficients["X2","t value"]
		T.12.TTT.1 <- summary( LM.12.TTT )$coefficients["X1","t value"]
		T.12.TTT.2 <- summary( LM.12.TTT )$coefficients["X2","t value"]
		 # TTF
		T.1.TTF <- summary( LM.1.TTF )$coefficients["X1","t value"]
		T.2.TTF <- summary( LM.2.TTF )$coefficients["X2","t value"]
		T.12.TTF.1 <- summary( LM.12.TTF )$coefficients["X1","t value"]
		T.12.TTF.2 <- summary( LM.12.TTF )$coefficients["X2","t value"]

		## Compile Values
		P.vals[i,] <- c( P.1.TTT, P.2.TTT, P.12.TTT.1, P.12.TTT.2, P.1.TTF, P.2.TTF, P.12.TTF.1, P.12.TTF.2 )
		T.vals[i,] <- c( T.1.TTT, T.2.TTT, T.12.TTT.1, T.12.TTT.2, T.1.TTF, T.2.TTF, T.12.TTF.1, T.12.TTF.2 )

		## Plot first one and then randomly
		if ( i==1 | RAND_UNIF[i]<.001 ) {
			COLS <- c("deepskyblue2","firebrick2")
			par(mfrow=c(1,2))
			plot( Y.TTT ~ X1, col=COLS[factor(X2)], pch="+" )
			abline( LM.1.TTT, col="black" )
			plot( Y.TTF ~ X1, col=COLS[factor(X2)], pch="+" )
			abline( LM.1.TTF, col="black" )
		}

		## Update
		print_time <- round( proc.time()-start_time, 3)[3]
		if ( i%%100==0 ) { print(paste( "Done with",i,"of",N.sims,"-",print_time )) }

	} # Close "i" Loop
	
	## Return Tables
	COMPILE <- list( P.vals, T.vals )
	names(COMPILE) <- c("P","T")
	return(COMPILE)

} # Close "SIM.DAT" function
# OUT <- SIM.DAT( N.pats, N.sims, 2 )

###########################################################
## Simulate w/ Different Variant Ranges ###################
###########################################################

## Specify Number of Patients & Simulations
N.pats <- 500
N.sims <- 10000

## Specify Range for Multiplier values
MULT_RANGE <- c( 0.1, 0.2, 0.5, 1, 2, 5, 10 )

## Create Compile Variables
M <- list()
PROPS <- array( , c(length(MULT_RANGE),8) )

## Loop through different Multipler values
m_time <- proc.time()
for ( m in 1:length(MULT_RANGE) ) {
	MULT <- MULT_RANGE[m]
	print(paste("#### Running Multiplier =",MULT ))
	print(paste("#### Loop",m,"of",length(MULT_RANGE),"-",round(proc.time()-m_time,3)[3] ))
	## Run Simulations using this Multiplier
	M[[m]] <- SIM.DAT( N.pats, N.sims, MULT )
	## Pull out Significance
	P.vals <- M[[m]]$P
	## Show what fraction of P-vals are <.05
	P.vals.sig <- P.vals <= .05
	P.vals.sig.count <- apply( P.vals.sig, 2, table )
	P.vals.sig.prop <- prop.table( P.vals.sig.count, 2 )
	PROPS[m,] <- P.vals.sig.prop[2,]
	 # Plot it
	COLS <- c("deepskyblue2","firebrick2")
	barplot( P.vals.sig.prop, col=COLS, las=2 )
	abline( h=0.95, lty=2,col="black",lwd=2 )

}
colnames(PROPS) <- colnames(P.vals.sig.prop)
rownames(PROPS) <- paste("M",MULT_RANGE,sep="_")

## Save Data
write.table( PROPS, paste(PathToSave,"PROPS.txt",sep="/"), sep="\t",row.names=T,col.names=T,quote=F )
save( M, file=paste(PathToSave,"M.Rdata",sep="/") )

## Plot it
XLIM <- range( log10(MULT_RANGE) )
YLIM <- range( PROPS ) # c( 0, max(PROPS) )
COLS.list <- c( "firebrick3","chocolate2","gold2","springgreen2","steelblue2","slateblue3" )
COLS <- colorRampPalette(COLS.list)(8)
LTYS <- rep( c(1,2), 4 ) # c(1,1,2,2,1,1,2,2)
png( paste(PathToSave,"1-FDR.png",sep="/"), height=800,width=1200, pointsize=26 )
plot( 0,0,type="n", xlim=XLIM, ylim=YLIM, xlab="-log10(MULT)", ylab="FDR", main="FDR vs Multiplier" )
abline( h=seq(0,1,.01), lty=2,col="grey50",lwd=1)
abline( h=.05, lty=1,col="black",lwd=1)
for ( col in 1:8 ) {
	points( log10(MULT_RANGE), PROPS[,col], col=COLS[col], type="o",pch=20,lty=LTYS[col],lwd=2 )
}
legend( "topleft", col=COLS, legend=colnames(PROPS), pch=20,lty=LTYS,lwd=2, ncol=2 )
dev.off()

###########################################################
## Plot Results ###########################################
###########################################################

P.vals <- OUT$P
T.vals <- OUT$T

## Show what fraction of P-vals are <.05
P.vals.sig <- P.vals <= .05
P.vals.sig.count <- apply( P.vals.sig, 2, table )
P.vals.sig.prop <- prop.table( P.vals.sig.count, 2 )
 # Plot it
COLS <- c("deepskyblue2","firebrick2")
barplot( P.vals.sig.prop, col=COLS, las=2 )
abline( h=0.95, lty=2,col="black",lwd=2 )




## Power ##
N <- 5000
P_VALS.0.h <- P_VALS.0.v <- P_VALS.1.h <- P_VALS.1.v <- numeric(N)
T_VALS.0.h <- T_VALS.0.v <- T_VALS.1.h <- T_VALS.1.v <- numeric(N)

start_time <- proc.time()
for ( i in 1:N ) {
	## Simulate X & Y Values
	XX <- rnorm( 100,0,2 )
	YY.0.h <- rnorm( 100,0,mean(XX-min(XX)) ) # Null - Homoscedasctic
	YY.0.v <- rnorm( 100,0,(XX-min(XX)) ) # Null - Violates Homoscedasticity
	BETA <- .3
	YY.1.h <- BETA*XX+rnorm( 100,0,mean(XX-min(XX)) ) # Assoc - Homoscedasctic
	YY.1.v <- BETA*XX+rnorm( 100,0,(XX-min(XX)) ) # Assoc - Violates Homoscedasticity
	## Model Linear Fit
	MOD.0.h <- lm( YY.0.h ~ XX )
	MOD.0.v <- lm( YY.0.v ~ XX )
	MOD.1.h <- lm( YY.1.h ~ XX )
	MOD.1.v <- lm( YY.1.v ~ XX )
	## Get P & T Values
	TAB.0.h <- summary(MOD.0.h)$coefficients
	T.0.h <- TAB.0.h["XX","t value"]
	P.0.h <- TAB.0.h["XX","Pr(>|t|)"]
	TAB.0.v <- summary(MOD.0.v)$coefficients
	T.0.v <- TAB.0.v["XX","t value"]
	P.0.v <- TAB.0.v["XX","Pr(>|t|)"]
	TAB.1.h <- summary(MOD.1.h)$coefficients
	T.1.h <- TAB.1.h["XX","t value"]
	P.1.h <- TAB.1.h["XX","Pr(>|t|)"]
	TAB.1.v <- summary(MOD.1.v)$coefficients
	T.1.v <- TAB.1.v["XX","t value"]
	P.1.v <- TAB.1.v["XX","Pr(>|t|)"]
	## Compile Values
	P_VALS.0.h[i] <- P.0.h
	T_VALS.0.h[i] <- T.0.h
	P_VALS.0.v[i] <- P.0.v
	T_VALS.0.v[i] <- T.0.v
	P_VALS.1.h[i] <- P.1.h
	T_VALS.1.h[i] <- T.1.h
	P_VALS.1.v[i] <- P.1.v
	T_VALS.1.v[i] <- T.1.v
	# T.txt <- formatC( T, format="e",digits=3 )
	# P.txt <- formatC( P, format="e",digits=3 )
	# XLIM <- range(XX)
	# YLIM <- range(YY)
	# plot( YY ~ XX )
	# abline( MOD )
	# text( quantile(XLIM,.1),quantile(YLIM,.9), paste("P=",P.txt,"; T=",T.txt,sep=""), pos=4 )	
	if ( i%%100==0 ) { print(paste(i,"of",N,"-",round(proc.time()-start_time,3)[3] )) }
}

length(which(P_VALS.0.h<.05)) / N
length(which(P_VALS.0.v<.05)) / N
length(which(P_VALS.1.h<.05)) / N
length(which(P_VALS.1.v<.05)) / N

COLS <- c("firebrick2","gold2","chartreuse2","blue3")
BRKS <- seq( -10,10,.25 )
par(mfrow=c(1,2))
hist( T_VALS.0.h, col=COLS[1], density=20,angle=45, breaks=BRKS )
hist( T_VALS.0.v, col=COLS[2], density=20,angle=-45, breaks=BRKS, add=T )
hist( T_VALS.1.h, col=COLS[3], density=20,angle=45, breaks=BRKS )
hist( T_VALS.1.v, col=COLS[4], density=20,angle=-45, breaks=BRKS, add=T )

plot( YY.0.h ~ XX, col=COLS[1], pch="+" )
points( YY.0.v ~ XX, col=COLS[2], pch="+" )
points( YY.1.h ~ XX, col=COLS[3], pch="+" )
points( YY.1.v ~ XX, col=COLS[4], pch="+" )
abline( MOD.0.h, col=COLS[1] )
abline( MOD.0.v, col=COLS[2] )
abline( MOD.1.h, col=COLS[3] )
abline( MOD.1.v, col=COLS[4] )






## Simulate Initial Phenotype Values
 # DAS
DAS.0.mn <- 6
DAS.0.sd <- 1
DAS.0.sim <- rnorm( N.pats, DAS.0.mn, DAS.0.sd )

