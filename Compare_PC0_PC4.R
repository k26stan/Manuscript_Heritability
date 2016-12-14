## Compare Results of 0 and 4 PC models ##
## June 23, 2016 ##
## Kristopher Standish ##

library( gplots )

###########################################################
## LOAD DATA ##############################################
###########################################################

## Set Date
DATE <- gsub("-","",Sys.Date())

 # 4 PC Model
PathToSing <- "/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Data/GCTA_Results/PC4_Sing/4-PERM_Compile.Rdata"
PathToDer <- "/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Data/GCTA_Results/PC4_Der/4-PERM_Compile.Rdata"
PathToResTabs <- "/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Plots/20160622_GCTA_PC4/TAB"
load(PathToSing) ; SING.4 <- COMPILE.full
load(PathToDer) ; DER.4 <- COMPILE.full

 # No PC Model
PathToSing <- "/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Data/GCTA_Results/PC0_Sing/4-PERM_Compile.Rdata"
PathToDer <- "/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Data/GCTA_Results/PC0_Der/4-PERM_Compile.Rdata"
PathToResTabs <- "/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Plots/20160623_GCTA_PC0/TAB"
load(PathToSing) ; SING.0 <- COMPILE.full
load(PathToDer) ; DER.0 <- COMPILE.full

 # PathToSave
PathToSave <- paste("/Users/kstandis/Dropbox/Schork/JNJ11/Manuscripts/Resp_Herit/Plots/",DATE,"_GCTA_PC0vPC4",sep="")
dir.create( PathToSave )



# load( "20160622_GCTA_MAF1_SNP_PC0_20160510_Derived_Phenos/4-PERM_Compile.Rdata" )
# DER.0 <- COMPILE.full
# load( "20160621_GCTA_MAF1_SNP_PC4_20160510_Derived_Phenos/4-PERM_Compile.Rdata" )
# DER.4 <- COMPILE.full
# # load( "20160621_GCTA_MAF1_SNP_PC10_20160510_Derived_Phenos/4-PERM_Compile.Rdata" )
# # DER.10 <- COMPILE.full

## 0 PC Model
 # Derived
VAR.d0 <- Reduce( rbind, lapply( DER.0$VAR, tail, 1 ))
SE.d0 <- Reduce( rbind, lapply( DER.0$SE, tail, 1 ))
rownames(VAR.d0) <- rownames(SE.d0) <- names(DER.0$VAR)
 # Single
VAR.s0 <- Reduce( rbind, lapply( SING.0$VAR, tail, 1 ))
SE.s0 <- Reduce( rbind, lapply( SING.0$SE, tail, 1 ))
rownames(VAR.s0) <- rownames(SE.s0) <- names(SING.0$VAR)
## 4 PC Model
 # Derived
VAR.d4 <- Reduce( rbind, lapply( DER.4$VAR, tail, 1 ))
SE.d4 <- Reduce( rbind, lapply( DER.4$SE, tail, 1 ))
rownames(VAR.d4) <- rownames(SE.d4) <- names(DER.4$VAR)
 # Single
VAR.s4 <- Reduce( rbind, lapply( SING.4$VAR, tail, 1 ))
SE.s4 <- Reduce( rbind, lapply( SING.4$SE, tail, 1 ))
rownames(VAR.s4) <- rownames(SE.s4) <- names(SING.4$VAR)

## Plot Model Comparison
 # Specify Colors
COLS <- c("firebrick1","gold1","chartreuse1","dodgerblue1","gold4","chartreuse4","dodgerblue4") # c("firebrick2","gold2","chartreuse2","deepskyblue2","slateblue2")
names(COLS) <- c("DAS","rSJC","lCRP","rTJC","SJC","CRP","TJC")
PHE.s <- sapply(strsplit(rownames(VAR.s4),"_"),"[",3)
PHE.d <- sapply(strsplit(rownames(VAR.d4),"_"),"[",3)
COLS.s <- COLS[PHE.s]
COLS.d <- COLS[PHE.d]
 # pch
PCHS <- 0:4
NAMES.s <- c("WAG4","WAG12","WAG20","WAG28","FL")
NAMES.d <- c("MNa","MNcd","PRC","Bwk","VARwk")
PHE.s <- sapply(strsplit(rownames(VAR.s4),"_"),"[",2)
PHE.d <- sapply(strsplit(rownames(VAR.d4),"_"),"[",2)
names(PCHS) <- NAMES.s
PCHS.s <- PCHS[PHE.s]
names(PCHS) <- NAMES.d
PCHS.d <- PCHS[PHE.d]

 # Correlation b/n Models
COR.d <- cor( VAR.d4[,"VgVp"], VAR.d0[,"VgVp"] )
COR.s <- cor( VAR.s4[,"VgVp"], VAR.s0[,"VgVp"] )
png( paste(PathToSave,"/20160621_Var_PC4vPC10.png",sep=""), height=1200, width=2400, pointsize=30 )
par(mfrow=c(1,2))
 # Single Phenos
plot( 0,0,type="n", xlim=c(0,1),ylim=c(0,1),
	main="Vg/Vp Estimates in PC0 & PC4 Models",xlab="PC0 Model",ylab="PC4 Model" )
abline( h=seq(0,1,.1), v=seq(0,1,.1), lty=3,col="grey50" )
abline( 0,1 )
points( VAR.s0[,"VgVp"], VAR.s4[,"VgVp"], pch=PCHS.s,col=adjustcolor(COLS.s,.7),cex=1.2,lwd=4 )
text( 1,.25, paste("R(s) =",round(COR.s,3)), pos=2, col="black" )
legend(.70,.2, pch=PCHS,col=adjustcolor("black",.7),pt.cex=1.2,legend=NAMES.s,title="Time Point",bg="white", cex=.8,ncol=2,pt.lwd=4 )
legend(.45,.2, pch=16,col=adjustcolor(COLS,.7),pt.cex=1.2,legend=names(COLS),title="Phenotypes",bg="white", cex=.8,ncol=2 )
 # Derived Phenos
plot( 0,0,type="n", xlim=c(0,1),ylim=c(0,1),
	main="Vg/Vp Estimates in PC0 & PC4 Models",xlab="PC0 Model",ylab="PC4 Model" )
abline( h=seq(0,1,.1), v=seq(0,1,.1), lty=3,col="grey50" )
abline( 0,1 )
points( VAR.d0[,"VgVp"], VAR.d4[,"VgVp"], pch=PCHS.d,col=adjustcolor(COLS.d,.7),cex=1.2,lwd=4 )
text( 1,.25, paste("R(d) =",round(COR.d,3)), pos=2, col="black" )
legend(.75,.2, pch=PCHS,col=adjustcolor("black",.7),pt.cex=1.2,legend=NAMES.d,title="Derived Stat",bg="white", cex=.8,ncol=2,pt.lwd=4 )
# legend(.5,.2, pch=16,col=adjustcolor(COLS,.7),pt.cex=1.2,legend=names(COLS),title="Phenotypes",bg="white", cex=.8,ncol=2 )
dev.off()


hist(VAR.s0[,4],col=adjustcolor(COLS[1],.5),breaks=seq(-.1,1.1,.05))
hist(VAR.s4[,4],col=adjustcolor(COLS[2],.5),breaks=seq(-.1,1.1,.05),add=T)






