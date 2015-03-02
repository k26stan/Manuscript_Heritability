## Quantify Placebo Effect using Time-Series Data ##
## February 26, 2015 ##
## Kristopher Standish ##

###########################################################
## LOAD DATA ##############################################
###########################################################

## Set Date
DATE <- "20150226"

## Set Paths to Data and to Save
PathToFT <- "/Users/kstandis/Data/Burn/20141229_Full_Table.txt"
PathToRep <- "/Users/kstandis/Data/Burn/20140731_Resp_v_Time.txt"
PathToSave <- paste("/Users/kstandis/Dropbox/Schork/JNJ11/Writing/Resp_Herit/Plots/",DATE,sep="")

## Load Real Data
FT <- read.table( PathToFT, sep="\t",header=T )
RP <- read.table( PathToRep, sep="\t",header=T )

## Game Plan







###########################################################
## END OF DOC #############################################
###########################################################
