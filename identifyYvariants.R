####packages####
library(data.table) #fread()
library(rgdal) #readOGR()
library(dismo) #getData() library(maptools)

library(ggplot2)
library(ggrepel)
library(ggExtra)
#source("~/Google Drive/Research/Scripts/multiplot_fnc.R") #multiplot

library(plotrix) #FUN=std.error
library(sqldf)

library(FactoMineR) #PCA
library(HIest)

library(elevatr) #get elevation

library(plyr)
library(dplyr)
library(tidyr) #separate

library(hzar)
if(require(doMC)){  ## If you have doMC, use foreach in parallel mode to speed up computation.
  registerDoMC()
} else {## Use foreach in sequential mode
  registerDoSEQ();
}
library(geosphere) #distm

library(corrplot)
source("~/Google Drive/Research/Scripts/baypass_2.2/utils/baypass_utils.R")
library(reshape2) #melt

#EEMS#

library(Rcpp)
library(RcppEigen)
library(raster)
library(rgeos)
library(sp)
library(scatterpie)

#setwd('~/Google Drive/Research/Scripts/eems-master/plotting/')
#install.packages("rEEMSplots", repos = NULL, type = "source")
library(rEEMSplots)

#library(rgdal)
library(rworldmap)
library(rworldxtra)

library(cowplot)

library(qualV)
require(splines)
source("~/Google Drive/Research/Scripts/linked-selection-master/linkedsel_functions.R")
est_rec_piece_func  <- function(x) {est_rec_piece(x)}

library(conStruct) #version(conStruct)
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

library(directlabels)

na.test <-  function (x) {
  w <- sapply(x, function(x)all(is.na(x)))
  if (any(w)) {
    stop(paste("All NA in columns", paste(which(w), collapse=", ")))
  }
}

library("viridis") 

library(sf)
library(corrplot)

library(ape)

sessionInfo()
'%ni%' <- Negate('%in%')
####file intake####
setwd('~/Google Drive/Research/Data2/RumexhybridGenome/')

maf="01"

#add altitude
#pop_pos <- fread('pickupPops_hyb_round.txt')
#location <-  cbind.data.frame("decimalLongitude"=pop_pos$Longitude,"decimalLatitude"=pop_pos$Latitude)
#projection_none <- "+proj=longlat +datum=WGS84"
#space <- SpatialPoints(location, proj4string = CRS(projection_none))
#alt <- get_elev_point(space)
#pop_pos <- cbind.data.frame(pop_pos,"elevation"=alt$elevation)
#write.table(pop_pos, file = "pickupPops_hyb_round_elev.txt", append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"))

pop_pos <- fread('pickupPops_hyb_round_elev.txt')
#length(pop_pos$Pop[pop_pos$Collec == "New"])
#length(pop_pos$Pop[pop_pos$Collec == "Old"])

#length(pop_pos_gen$Pop[pop_pos_gen$Collec == "New"])

pop_pos_n <- pop_pos[pop_pos$Collec == "New",c(4,8)]
pop_pos_o <- pop_pos[pop_pos$Collec == "Old",c(8,8)]
names(pop_pos_n) <- c("mixName","Name")
names(pop_pos_o) <- c("mixName","Name")

mixname <- rbind.data.frame(pop_pos_n,pop_pos_o)
pop_pos <- left_join(pop_pos,mixname,by=c("Name"="Name"))

#inds
indInfo <- fread('HiC_hyb.indv.txt')
indInfo$flower[indInfo$flower == "U"] <- "Inc"
indPopInfo <- sqldf('select indInfo.*, pop_pos.*   from indInfo left join pop_pos on indInfo.statePop = pop_pos.Name')
indPopInfo <- indPopInfo[,c(1:8,11:13,15,16)]
inds <- fread(paste("GBS.mis60.maf",maf,".012.indv",sep=""),header=FALSE)
#length(indPopInfo$name[indPopInfo$Collec == "New"])
#length(indPopInfo$name[indPopInfo$Collec == "Old"])
#length(indPopInfo$name)


indsInfo <- sqldf('select inds.*, indPopInfo.* from inds left join indPopInfo on inds.V1 = indPopInfo.name')
indsInfo$Type = factor(indsInfo$Type, levels=c('XY','Hybrid','XYY'))
popSort <- as.data.frame(unique(indsInfo$statePop))
names(popSort) <- c("pop")
pop_pos_gen <- pop_pos[pop_pos$Name %in% popSort$pop]

pop_tally <- indsInfo %>% group_by(statePop) %>%  tally()
pop_pos_gen <- left_join(pop_pos_gen,pop_tally,by=c("Name"="statePop"))

pheno_pos <- fread('populations_Rumex_hastatulus.csv')

s012 <- fread(paste("GBS.mis60.maf",maf,".012",sep=""),na.strings = "-1")

s012 <- s012[,-1]

#pos

loc <- fread(paste("GBS.mis60.maf",maf,".012.chrom.pos",sep=""))
names(loc) <- c("scaffold","sPos","LG","lPos")

loc$LG[loc$LG == "L.10"] <- "X"
loc$LG[loc$LG == "L.7"] <- "A1"
loc$LG[loc$LG == "L.8"] <- "A2"
loc$LG[loc$LG == "L.5"] <- "A3"
loc$LG[loc$LG == "L.3"] <- "A4"
loc$LG[is.na(loc$LG)] <- "U"

maxA3pos <- 175014771
loc$lPos[loc$LG == "A3"] <- maxA3pos - loc$lPos[loc$LG == "A3"] 

loc$mb <- loc$lPos / 1000000
loc$SNP <- seq(1,length(loc$lPos),1)
loc$chrom = factor(loc$LG, levels=c('A1','A2','A4','X','A3','U'))

#list sites where all individuals are heterozygotes
colvar0<-apply(s012,2,function(x) var(x,na.rm=T)==0)
removed_sites <- as.numeric(gsub("V", "", names(s012)[colvar0])) - 1
#s012HI <- s012[,-..fixedHet] #remove these sites; specific to these data
#locHI <- loc[-fixedHet,]
#s012HIm <- as.matrix(s012HI)
#locHI$SNP <- seq(1,length(locHI$lPos),1)

indsInfo$hybrid_cat[indsInfo$Collec == "Old" & indsInfo$Type == "XYY"] <- "XYY"
indsInfo$hybrid_cat[indsInfo$Collec == "Old" & indsInfo$Type == "XY"] <- "XY"
indsInfo$hybrid_cat[indsInfo$Collec == "New" & indsInfo$riverside == "E"] <- "Hybrid_E"
indsInfo$hybrid_cat[indsInfo$Collec == "New" & indsInfo$riverside == "W"] <- "Hybrid_W"

for(hybrid_cat in c("XYY","XY","Hybrid_E","Hybrid_W")){
  # for(chrom in c("X","A3")){
  
  # chrom_cols  <- loc$SNP[loc$LG %in% chrom ] 
  
  # s012_chrom <-  as.data.frame( s012[,..chrom_cols])
  #subset XX/XYY inds
  s012_chrom_pop <-  s012[indsInfo$hybrid_cat == hybrid_cat,]
  
  #ID and remove sites missing in all XX/XYY
  s012_chrom_pop_NAs <- sapply(s012_chrom_pop, function(x)all(is.na(x)))
  
  removed_sites <- unique(c(removed_sites,as.numeric(gsub("V", "", names(s012_chrom_pop_NAs[s012_chrom_pop_NAs == T]))) - 1))
  
  #XYY_NAcols <- as.numeric(gsub("V", "", names( s012_XYY_NAs[s012_XYY_NAs == T]))) - 1
  # cat(hybrid_cat," ",chrom," ",removed_sites,"\n")
  
  #}
}

#mostly_miss <- as.numeric(gsub("V", "", hybrid_cat)) - 1
s012HI <- s012[,-..removed_sites] #remove these sites; specific to these data
locHI <- loc[loc$SNP %ni% removed_sites,]


####
#Case 1: #THIS is actually a fixed difference on the neoX or X, which we would want to include

#all females in TX and males in TX 0/0

#list of columns where everyone is same
allTX <- vapply(s012[indsInfo$hybrid_cat == "XY",], function(x) length(na.omit(unique(x))) == 1, logical(1L))

#subset by list of columns where everyone is same
s012_allTX <- s012[indsInfo$hybrid_cat == "XY",..allTX]

#calculate product: 0 is 0, 1 is 1 and >1 is 2
s012_allTX_prod <- apply(s012_allTX, 2, prod,na.rm=TRUE)
s012_allTX_prod[s012_allTX_prod > 1] <- 2
s012_allTX_prod_names <- cbind.data.frame("SNP"=names(s012_allTX_prod),"TX"=s012_allTX_prod)


#All females in NC 1/1
allNC_F <- vapply(s012[indsInfo$hybrid_cat == "XYY" & indsInfo$flower == "F",], function(x) length(na.omit(unique(x))) == 1, logical(1L))
s012_allNC_F <- s012[indsInfo$hybrid_cat == "XY",..allNC_F]
s012_allNC_F_prod <- apply(s012_allNC_F, 2, prod,na.rm=TRUE)
s012_allNC_F_prod[s012_allNC_F_prod > 1] <- 2
s012_allNC_F_names <- cbind.data.frame("SNP"=names(s012_allNC_F_prod),"NC_F"=s012_allNC_F_prod)

#all males 1/0
allNC_M <- vapply(s012[indsInfo$hybrid_cat == "XYY" & indsInfo$flower == "M",], function(x) length(na.omit(unique(x))) == 1, logical(1L))
s012_allNC_M <- s012[indsInfo$hybrid_cat == "XY",..allNC_M]
s012_allNC_M_prod <- apply(s012_allNC_M, 2, prod,na.rm=TRUE)
s012_allNC_M_prod[s012_allNC_M_prod > 1] <- 2
s012_allNC_M_prod_names <- cbind.data.frame("SNP"=names(s012_allNC_M_prod),"NC_M"=s012_allNC_M_prod)


joined_SNPs <- merge(s012_allTX_prod_names,s012_allNC_F_names) %>% merge(s012_allNC_M_prod_names)

joined_SNPs %>% filter((joined_SNPs$TX == 0 & joined_SNPs$NC_F == 2 & joined_SNPs$NC_M == 1) | 
              (joined_SNPs$TX == 2 & joined_SNPs$NC_F == 0 & joined_SNPs$NC_M == 1) )

#joined_SNPs %>% filter((joined_SNPs$TX == 0 | joined_SNPs$TX == 2 )



#log(268435456,2)

