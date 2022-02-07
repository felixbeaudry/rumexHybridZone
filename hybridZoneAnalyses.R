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
library(rEEMSplots) #fail

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
today<-format(Sys.Date(),format="%y%m%d")


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


#write.table(indsInfo$name[indsInfo$hybrid_cat %in% c("XY","Hybrid_W")], file = "GBS.XY.inds", append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE,    col.names = FALSE, qmethod = c("escape", "double"))
#write.table(indsInfo$name[indsInfo$hybrid_cat %in% c("Hybrid_E")], file = "GBS.admx.inds", append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE,    col.names = FALSE, qmethod = c("escape", "double"))
#write.table(indsInfo$name[indsInfo$hybrid_cat %in% c("XYY")], file = "GBS.XYY.inds", append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE,    col.names = FALSE, qmethod = c("escape", "double"))


popSort <- as.data.frame(unique(indsInfo$statePop))
names(popSort) <- c("pop")
pop_pos_gen <- pop_pos[pop_pos$Name %in% popSort$pop]

pop_tally <- indsInfo %>% group_by(statePop) %>%  tally()
pop_pos_gen <- left_join(pop_pos_gen,pop_tally,by=c("Name"="statePop"))
pop_pos_gen <- pop_pos_gen[order(pop_pos_gen$Name),]
pop_pos_gen$strc_popNum <- seq(1,45,1)
pop_pos_gen <- pop_pos_gen[order(pop_pos_gen$mixName),]
pop_pos_gen$constrc_popNum <- seq(1,45,1)

indsInfo <- left_join(indsInfo,pop_pos_gen[,c(8,11,13:14)],by=c("statePop"="Name"))

pheno_pos <- fread('populations_Rumex_hastatulus.csv')

s012 <- fread(paste("GBS.mis60.maf",maf,".noY.012",sep=""),na.strings = "-1")
s012 <- s012[,-1]

#pos
#GBS.mis60.maf01.012.chrom.pos
loc <- fread(paste("GBS.mis60.maf",maf,".noY.012.chrom.pos",sep=""))
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


##Download and load geographic shape files##

#download.file("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/50m/physical/ne_50m_rivers_lake_centerlines.zip",  destfile = 'rivers_lake.zip')
# unzip(zipfile = "rivers_lake.zip", exdir = 'ne-rivers_lake-10m')
# download.file("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/50m/cultural/ne_50m_admin_1_states_provinces.zip",   destfile = 'states_prov.zip')
# unzip(zipfile = "states_prov.zip",  exdir = 'ne-states-50m')
#download.file("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_0_countries.zip",   destfile = 'countries.zip')
#unzip(zipfile = "countries.zip",  exdir = 'ne-countries-50m')

states <- readOGR("ne-states-50m/ne_50m_admin_1_states_provinces.shp")
rivers <- readOGR("ne-rivers_lake-10m/ne_50m_rivers_lake_centerlines.shp")
countries <- readOGR("ne-countries-50m/ne_10m_admin_0_countries.shp")

smith69_hybrids <- readOGR("smith69_hybrids.shp", stringsAsFactors=FALSE)

smith69_hybrids_ll <- spTransform(smith69_hybrids, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
smith69_hybrids_ll_df <- as.data.frame(smith69_hybrids_ll)

####Fig S1####

pdf(paste("fig1_map.hybrids.",today,".pdf",sep=''),width=5.5,height=4)

ggplot() + 
  geom_path(data = countries, aes(x = long, y = lat, group = group),color="gray",size=0.5) +
  geom_path(data = rivers, aes(x = long, y = lat, group = group),color="lightblue",size=0.5) +
  
  geom_point(data = pop_pos[pop_pos$Collec == "Old" & pop_pos$Type == "XY",], aes(x=Longitude,y=Latitude,color=Type),shape=19,color="darkblue")+ #allopatric
  geom_point(data = pop_pos[pop_pos$Collec == "Old" & pop_pos$Type == "XYY",], aes(x=Longitude,y=Latitude),shape=19,color="darkgreen")+ #allopatric
  
  
  geom_point(data = pop_pos[pop_pos$Collec == "New",], aes(x=Longitude,y=Latitude,color=Type),shape=19,color="#ffd700")+ #new sampling sites - alive
  geom_point(data = pop_pos[pop_pos$Collec == "New" & pop_pos$Pop %ni% pop_pos_gen$Pop[pop_pos_gen$Type == "Hybrid"],], aes(x=Longitude,y=Latitude),shape=4)+ #new sampling sites - alive
  
  geom_point(data = smith69_hybrids_ll_df,aes(x=coords.x1,y=coords.x2),shape=19,color="gray") + #Smith hybrids
  
  coord_sf() + 
  xlim(-97.6,-85)+ ylim(28,36.5)+ #all samples
  
  theme_bw(base_size = 12) + guides(color = FALSE) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Longitude",y="Latitude",size="1st flower") +
  scale_color_manual(values=c("gray", "darkblue","darkgreen")) #+  scale_shape_manual(values=c(16,1))

dev.off()


#MLRA
MLRA_shapes <- st_read("nrcs142p2_052440/mlra_v42.shp")

MLRA_shapes_SE = st_crop(MLRA_shapes, xmin=-100, xmax=-75, ymin=28, ymax=38)
MLRA_shapes_SE_top5 <- MLRA_shapes_SE[MLRA_shapes_SE$LRRSYM %in% c("J","N","O","P","T"),]

pdf(paste("figS1_map.MLRA.",today,".pdf",sep=''),width=5.5,height=4)

ggplot() + 
  geom_sf(data = MLRA_shapes_SE_top5,alpha=0.2, lwd=0, color=NA, size=0, aes( fill =as.factor(LRRSYM))) + 
  
  geom_path(data = countries, aes(x = long, y = lat, group = group),color="lightgray",size=0.5) +
   geom_path(data = rivers, aes(x = long, y = lat, group = group),color="lightblue",size=0.5) +
  
  geom_point(data = pop_pos[pop_pos$Collec == "Old" & pop_pos$Type == "XY",], aes(x=Longitude,y=Latitude,color=Type),shape=19,color="darkblue")+ #allopatric
  geom_point(data = pop_pos[pop_pos$Collec == "Old" & pop_pos$Type == "XYY",], aes(x=Longitude,y=Latitude),shape=19,color="darkgreen")+ #allopatric
  
  
  geom_point(data = pop_pos[pop_pos$Collec == "New",], aes(x=Longitude,y=Latitude,color=Type),shape=19,color="#ffd700")+ #new sampling sites - alive
  geom_point(data = pop_pos[pop_pos$Collec == "New" & pop_pos$Pop %ni% pop_pos_gen$Pop[pop_pos_gen$Type == "Hybrid"],], aes(x=Longitude,y=Latitude),shape=4)+ #new sampling sites - alive
  
  geom_point(data = smith69_hybrids_ll_df,aes(x=coords.x1,y=coords.x2),shape=19,color="gray") + #Smith hybrids
  
  coord_sf() + 
  xlim(-98,-76)+ ylim(28,36.5)+ #all samples
  
  theme_bw(base_size = 8) + guides(color = FALSE) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Longitude",y="Latitude",fill="MLRA") +
  scale_color_manual(values=c("gray", "darkblue","darkgreen")) +
  #guides(fill=FALSE) +
  scale_fill_manual(values=c("#ffd700",
                             "#ffb14e",
                             "#fa8775",
                             "#ea5f94",
                             "#9d02d7")) + 
 # theme(
 #   legend.position = c(.08, .70),
 #   legend.justification = c("right", "bottom"),
 #   legend.box.just = "right",
 #   legend.margin = margin(6, 6, 6, 6)
 # )  
  theme(legend.position="bottom")

dev.off()


####Fig 1A PCA####
hybPCA <- PCA(s012, graph = FALSE) #missMDA assumes continuous variables
hybPCAcoord <- cbind.data.frame(hybPCA$ind$coord,indsInfo[,-1])
#hybPCA <- PCA(s012HI[indsInfo$flower =="F",], graph = FALSE) #missMDA assumes continuous variables
#hybPCAcoord <- cbind.data.frame(hybPCA$ind$coord,indsInfo[indsInfo$flower =="F",-1])

hybPCA_plot <- 
ggplot(hybPCAcoord,aes(y=Dim.2*hybPCA$eig[2,2], x=Dim.1*hybPCA$eig[1,2],color=Longitude,shape=Collec)) + 
  geom_point(size=3) + 
  labs(x=paste("PC1 (",round(hybPCA$eig[1,2],digits=2),"%)",sep=""),y=paste("PC2 (",round(hybPCA$eig[2,2],digits=2),"%)",sep=""),color="Long.",shape="Pop.") +
  theme_bw(base_size = 18) + 
  guides(color=FALSE,shape=FALSE)  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
  scale_colour_gradient2(low = "darkblue", high = "darkgreen", mid = "lightgray", midpoint = -91, na.value = NA) +
  theme(
    legend.position = c(.99, .99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )  +  scale_shape_manual(values=c(16,1)) +  coord_sf() 
 #theme(aspect.ratio=1)


cor.test(hybPCAcoord$Dim.1, hybPCAcoord$Longitude, method=c("pearson"))

####STRUCTURE####
strcPop <- indsInfo[order(indsInfo$mixName),c(1,16)] #fread('HiC_hyb.pop.num')


for (k in c(2:3)){
  kfileName <-  paste('K',k,'.ind.clumpp',sep='')
  Kfile <- fread(kfileName)
  for(j in c(1:5)){
    Kfile <- Kfile[,-1]
  }
  indpopsk <- cbind.data.frame(strcPop,Kfile)
  nameVector <- c("inds","Pop")
  for(i in (1:k)){
    nameVector <- c(nameVector,paste("k",i,sep=''))
  }
  names(indpopsk) <- nameVector
  
  indpopsk.melt <- melt(indpopsk,id.vars = c("inds","Pop"),verbose=FALSE)
  if(k==2){
    indpopsk.melt.ks <- cbind.data.frame(indpopsk.melt,"k"=k)
  }
  else{
    indpopsk.melt.ks <- 
      rbind.data.frame(indpopsk.melt.ks,
                       cbind.data.frame(indpopsk.melt,"k"=k)
      )
  }
  
}
indpopskmeltks <- indpopsk.melt.ks
strc_Pop <- sqldf('select indpopskmeltks.*, indsInfo.* from indpopskmeltks
                         left join indsInfo on indpopskmeltks.inds = indsInfo.name')
strc_Pop <- strc_Pop[,c(1,3:5,8,10,12:14)]

strc_Pop <- strc_Pop[order(strc_Pop$Longitude),]

strc_Pop$Statepop <- factor(strc_Pop$statePop,levels=c(unique(strc_Pop$statePop)))

#strc_Pop$var_alt <- strc_Pop$variable

strc_Pop$value_flip <- strc_Pop$value
strc_Pop$value_flip[strc_Pop$k == 2] <- 1 - strc_Pop$value[strc_Pop$k == 2]

strc_Pop$var_alt[strc_Pop$k == 2 & strc_Pop$variable == "k1" ] <- "TX"
strc_Pop$var_alt[strc_Pop$k == 2 & strc_Pop$variable == "k2" ] <- "NC"

strc_Pop$var_alt[strc_Pop$k == 3 & strc_Pop$variable == "k2" ] <- "TXl"
strc_Pop$var_alt[strc_Pop$k == 3 & strc_Pop$variable == "k3" ] <- "TXh"
strc_Pop$var_alt[strc_Pop$k == 3 & strc_Pop$variable == "k1" ] <- "NC"

####Fig S3####
pdf(paste("figS3_structure.",today,".pdf",sep=''),width=5.5,height=4)

ggplot(strc_Pop,aes(x=inds,y=value_flip,fill=var_alt)) +
  geom_bar( stat="identity" ) +
  labs(x="Population",y="Hybrid Index",fill="") +
  facet_grid(as.factor(k) ~ Statepop,scales = "free",space="free")  +
  scale_fill_manual(values=c( 
    "darkgreen", "lightblue",  "lightblue","grey"))+ 
  guides(color=FALSE,fill=FALSE) + 
  theme(panel.spacing = unit(0, "lines"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text.x = element_blank(),
        panel.border = element_blank(),
        text = element_text(size=18),
        panel.background = element_rect(fill = "white")) 

dev.off()

#strc_Pop_k2 <- 
ggplot(strc_Pop %>% filter(k==2),aes(x=inds,y=value_flip,fill=var_alt)) +
  geom_bar( stat="identity" ) +
  labs(x="Individuals grouped by Population, ordered East to West",y="Hybrid Index",fill="") +
  facet_grid(as.factor(k) ~ Statepop,scales = "free",space="free")  +
  scale_fill_manual(values=c( 
    "darkgreen", "lightblue",  "lightblue","grey"))+ 
  guides(color=FALSE,fill=FALSE) + 
  theme(panel.spacing = unit(0, "lines"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text.x = element_blank(),
        panel.border = element_blank(),
        text = element_text(size=18),
        panel.background = element_rect(fill = "white")) 

#PCA_STRC <- plot_grid(  hybPCA_plot, strc_Pop_k2,  labels = c('A', 'B'), label_size = 12,ncol = 1, align = 'vh',axis='tbrl')


#HI distribution fit 


PC_strc <- left_join(hybPCAcoord,strc_Pop[strc_Pop$k == 2 & strc_Pop$variable == "k1",])

cor.test(PC_strc$value, PC_strc$Dim.1, method=c("pearson"))

pc1_hi_plot <- 
  ggplot(PC_strc,aes(x=value,y=Dim.1,color=Longitude)) + geom_point(size=2) + 
  theme_bw(base_size = 8) + 
  #guides(color=FALSE)+ 
  labs(y="Individual PC1",x="Hybrid Index (CONSTRUCT)") +
  scale_colour_gradient2(low = "darkblue", high = "darkgreen", mid = "lightgray", midpoint = -91, na.value = NA) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(
    legend.position = c(.05, .95),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(6, 6, 6, 6)
  ) 

pdf(paste("figS4_HI2PC_A.",today,".pdf",sep=''),width=4,height=4)

ggExtra::ggMarginal(pc1_hi_plot, type = "histogram")
dev.off()

hybrid_pops <-  pop_pos$Name[pop_pos$Type == "Hybrid" & pop_pos$riverside == "E"]
hybrid_pops_HI <- strc_Pop[strc_Pop$statePop %in% hybrid_pops,]

min(hybrid_pops_HI$value[hybrid_pops_HI$k == 2 & hybrid_pops_HI$variable == "k1"])
max(hybrid_pops_HI$value[hybrid_pops_HI$k == 2 & hybrid_pops_HI$variable == "k1"])

ggplot(hybrid_pops_HI[hybrid_pops_HI$k == 2 & hybrid_pops_HI$variable == "k1",],aes(x=value)) +
  geom_density() + geom_histogram(aes(y=..density..)) + theme_bw()

####HIest####
s012HIest <- replace(s012,is.na(s012),-9)
# parental allele frequencies (assumed diagnostic)

BS.est <-HIC(s012HIest)

#BS.test.ind <- cbind.data.frame(BS.est,indsInfo[indsInfo$flower == "F",])
BS.test.ind <- cbind.data.frame(BS.est,indsInfo)

dt.triangle <- data.table(group = c(1,1,1), polygon.x = c(0,0.5,1), polygon.y = c(0,1,0))

pdf(paste("figS3_fullTriangle.",today,".pdf",sep=''),width=4,height=4)

#triangle_plot <- 
ggplot() +
  geom_polygon(
    data = dt.triangle
    ,aes(
      x=polygon.x
      ,y=polygon.y
      ,group=group
    ), fill="white", color="black"
  ) +
  geom_point(data=BS.test.ind,aes(x=S,y=H,color=Longitude,shape=Collec),size=0.75) + 
  #  geom_text(data=BS.test.ind,aes(x=S,y=H,color=Longitude,label=name),size=3) + 
  
  xlim(0,1) + ylim(0,1)+ theme(aspect.ratio=1) +
  theme_bw(base_size = 8) + #guides(shape=FALSE)+
  labs(x="Ancestry",y="Heterozygosity",shape="Sex") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_colour_gradient2(low = "darkblue", high = "darkgreen", mid = "lightgray", midpoint = -91, na.value = NA) +
 # scale_shape_manual(values=c(16,18,15)) +
  theme(
    legend.position = c(.99, .99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)  
  ) +  scale_shape_manual(values=c(16,1))
dev.off()

# parental allele frequencies (assumed diagnostic)
BS.P <- data.frame(Locus=names(s012HIest),Allele="XY",P1=0,P2=2)
# calculate likelihoods for early generation hybrid classes
BS.class <- HIclass(s012HIest,BS.P,type="allele.count")
# compare classification with maximum likelihood estimates
BS.test <- HItest(BS.class,BS.est,thresholds=c(2,2))

##classification efficacy tests##
table(BS.test$c1)
# x are TRUE, meaning the best classification is at least x log-likelihood units better than the next best
table(BS.test$c2)
# x are TRUE, meaning the MLE S and H are within x log-likelihood units of the best classification, i.e., the simple classification is rejected in all but x cases
table(BS.test$Best.class,BS.test$c2)
# individuals were classified as F2-like (class 3) or backcross to CTS (class 4), but only x of the F2's were credible
BS.test[BS.test$c2,]
# in only x case was the F2 classification a better fit (based on AIC) than the continuous model.
# equivalent to the AIC criterion:
BS.test <- HItest(BS.class,BS.est,thresholds=c(2,1))



loc_hybPCA  <- cbind.data.frame(loc,hybPCA$var$contrib)
ggplot(loc_hybPCA,aes(x=mb,y=Dim.1)) + geom_point() + facet_grid(. ~ LG, space="free",scales = "free")
ggplot(loc_hybPCA,aes(x=mb,y=Dim.2)) + geom_point() + facet_grid(. ~ LG, space="free",scales = "free")
ggplot(loc_hybPCA,aes(x=mb,y=Dim.3)) + geom_point() + facet_grid(. ~ LG, space="free",scales = "free")
ggplot(loc_hybPCA,aes(x=mb,y=Dim.4)) + geom_point() + facet_grid(. ~ LG, space="free",scales = "free")


loc_hybPCA_order <- loc_hybPCA[order(loc_hybPCA$Dim.1,decreasing = T),]

ggplot(loc_hybPCA_order[c(1:100),],aes(x=chrom)) + geom_bar()

loc_hybPCA_order[c(1:100),] %>% group_by(chrom) %>% tally()
loc_hybPCA_order %>% group_by(chrom) %>% tally()

#(2000 + 1458)/13613
#2000/13613

dim1_SNPs <- loc_hybPCA_order$SNP[c(1:100)]

dim1_SNPs <- dim1_SNPs+1
dim1_SNPs <- paste("V",dim1_SNPs,sep="")

#s012HI_top <- s012HI[indsInfo$flower == "F",..dim1_SNPs]
s012HI_top <- s012[,..dim1_SNPs]

s012HIest <- replace(s012HI_top,is.na(s012HI_top),-9)
# parental allele frequencies (assumed diagnostic)

BS.est <-HIC(s012HIest)

BS.test.ind <- cbind.data.frame(BS.est,indsInfo)

#triangle_plot <- 
ggplot() +
  geom_polygon(
    data = dt.triangle
    ,aes(
      x=polygon.x
      ,y=polygon.y
      ,group=group
    ), fill="white", color="black"
  ) + theme(aspect.ratio=1) +
  geom_point(data=BS.test.ind,aes(x=S,y=H,color=Longitude,shape=Collec),size=3) + 
#  geom_text(data=BS.test.ind,aes(x=S,y=H,color=Longitude,label=name),size=3) + 
  
  xlim(0,1) + ylim(0,1)+
  theme_bw(base_size = 18) + #guides(shape=FALSE)+
  labs(x="Ancestry",y="Heterozygosity",shape="Sex") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_colour_gradient2(low = "darkblue", high = "darkgreen", mid = "lightgray", midpoint = -91, na.value = NA) +
  
  theme(
    legend.position = c(.99, .99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  ) +  scale_shape_manual(values=c(16,1)) +  coord_sf() #theme(aspect.ratio=1) #
#



# parental allele frequencies (assumed diagnostic)
BS.P <- data.frame(Locus=names(s012HI_top),Allele="XY",P1=0,P2=2)
# calculate likelihoods for early generation hybrid classes
BS.class <- HIclass(s012HI_top,BS.P,type="allele.count")
# compare classification with maximum likelihood estimates
BS.test <- HItest(BS.class,BS.est,thresholds=c(2,2))

##classification efficacy tests##
#table(BS.test$c1)
# x are TRUE, meaning the best classification is at least x log-likelihood units better than the next best
#table(BS.test$c2)
# x are TRUE, meaning the MLE S and H are within x log-likelihood units of the best classification, i.e., the simple classification is rejected in all but x cases
#table(BS.test$Best.class,BS.test$c2)
# individuals were classified as F2-like (class 3) or backcross to CTS (class 4), but only x of the F2's were credible
#BS.test[BS.test$c2,]
# in only x case was the F2 classification a better fit (based on AIC) than the continuous model.
# equivalent to the AIC criterion:
#BS.test <- HItest(BS.class,BS.est,thresholds=c(2,1))

hybPCA_plot <- 
 # ggplot(hybPCAcoord,aes(y=Dim.2*hybPCA$eig[2,2], x=Dim.1*hybPCA$eig[1,2],color=Longitude,shape=Collec)) + 
   ggplot(hybPCAcoord,aes(y=Dim.2, x=Dim.1,color=Longitude,shape=Collec)) + 
  
    geom_point(size=2) + 
  labs(x=paste("PC1 (",round(hybPCA$eig[1,2],digits=2),"%)",sep=""),y=paste("PC2 (",round(hybPCA$eig[2,2],digits=2),"%)",sep=""),color="Long.",shape="Pop.") +
  theme_bw(base_size = 8) + 
  guides(color=FALSE,shape=FALSE)  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
  scale_colour_gradient2(low = "darkblue", high = "darkgreen", mid = "lightgray", midpoint = -91, na.value = NA) +
  theme(
    legend.position = c(.99, .99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )  +  scale_shape_manual(values=c(16,1)) +  coord_sf() 

strc_Pop_k2 <- 
  ggplot(strc_Pop %>% filter(k==2),aes(x=inds,y=value_flip,fill=var_alt)) +
  geom_bar( stat="identity" ) +
  labs(x="Individuals grouped by Population, ordered East to West",y="Hybrid Index",fill="") +
  facet_grid(as.factor(k) ~ Statepop,scales = "free",space="free")  +
  scale_fill_manual(values=c( 
    "darkgreen", "lightblue",  "lightblue","grey"))+ 
  guides(color=FALSE,fill=FALSE) + 
  theme(panel.spacing = unit(0, "lines"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text.x = element_blank(),
        panel.border = element_blank(),
        text = element_text(size=8),
        panel.background = element_rect(fill = "white")) 

PCA_STRC <- plot_grid(  hybPCA_plot, strc_Pop_k2,  labels = c('A', 'B'), label_size = 12,ncol = 1, align = 'vh',axis='tbrl')

BS.test.ind$Collec[BS.test.ind$Collec == "New"] <- "Sympatric"
BS.test.ind$Collec[BS.test.ind$Collec == "Old"] <- "Allopatric"

triangle_plot <- 
ggplot() +
  geom_polygon(
    data = dt.triangle
    ,aes(
      x=polygon.x
      ,y=polygon.y
      ,group=group
    ), fill="white", color="black"
  ) + theme(aspect.ratio=1) +
  geom_point(data=BS.test.ind,aes(x=S,y=H,color=Longitude,shape=Collec),size=2) + 
  #  geom_text(data=BS.test.ind,aes(x=S,y=H,color=Longitude,label=name),size=3) + 
  
  xlim(0,1) + ylim(0,1)+
  theme_bw(base_size = 8) + #guides(shape=FALSE)+
  labs(x="Ancestry",y="Heterozygosity",shape="Population") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_colour_gradient2(low = "darkblue", high = "darkgreen", mid = "lightgray", midpoint = -91, na.value = NA) +
  
  theme(
    legend.position = c(.99, .99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  ) +  scale_shape_manual(values=c(1,16)) +  coord_sf()

pdf(paste("fig2_individuals.",today,".pdf",sep=''),width=7,height=4)

plot_grid(  hybPCA_plot, triangle_plot,  labels = c('A', 'B'), label_size = 12,ncol = 2, align = 'vh',axis='tbrl')

#plot_grid(  hybPCA_plot, triangle_plot,  labels = c('A', 'B'), label_size = 12,ncol = 2)

dev.off()



####population level clustering####
constrcPop <-  indsInfo[order(indsInfo$mixName),c(1,17)]  
#constrcPop <- fread("GBS.pop.num")
#constrcPop_f <- constrcPop[indsInfo$flower == "F",] #loss of BRAN & WEIR (no females)

#constrcPopf_num <- cbind.data.frame("og"=unique(constrcPop$V2),"new"=seq(1,length(unique(constrcPop$V2)),1))
#constrcPop_f <- left_join(constrcPop_f,constrcPopf_num,by=c("V2"="og"))

#constrcPop_f$V1
pop.data.matrix <- matrix(NA,nrow=length(unique(constrcPop$constrc_popNum)),ncol=ncol(s012))
for(i in 1:nrow(pop.data.matrix)){
  pop.data.matrix[i,] <- colMeans(
    s012[
      which(constrcPop$constrc_popNum==i),,
      drop=FALSE
      ],na.rm=TRUE
  )/2
}

pop.data.matrix[is.nan(pop.data.matrix)] <- NA

pop_pca <- PCA(pop.data.matrix, graph = FALSE) #missMDA assumes continuous variables
#pop_pos_genf <- pop_pos_gen[pop_pos_gen$Pop %ni% c("BRAN" , "WEIR"),]
 # pop_pos_genf$mixName[order(pop_pos_genf$mixName) ]
pop_PCAcoord <- cbind.data.frame(pop_pca$ind$coord,pop_pos_gen[order(pop_pos_gen$mixName) ,])

pop_PCAcoord$pc1_norm <- (pop_PCAcoord$Dim.1 - min(pop_PCAcoord$Dim.1))/ (max(pop_PCAcoord$Dim.1) - min(pop_PCAcoord$Dim.1))
pop_PCAcoord$pc1_oppo <- 1 - pop_PCAcoord$pc1_norm


#poppc1_long_plot <- 
  ggplot(pop_PCAcoord,aes(x=Longitude, y=pc1_norm)) + 
 # ggplot(pop_PCAcoord,aes(x=Dim.1, y=Dim.2,color=Type)) + 
  geom_rect( mapping=aes(xmin=-91, xmax=-88.5, ymin=0, ymax=1), fill="lightblue",color=NA, alpha=0.5) +
  
  geom_point(size=3) + 
#   geom_text(aes(label=Pop)) +
    geom_smooth(method = "lm", se=FALSE, color="black") +
    guides(color=FALSE) +
  labs(x="Longitude",y="PC1 Normalized") +
  theme_bw(base_size = 18) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_colour_gradient2(low = "darkblue", high = "darkgreen", mid = "lightgray", midpoint = -90, na.value = NA) +
  xlim(-98,-76)

cor.test(pop_PCAcoord$Longitude, pop_PCAcoord$pc1_norm, method=c("pearson"))


####CONSTRUCT####
  # how to run a conStruct analysis vignette(topic="run-conStruct",package="conStruct")
  latlon_sort <- pop_pos_gen[order(pop_pos_gen$mixName),]
  
  #for(i in   1:6  ){
  #  my.run <- conStruct(spatial = FALSE, #non-spatial
  #                      K = as.numeric(i), 
  #                      freqs = pop.data.matrix,
  #                      geoDist = as.matrix(fields::rdist.earth(x1 =  latlon_sort[,c(2,1)])), 
  #                      coords = as.matrix(latlon_sort[,c(2,1)]),
  #                      prefix = paste("nspk",i,sep=""),
  #                      make.figs=FALSE)
  #}
  
  #for(i in   1:6){
  #  my.run <- conStruct(spatial = TRUE, #spatial
  #                      K = as.numeric(i), 
  #                      freqs = pop.data.matrix,
  #                      geoDist = as.matrix(fields::rdist.earth(x1 =  latlon_sort[,c(2,1)])), 
  #                      coords = as.matrix(latlon_sort[,c(2,1)]),
  #                      prefix = paste("spk",i,sep=""),
  #                      make.figs=FALSE)
  #}
  
  # how to visualize the output of a conStruct model
  #vignette(topic="visualize-results",package="conStruct")
  
  load("nspK3_conStruct.results.Robj")
  load("nspK3_data.block.Robj")
  
  admix.props <- as.data.frame(conStruct.results$chain_1$MAP$admix.proportions)
  names(admix.props) <- c("K1","K2","K3")
  
  admix.melt <-
    melt(admix.props,verbose=FALSE)
  
  admix.props.df <- cbind.data.frame( admix.props,pop_pos_gen[order(pop_pos_gen$mixName),])
  admix.melt.df <- cbind.data.frame( admix.melt,pop_pos_gen[order(pop_pos_gen$mixName),])
  
  admix.melt.df$factorName = factor(admix.melt.df$mixName, levels=pop_pos_gen$mixName[order(pop_pos_gen$Longitude)])
  
construct_pie_nsp_plot <- 
  ggplot() + 
    geom_rect( mapping=aes(xmin=-91, xmax=-88.5, ymin=30, ymax=36), color=NA, alpha=0.3) +
    
    geom_path(data = countries, aes(x = long, y = lat, group = group),color="lightgray") +
    geom_path(data = rivers, aes(x = long, y = lat, group = group),color="lightblue") +
 #   xlim(-94.5,-84.5)+ ylim(30,34)+ #all samples
     xlim(-98,-76)+ ylim(30,36)+ #all samples
    
    geom_scatterpie(data= admix.props.df, cols=c("K1","K2","K3") , aes(x=Longitude, y=Latitude, r = 0.2),color=NA) +
    
    theme_bw(base_size = 18) + guides(fill = FALSE) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(x="",y="Latitude",fill="") +
    scale_fill_manual(values=c( "darkgreen","darkblue","gray")) 
    
  plot_grid(  construct_pie_nsp_plot, construct_pie_sp_plot,  labels = c('A', 'B'), label_size = 12,ncol = 1, align = 'vh',axis='tbrl')
    

  # how to compare and select between different conStruct models
  #vignette(topic="model-comparison",package="conStruct")
  
  # Loop through output files generated by conStruct runs with K=1 through 5 and calculate the   layer contributions for each layer in each run  
  
  layer.contributions <- matrix(NA,nrow=6,ncol=6)
  
  # load the conStruct.results.Robj and data.block.Robj
  #   files saved at the end of a conStruct run
  load("nspK1_conStruct.results.Robj")
  load("nspK1_data.block.Robj")
  
  # calculate layer contributions
  layer.contributions[,1] <- c(calculate.layer.contribution(conStruct.results[[1]],data.block),rep(0,5))
  tmp <- conStruct.results[[1]]$MAP$admix.proportions
  
  for(i in 2:6){
    # load the conStruct.results.Robj and data.block.Robj
    #   files saved at the end of a conStruct run
    load(sprintf("nspK%s_conStruct.results.Robj",i))
    load(sprintf("nspK%s_data.block.Robj",i))
    
    # match layers up across runs to keep plotting colors consistent
    #   for the same layers in different runs
    tmp.order <- match.layers.x.runs(tmp,conStruct.results[[1]]$MAP$admix.proportions)  
    
    # calculate layer contributions
    layer.contributions[,i] <- c(calculate.layer.contribution(conStruct.results=conStruct.results[[1]],
                                                              data.block=data.block,
                                                              layer.order=tmp.order),
                                 rep(0,6-i))
    tmp <- conStruct.results[[1]]$MAP$admix.proportions[,tmp.order]
  }
  
  ####Fig S3: CONSTRUCT model selection###
  pdf(paste("figS7D_construct_partition_nsp.",today,".pdf",sep=''),width=7,height=4)
  

  barplot(layer.contributions,
          #col=c("blue", "red", "goldenrod1", "forestgreen", "darkorchid1"),
          xlab="",
          ylab="layer contributions",
          names.arg=paste0("K=",1:6))

 dev.off()
  

  
#join PCA to CONSTRUCT
popPCA_contrst <- left_join(admix.props.df,pop_PCAcoord[,c(1:5,16)])

pca_K2_plot <- 
  ggplot(popPCA_contrst,aes(x=K2,y=Dim.1,color=Longitude)) + 
    geom_point(size=2) + 
 #   geom_text(aes(label=Pop))+
  theme_bw(base_size = 8) + 
  #guides(color=FALSE)+ 
  labs(y="PC1",x="Hybrid Index (CONSTRUCT)") +
  scale_colour_gradient2(low = "darkblue", high = "darkgreen", mid = "lightgray", midpoint = -91, na.value = NA) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(
    legend.position = c(.05, .95),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(6, 6, 6, 6)
  ) 

pdf(paste("figS4_HI2PC_B.",today,".pdf",sep=''),width=4,height=4)

ggExtra::ggMarginal(pca_K2_plot, type = "histogram")
dev.off()


cor.test(popPCA_contrst$Dim.1, popPCA_contrst$K2, method=c("pearson"))





####HZAR####


cnstrc_dist <-  diag(distm(cbind(min(admix.props.df$Longitude),admix.props.df$Latitude),cbind(admix.props.df$Longitude,admix.props.df$Latitude)))

cnstrc_hzar_table <- cbind.data.frame("Locality_ID"=admix.props.df$Pop,"Locality"=admix.props.df$State,"dist"=cnstrc_dist,"strc1"=admix.props.df$K2,"strc2"=admix.props.df$K1,"Nsamples"=admix.props.df$n)

cnstrc_cline <- hzar.doMolecularData1DPops(cnstrc_hzar_table$dist, cnstrc_hzar_table$strc1, cnstrc_hzar_table$Nsamples)
hzar.plot.obsData(cnstrc_cline)

##start fitting
chainLength=1e5;                       

## Make each model run off a separate seed
mainSeed=
  list(A=c(596,528,124,978,544,99),
       B=c(528,124,978,544,99,596),
       C=c(124,978,544,99,596,528))

## Blank out space in memory to hold molecular analysis
if(length(apropos("^cnstrc_cline$",ignore.case=FALSE)) == 0 ||
   !is.list(cnstrc_cline) ) cnstrc_cline <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
cnstrc_cline$WGS <- list();
## Space to hold the observed data
cnstrc_cline$WGS$obs <- list();
cnstrc_cline$WGS$obs <- hzar.doMolecularData1DPops(cnstrc_hzar_table$dist, cnstrc_hzar_table$strc1, cnstrc_hzar_table$Nsamples)

## Space to hold the models to fit
cnstrc_cline$WGS$models <- list();
## Space to hold the compiled fit requests
cnstrc_cline$WGS$fitRs <- list();
## Space to hold the output data chains
cnstrc_cline$WGS$runs <- list();
## Space to hold the analysed data
cnstrc_cline$WGS$analysis <- list();

cnstrc_cline.loadWGSmodel <- function(scaling,tails,
                                  id=paste(scaling,tails,sep="."))
  cnstrc_cline$WGS$models[[id]] <<- hzar.makeCline1DFreq(cnstrc_cline$WGS$obs, scaling, tails)

cnstrc_cline.loadWGSmodel("fixed","none","modelI");
cnstrc_cline.loadWGSmodel("fixed" ,"mirror","modelII");
cnstrc_cline.loadWGSmodel("fixed","both","modelIII");

cnstrc_cline.loadWGSmodel("free" ,"none","modelIV");
cnstrc_cline.loadWGSmodel("free" ,"mirror","modelV");
cnstrc_cline.loadWGSmodel("free" ,"both","modelVI");

print(cnstrc_cline$WGS$models)

## Compile each of the models to prepare for fitting
cnstrc_cline$WGS$fitRs$init <- sapply(cnstrc_cline$WGS$models,
                                  hzar.first.fitRequest.old.ML,
                                  obsData=cnstrc_cline$WGS$obs,
                                  verbose=FALSE,
                                  simplify=FALSE)

## Run just one of the models for an initial chain
cnstrc_cline$WGS$runs$init <- list()
cnstrc_cline$WGS$runs$init$modelI <-  hzar.doFit(cnstrc_cline$WGS$fitRs$init$modelI)
cnstrc_cline$WGS$runs$init$modelII <-  hzar.doFit(cnstrc_cline$WGS$fitRs$init$modelII)
cnstrc_cline$WGS$runs$init$modelIII <-  hzar.doFit(cnstrc_cline$WGS$fitRs$init$modelIII)
cnstrc_cline$WGS$runs$init$modelIV <-  hzar.doFit(cnstrc_cline$WGS$fitRs$init$modelIV)
cnstrc_cline$WGS$runs$init$modelV <-  hzar.doFit(cnstrc_cline$WGS$fitRs$init$modelV)
cnstrc_cline$WGS$runs$init$modelVI <-  hzar.doFit(cnstrc_cline$WGS$fitRs$init$modelVI)

## Plot the trace
plot(hzar.mcmc.bindLL(cnstrc_cline$WGS$runs$init$modelI))
plot(hzar.mcmc.bindLL(cnstrc_cline$WGS$runs$init$modelII))
plot(hzar.mcmc.bindLL(cnstrc_cline$WGS$runs$init$modelIII))
plot(hzar.mcmc.bindLL(cnstrc_cline$WGS$runs$init$modelIV))
plot(hzar.mcmc.bindLL(cnstrc_cline$WGS$runs$init$modelV))
plot(hzar.mcmc.bindLL(cnstrc_cline$WGS$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
cnstrc_cline$WGS$fitRs$chains <-
  lapply(cnstrc_cline$WGS$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
cnstrc_cline$WGS$fitRs$chains <-
  hzar.multiFitRequest(cnstrc_cline$WGS$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Just to be thorough, randomize the initial value for each fit
for(chain in c(1:18)){
  cnstrc_cline$WGS$fitRs$chains[[chain]]$modelParam$init["center"] <- runif(1,0,1809877)
  cnstrc_cline$WGS$fitRs$chains[[chain]]$modelParam$init["width"] <- runif(1,0,1809877)
}

#free freqs
for(chain in c(10:18)){
  cnstrc_cline$WGS$fitRs$chains[[chain]]$modelParam$init["pMin"] <- runif(1,0,1)
  cnstrc_cline$WGS$fitRs$chains[[chain]]$modelParam$init["pMax"] <- runif(1,0,1)
}

#mirror tails
for(chain in c(4:6,13:15)){
  cnstrc_cline$WGS$fitRs$chains[[chain]]$modelParam$init["deltaM"] <- runif(1,0,1809877)
  cnstrc_cline$WGS$fitRs$chains[[chain]]$modelParam$init["tauM"] <- runif(1,0,1)
}

#both tails
for(chain in c(7:9,16:18)){
  cnstrc_cline$WGS$fitRs$chains[[chain]]$modelParam$init["deltaL"] <- runif(1,0,1809877)
  cnstrc_cline$WGS$fitRs$chains[[chain]]$modelParam$init["deltaR"] <- runif(1,0,1809877)
  
  cnstrc_cline$WGS$fitRs$chains[[chain]]$modelParam$init["tauL"] <- runif(1,0,1)
  cnstrc_cline$WGS$fitRs$chains[[chain]]$modelParam$init["tauR"] <- runif(1,0,1)
  
}

cnstrc_cline$WGS$runs$chains <-  hzar.doChain.multi(cnstrc_cline$WGS$fitRs$chains,
                                                doPar=TRUE,
                                                inOrder=FALSE,
                                                count=3)

today<-format(Sys.Date(),format="%d%b%Y")
save(cnstrc_cline,file=paste("cnstrc_cline_WGS_",today,".rdata",sep=''))
load('cnstrc_cline_WGS_09Sep2021.rdata')

#did the runs converge?
summary(do.call(mcmc.list,
                lapply(cnstrc_cline$WGS$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )
summary(do.call(mcmc.list,
                lapply(cnstrc_cline$WGS$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )


summary(do.call(mcmc.list,
                lapply(cnstrc_cline$WGS$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )
summary(do.call(mcmc.list,
                lapply(cnstrc_cline$WGS$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )
summary(do.call(mcmc.list,
                lapply(cnstrc_cline$WGS$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )
summary(do.call(mcmc.list,
                lapply(cnstrc_cline$WGS$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Start aggregation of data for analysis
## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
cnstrc_cline$WGS$analysis$initDGs <- list(  nullModel =  hzar.dataGroup.null(cnstrc_cline$WGS$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
cnstrc_cline$WGS$analysis$initDGs$modelI <-  hzar.dataGroup.add(cnstrc_cline$WGS$runs$init$modelI)
cnstrc_cline$WGS$analysis$initDGs$modelII <- hzar.dataGroup.add(cnstrc_cline$WGS$runs$init$modelII)
cnstrc_cline$WGS$analysis$initDGs$modelIII <- hzar.dataGroup.add(cnstrc_cline$WGS$runs$init$modelIII)
cnstrc_cline$WGS$analysis$initDGs$modelIV <- hzar.dataGroup.add(cnstrc_cline$WGS$runs$init$modelIV)
cnstrc_cline$WGS$analysis$initDGs$modelV <- hzar.dataGroup.add(cnstrc_cline$WGS$runs$init$modelV)
cnstrc_cline$WGS$analysis$initDGs$modelVI <- hzar.dataGroup.add(cnstrc_cline$WGS$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created
cnstrc_cline$WGS$analysis$oDG <-  hzar.make.obsDataGroup(cnstrc_cline$WGS$analysis$initDGs)
cnstrc_cline$WGS$analysis$oDG <- hzar.copyModelLabels(cnstrc_cline$WGS$analysis$initDGs,
                                                      cnstrc_cline$WGS$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
cnstrc_cline$WGS$analysis$oDG <-  hzar.make.obsDataGroup(lapply(cnstrc_cline$WGS$runs$chains,
                                                            hzar.dataGroup.add),
                                                         cnstrc_cline$WGS$analysis$oDG);

## Do model selection based on the AICc scores
print(cnstrc_cline$WGS$analysis$AICcTable <- hzar.AICc.hzar.obsDataGroup(cnstrc_cline$WGS$analysis$oDG));

## Print out the model with the minimum AICc score
print(cnstrc_cline$WGS$analysis$model.name <-
        rownames(cnstrc_cline$WGS$analysis$AICcTable
        )[[ which.min(cnstrc_cline$WGS$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
cnstrc_cline$WGS$analysis$model.selected <- cnstrc_cline$WGS$analysis$oDG$data.groups[[cnstrc_cline$WGS$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(cnstrc_cline$WGS$analysis$model.selected,   names(cnstrc_cline$WGS$analysis$model.selected$data.param)));

(( 684713.9 /1809877)*19.99)-96.86
(( 752159.5 /1809877)*19.99)-96.86
(( 844511.1 /1809877)*19.99)-96.86


## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(cnstrc_cline$WGS$analysis$model.selected))


## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(cnstrc_cline$WGS$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(cnstrc_cline$WGS$analysis$model.selected);

###PCA cline fit###
cnstrc_hzar_table <- cbind.data.frame("Locality_ID"=admix.props.df$Pop,"Locality"=admix.props.df$State,"dist"=cnstrc_dist,"strc1"=admix.props.df$K2,"strc2"=admix.props.df$K1,"Nsamples"=admix.props.df$n)


pca_dist <-  diag(distm(cbind(min(pop_PCAcoord$Longitude),pop_PCAcoord$Latitude),cbind(pop_PCAcoord$Longitude,pop_PCAcoord$Latitude)))

pca_hzar_table <- cbind.data.frame("Locality_ID"=pop_PCAcoord$Pop,"Locality"=pop_PCAcoord$State,"dist"=pca_dist,"strc1"=pop_PCAcoord$pc1_norm,"strc2"=pop_PCAcoord$pc1_norm,"Nsamples"=pop_PCAcoord$n)

pca_cline <- hzar.doMolecularData1DPops(pca_hzar_table$dist, pca_hzar_table$strc1, pca_hzar_table$Nsamples)
hzar.plot.obsData(pca_cline)

##start fitting
chainLength=1e5;                       

## Make each model run off a separate seed
mainSeed=
  list(A=c(596,528,124,978,544,99),
       B=c(528,124,978,544,99,596),
       C=c(124,978,544,99,596,528))

## Blank out space in memory to hold molecular analysis
if(length(apropos("^pcacline$",ignore.case=FALSE)) == 0 ||
   !is.list(pcacline) ) pcacline <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
pcacline$WGS <- list();
## Space to hold the observed data
pcacline$WGS$obs <- list();
pcacline$WGS$obs <- hzar.doMolecularData1DPops(pca_hzar_table$dist, pca_hzar_table$strc1, pca_hzar_table$Nsamples)

## Space to hold the models to fit
pcacline$WGS$models <- list();
## Space to hold the compiled fit requests
pcacline$WGS$fitRs <- list();
## Space to hold the output data chains
pcacline$WGS$runs <- list();
## Space to hold the analysed data
pcacline$WGS$analysis <- list();

pcacline.loadWGSmodel <- function(scaling,tails,
                                  id=paste(scaling,tails,sep="."))
  pcacline$WGS$models[[id]] <<- hzar.makeCline1DFreq(pcacline$WGS$obs, scaling, tails)

pcacline.loadWGSmodel("fixed","none","modelI");
pcacline.loadWGSmodel("fixed" ,"mirror","modelII");
pcacline.loadWGSmodel("fixed","both","modelIII");

pcacline.loadWGSmodel("free" ,"none","modelIV");
pcacline.loadWGSmodel("free" ,"mirror","modelV");
pcacline.loadWGSmodel("free" ,"both","modelVI");

print(pcacline$WGS$models)

## Compile each of the models to prepare for fitting
pcacline$WGS$fitRs$init <- sapply(pcacline$WGS$models,
                                  hzar.first.fitRequest.old.ML,
                                  obsData=pcacline$WGS$obs,
                                  verbose=FALSE,
                                  simplify=FALSE)

## Run just one of the models for an initial chain
pcacline$WGS$runs$init <- list()
pcacline$WGS$runs$init$modelI <-  hzar.doFit(pcacline$WGS$fitRs$init$modelI)
pcacline$WGS$runs$init$modelII <-  hzar.doFit(pcacline$WGS$fitRs$init$modelII)
pcacline$WGS$runs$init$modelIII <-  hzar.doFit(pcacline$WGS$fitRs$init$modelIII)
pcacline$WGS$runs$init$modelIV <-  hzar.doFit(pcacline$WGS$fitRs$init$modelIV)
pcacline$WGS$runs$init$modelV <-  hzar.doFit(pcacline$WGS$fitRs$init$modelV)
pcacline$WGS$runs$init$modelVI <-  hzar.doFit(pcacline$WGS$fitRs$init$modelVI)

## Plot the trace
plot(hzar.mcmc.bindLL(pcacline$WGS$runs$init$modelI))
plot(hzar.mcmc.bindLL(pcacline$WGS$runs$init$modelII))
plot(hzar.mcmc.bindLL(pcacline$WGS$runs$init$modelIII))
plot(hzar.mcmc.bindLL(pcacline$WGS$runs$init$modelIV))
plot(hzar.mcmc.bindLL(pcacline$WGS$runs$init$modelV))
plot(hzar.mcmc.bindLL(pcacline$WGS$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
pcacline$WGS$fitRs$chains <-
  lapply(pcacline$WGS$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
pcacline$WGS$fitRs$chains <-
  hzar.multiFitRequest(pcacline$WGS$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Just to be thorough, randomize the initial value for each fit
for(chain in c(1:18)){
  pcacline$WGS$fitRs$chains[[chain]]$modelParam$init["center"] <- runif(1,0,1809877)
  pcacline$WGS$fitRs$chains[[chain]]$modelParam$init["width"] <- runif(1,0,1809877)
}

#free freqs
for(chain in c(10:18)){
  pcacline$WGS$fitRs$chains[[chain]]$modelParam$init["pMin"] <- runif(1,0,1)
  pcacline$WGS$fitRs$chains[[chain]]$modelParam$init["pMax"] <- runif(1,0,1)
}

#mirror tails
for(chain in c(4:6,13:15)){
  pcacline$WGS$fitRs$chains[[chain]]$modelParam$init["deltaM"] <- runif(1,0,1809877)
  pcacline$WGS$fitRs$chains[[chain]]$modelParam$init["tauM"] <- runif(1,0,1)
}

#both tails
for(chain in c(7:9,16:18)){
  pcacline$WGS$fitRs$chains[[chain]]$modelParam$init["deltaL"] <- runif(1,0,1809877)
  pcacline$WGS$fitRs$chains[[chain]]$modelParam$init["deltaR"] <- runif(1,0,1809877)
  
  pcacline$WGS$fitRs$chains[[chain]]$modelParam$init["tauL"] <- runif(1,0,1)
  pcacline$WGS$fitRs$chains[[chain]]$modelParam$init["tauR"] <- runif(1,0,1)
  
}

pcacline$WGS$runs$chains <-  hzar.doChain.multi(pcacline$WGS$fitRs$chains,
                                                doPar=TRUE,
                                                inOrder=FALSE,
                                                count=3)

today<-format(Sys.Date(),format="%d%b%Y")
save(pcacline,file=paste("pcacline_WGS_",today,".rdata",sep=''))
load('pcacline_WGS_01Sep2021.rdata')

#did the runs converge?
summary(do.call(mcmc.list,
                lapply(pcacline$WGS$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )
summary(do.call(mcmc.list,
                lapply(pcacline$WGS$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )
summary(do.call(mcmc.list,
                lapply(pcacline$WGS$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )
summary(do.call(mcmc.list,
                lapply(pcacline$WGS$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )
summary(do.call(mcmc.list,
                lapply(pcacline$WGS$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )
summary(do.call(mcmc.list,
                lapply(pcacline$WGS$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Start aggregation of data for analysis
## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
pcacline$WGS$analysis$initDGs <- list(  nullModel =  hzar.dataGroup.null(pcacline$WGS$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
pcacline$WGS$analysis$initDGs$modelI <-  hzar.dataGroup.add(pcacline$WGS$runs$init$modelI)
pcacline$WGS$analysis$initDGs$modelII <- hzar.dataGroup.add(pcacline$WGS$runs$init$modelII)
pcacline$WGS$analysis$initDGs$modelIII <- hzar.dataGroup.add(pcacline$WGS$runs$init$modelIII)
pcacline$WGS$analysis$initDGs$modelIV <- hzar.dataGroup.add(pcacline$WGS$runs$init$modelIV)
pcacline$WGS$analysis$initDGs$modelV <- hzar.dataGroup.add(pcacline$WGS$runs$init$modelV)
pcacline$WGS$analysis$initDGs$modelVI <- hzar.dataGroup.add(pcacline$WGS$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created
pcacline$WGS$analysis$oDG <-  hzar.make.obsDataGroup(pcacline$WGS$analysis$initDGs)
pcacline$WGS$analysis$oDG <- hzar.copyModelLabels(pcacline$WGS$analysis$initDGs,
                                                  pcacline$WGS$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
pcacline$WGS$analysis$oDG <-  hzar.make.obsDataGroup(lapply(pcacline$WGS$runs$chains,
                                                            hzar.dataGroup.add),
                                                     pcacline$WGS$analysis$oDG);

## Do model selection based on the AICc scores
print(pcacline$WGS$analysis$AICcTable <- hzar.AICc.hzar.obsDataGroup(pcacline$WGS$analysis$oDG));

## Print out the model with the minimum AICc score
print(pcacline$WGS$analysis$model.name <-
        rownames(pcacline$WGS$analysis$AICcTable
        )[[ which.min(pcacline$WGS$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
pcacline$WGS$analysis$model.selected <- pcacline$WGS$analysis$oDG$data.groups[[pcacline$WGS$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(pcacline$WGS$analysis$model.selected,   names(pcacline$WGS$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(pcacline$WGS$analysis$model.selected))

((668744.7/1809877)*19.99)-96.86
((710441.5/1809877)*19.99)-96.86
((799218.8/1809877)*19.99)-96.86


## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(pcacline$WGS$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(pcacline$WGS$analysis$model.selected);


##
#xSeries <- seq(from = par("usr")[1], to = par("usr")[2], length.out = 109)
xSeries <- seq(from = min(pca_hzar_table$dist), to =max(pca_hzar_table$dist), length.out = 109)
if (par("xaxs") == "r") xSeries <- xSeries[2:108]

pca_fzCline = hzar.getCredParamRed(pca_cline_model_Data)
pca_fzCoor <- pca_fzCline$fzCline(xSeries)




####FIG 3 - CLINE ####
#save(cnstrc_cline,file="cnstrc_cline")
#save(pcacline,file="pcacline")
load("cnstrc_cline")
load("pcacline")

xSeries <- seq(from = par("usr")[1], to = par("usr")[2], length.out = 109)
xSeries <- seq(from = cnstrc_cline$WGS$obs$frame$dist[order(cnstrc_cline$WGS$obs$frame$dist)][1], to =cnstrc_cline$WGS$obs$frame$dist[order(cnstrc_cline$WGS$obs$frame$dist)][length(cnstrc_cline$WGS$obs$frame$dist)], length.out = 109)
if (par("xaxs") == "r") xSeries <- xSeries[2:108]

cnstrc_fzCline = hzar.getCredParamRed(cnstrc_cline$WGS$analysis$model.selected)
cnstrc_fzCoor <- cnstrc_fzCline$fzCline(xSeries)

pcacline_fzCline = hzar.getCredParamRed(pcacline$WGS$analysis$model.selected)
pcacline_fzCoor <- pcacline_fzCline$fzCline(xSeries)




conscnstrc_cline_res <- cbind.data.frame(cnstrc_cline$WGS$obs$frame,"Longitude"=pop_pos_gen$Longitude)
pcacline_res <- cbind.data.frame(pcacline$WGS$obs$frame,"Longitude"=pop_pos_gen$Longitude)

pdf(paste("fig3_cline.",today,".pdf",sep=''),width=5.5,height=4)

ggplot() + 
  geom_rect( mapping=aes(xmin=-91, xmax=-88.5, ymin=0, ymax=1), fill="lightgray",color=NA, alpha=0.3) +
  
  geom_point(data=conscnstrc_cline_res,aes(x = Longitude , y= obsFreq),color="purple") + 
  geom_line(aes(x = ((xSeries/1809877)*19.99)-96.86, y = cnstrc_cline$WGS$analysis$model.selected$ML.cline$clineFunc(xSeries)) ,color="purple")+ 
  geom_ribbon(data=cnstrc_fzCoor,aes(x = ((x/1809877)*19.99)-96.86, ymin = yMin, ymax = yMax),alpha=0.25,fill="purple") + 
  
  geom_point(data=pcacline_res,aes(x = Longitude , y= obsFreq),color="red") + 
  geom_line(aes(x = ((xSeries/1809877)*19.99)-96.86, y = pcacline$WGS$analysis$model.selected$ML.cline$clineFunc(xSeries)),color="red" )+ 
  geom_ribbon(data=pcacline_fzCoor,aes(x = ((x/1809877)*19.99)-96.86, ymin = yMin, ymax = yMax),alpha=0.25,fill="red") + 
  
  theme_bw(base_size = 8) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Longitude",y="") +   xlim(-98,-76)
dev.off()






Construct_mcmc_runs <- 
  rbind.data.frame(
    as.data.frame(
      do.call(mcmc.list,
              lapply(cnstrc_cline$WGS$runs$chains[10:12],
                     function(x) hzar.mcmc.bindLL(x[[3]]) ))[[1]][,c(1:2)]
    ),
    as.data.frame(
      do.call(mcmc.list,
              lapply(cnstrc_cline$WGS$runs$chains[10:12],
                     function(x) hzar.mcmc.bindLL(x[[3]]) ))[[2]][,c(1:2)]
    ),
    
    as.data.frame(
      do.call(mcmc.list,
              lapply(cnstrc_cline$WGS$runs$chains[10:12],
                     function(x) hzar.mcmc.bindLL(x[[3]]) ))[[3]][,c(1:2)]
    )
  )

pca_mcmc_runs <- 
  rbind.data.frame(
    as.data.frame(
      do.call(mcmc.list,
              lapply(pcacline$WGS$runs$chains[10:12],
                     function(x) hzar.mcmc.bindLL(x[[3]]) ))[[1]][,c(1:2)]
    ),
    as.data.frame(
      do.call(mcmc.list,
              lapply(pcacline$WGS$runs$chains[10:12],
                     function(x) hzar.mcmc.bindLL(x[[3]]) ))[[2]][,c(1:2)]
    ),
    
    as.data.frame(
      do.call(mcmc.list,
              lapply(pcacline$WGS$runs$chains[10:12],
                     function(x) hzar.mcmc.bindLL(x[[3]]) ))[[3]][,c(1:2)]
    )
  )

mcmc_runs <- 
  rbind.data.frame(
    cbind.data.frame(pca_mcmc_runs,"base"="PCA"),
    cbind.data.frame(Construct_mcmc_runs,"base"="CONSTRUCT")
  )

pdf(paste("figS7_clineMCMC.",today,".pdf",sep=''),width=5.5,height=4)

#mcmc_runs_plot <-  
ggplot(mcmc_runs,aes(x=((center/1809877)*19.99)-96.86,y=width/1000)) + 
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  facet_grid(base ~ .) + guides(fill=FALSE) +
  labs(x="Center (Longitude)",y="Width (km)",fill="Chain density") +
  scale_fill_viridis(option = "B",direction = -1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + theme_bw(base_size=18) +
  theme(strip.background =element_rect(fill="white")) 

dev.off()


####Fig 2####
poppc1_long_plot <- 
  ggplot(pop_PCAcoord,aes(x=Longitude, y=pc1_norm)) + 
  # ggplot(pop_PCAcoord,aes(x=Dim.1, y=Dim.2,color=Type)) + 
  geom_rect( mapping=aes(xmin=-91, xmax=-88.5, ymin=0, ymax=1), fill="lightgray",color=NA, alpha=0.5) +
  
  geom_point() + 
  #   geom_text(aes(label=Pop)) +
  guides(color=FALSE) +
  labs(x="",y="PC1 Normalized") +
  theme_bw(base_size = 8) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_colour_gradient2(low = "darkblue", high = "darkgreen", mid = "lightgray", midpoint = -90, na.value = NA) +
  xlim(-98,-76)

construct_pie_plot <- 
  ggplot() + 
  geom_rect( mapping=aes(xmin=-91, xmax=-88.5, ymin=30, ymax=36), color=NA, alpha=0.3) +
  
  geom_path(data = countries, aes(x = long, y = lat, group = group),color="lightgray") +
  geom_path(data = rivers, aes(x = long, y = lat, group = group),color="lightblue") +
  #   xlim(-94.5,-84.5)+ ylim(30,34)+ #all samples
  xlim(-98,-76)+ ylim(30,36)+ #all samples
  
  geom_scatterpie(data= admix.props.df, cols=c("K1","K2","K3") , aes(x=Longitude, y=Latitude, r = 0.3),color=NA) +
  
  theme_bw(base_size = 8) + guides(fill = FALSE) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="",y="Latitude",fill="") +
  scale_fill_manual(values=c( "darkblue","darkgreen","gray")) #+coord_sf() 


twoCline_plot <- 

ggplot() + 
  geom_rect( mapping=aes(xmin=-91, xmax=-88.5, ymin=0, ymax=1), fill="lightgray",color=NA, alpha=0.3) +
  
  geom_point(data=conscnstrc_cline_res,aes(x = Longitude , y= obsFreq),color="purple") + 
  geom_line(aes(x = ((xSeries/1809877)*19.99)-96.86, y = cnstrc_cline$WGS$analysis$model.selected$ML.cline$clineFunc(xSeries)) ,color="purple")+ 
  geom_ribbon(data=cnstrc_fzCoor,aes(x = ((x/1809877)*19.99)-96.86, ymin = yMin, ymax = yMax),alpha=0.25,fill="purple") + 
  
  geom_point(data=pcacline_res,aes(x = Longitude , y= obsFreq),color="red") + 
  geom_line(aes(x = ((xSeries/1809877)*19.99)-96.86, y = pcacline$WGS$analysis$model.selected$ML.cline$clineFunc(xSeries)),color="red" )+ 
  geom_ribbon(data=pcacline_fzCoor,aes(x = ((x/1809877)*19.99)-96.86, ymin = yMin, ymax = yMax),alpha=0.25,fill="red") + 
  
  theme_bw(base_size = 8) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Longitude",y="") +   xlim(-98,-76)


pdf(paste("fig2_cline.",today,".pdf",sep=''),width=5,height=4)


plot_grid(  poppc1_long_plot,construct_pie_plot, twoCline_plot,  labels = c('A', 'B','C'), label_size = 12,ncol = 1, align = 'vh',axis='tbrl')

dev.off()

###Fig 2 inset##
pdf(paste("fig2_inset.",today,".pdf",sep=''),width=1.5,height=1.5)

ggplot() + 
  geom_path(data = countries, aes(x = long, y = lat, group = group),color="lightgray",size=0.2) +
  geom_rect( mapping=aes(xmin=-98, xmax=-76, ymin=30, ymax=36), fill=NA,color="black", alpha=0.5,size=0.2) +
  
  coord_sf() + 
  
  xlim(-130,-70)+ ylim(10,60)+ #all samples
  
  theme_bw(base_size = 8) + guides(color = FALSE) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Longitude",y="Latitude") 

dev.off()


####FST####

fst <- fread('fst.feb1.txt')
fst <- fread('fst.weighted.txt')

fst <- fst[fst$V2 != fst$V3,]

fst <- fst[!duplicated(data.frame(t(apply(fst[,c(2,3)],1,sort)))),]
fst$V1[fst$V1 < 0 ] <- 0


fst_latlon_tmp <- left_join(fst,pop_pos_gen[,c(1,2,11)],by=c("V2"="mixName"))
names(fst_latlon_tmp)[c(4,5)] <- c("pop1_lat","pop1_long")

fst_latlon <- left_join(fst_latlon_tmp,pop_pos_gen[,c(1,2,11)],by=c("V3"="mixName"))
names(fst_latlon)[c(6,7)] <- c("pop2_lat","pop2_long")


fst_latlon$dist <-  diag(distm(cbind(fst_latlon$pop1_long,fst_latlon$pop1_lat),cbind(fst_latlon$pop2_long,fst_latlon$pop2_lat)))

fst_latlon$fstfst <- fst_latlon$V1/(1-fst_latlon$V1)

Xyy_pops <- unique(indPopInfo$statePop[indPopInfo$Collec == "Old" & indPopInfo$Type == "XYY"])
Xy_pops <- unique(indPopInfo$statePop[indPopInfo$Collec == "Old" & indPopInfo$Type == "XY"])
Xyy_hybs <- unique(indPopInfo$pop[indPopInfo$Collec == "New" & indPopInfo$riverside == "E"])
Xy_hybs <- unique(indPopInfo$pop[indPopInfo$Collec == "New" & indPopInfo$riverside == "W"])


fst_latlon$comp[fst_latlon$V2 %in% Xyy_pops & fst_latlon$V3 %in% Xyy_pops] <- "Xyy_pops"
fst_latlon$comp[fst_latlon$V2 %in% Xy_pops & fst_latlon$V3 %in% Xy_pops] <- "Xy_pops"
fst_latlon$comp[fst_latlon$V2 %in% Xyy_hybs & fst_latlon$V3 %in% Xyy_hybs] <- "Xyy_hybs"
fst_latlon$comp[fst_latlon$V2 %in% Xy_hybs & fst_latlon$V3 %in% Xy_hybs] <- "Xy_hybs"

ggplot(fst_latlon[complete.cases(fst_latlon),],aes(x=dist/1000,y=fstfst,color=comp)) + 
  geom_point() + geom_smooth(method = 'lm')+
  labs(x="Distance (Km)",y="Fst/1-Fst") +
  theme_bw(base_size = 18)

lm(formula =  fstfst ~ dist,data = fst_latlon[fst_latlon$comp == "Xyy_pops",])
#lm(formula =  dist ~ fstfst  ,data = fst_latlon[fst_latlon$comp == "Xyy_pops",])

#calculating selection

sqrt(8*((1/1.068e-07)/(4*6))^2) / 283710  #6 is approx density from Pickup 2013


####rec. rate####
#import linkage map
CMMB <- fread('CHRR_linkage_map.tsv')
CMMB$LG[CMMB$LG == "L.10"] <- "X"
CMMB$LG[CMMB$LG == "L.7"] <- "A1"
CMMB$LG[CMMB$LG == "L.8"] <- "A2"
CMMB$LG[CMMB$LG == "L.5"] <- "A3"
CMMB$LG[CMMB$LG == "L.3"] <- "A4"

CMMB <- rbind.data.frame(
  CMMB[CMMB$LG != "A3",],
  CMMB %>%  filter(LG=="A3") %>% mutate(LG_BP = maxA3pos - LG_BP)
)

minCm <- min(CMMB$cM[CMMB$LG == "A3"])
CMMB$cM[CMMB$LG == "A3" ] <- CMMB$cM[CMMB$LG == "A3" ] - minCm

maxCM <- max(CMMB$cM[CMMB$LG == "A3"])
CMMB$cM[CMMB$LG == "A3" ] <- maxCM - CMMB$cM[CMMB$LG == "A3" ] 

CMMB$mb <- CMMB$LG_BP / 1000000
CMMB$chrom = factor(CMMB$LG, levels=c('A1','A2','A4','X','A3'))

#marey_plot <- 
  ggplot(CMMB,aes(x=mb,y=cM)) + geom_point()  +  facet_grid(. ~ chrom,scales = "free") +theme_bw()

#fit linkage map models
CMMB_reform <- CMMB[,c(1:3,9)]
CMMB_reform <- cbind.data.frame(CMMB_reform,"sp"="Rumhas")
names(CMMB_reform) <- c("chr", "marker", "cm", "mb", "sp")
CMMB_clean <- cleanup_maps(CMMB_reform)
CMMB_mmRem <- remove_mismapped(CMMB_clean) # remove mismapped markers from all_maps_clean
rec_piecewise <-  est_rec_piece_func(CMMB_mmRem[, c("chr", "cm", "mb")])

#recombination estimation for SNPs
loc_comp <- loc[complete.cases(loc) & loc$LG != "U", ]
loc_comp$LG <- paste("chr",loc_comp$LG,sep="")

loop = 0
window.size = wSize = 1 
for(z in unique(loc_comp$LG)){
  loc_scaf_loop <- loc_comp[loc_comp$LG == z,]
  loc_scaf_loop <- loc_scaf_loop[order(lPos),]
  
  mb<-as.numeric(loc_scaf_loop$lPos)/1000000
  Y.piece.data = data.frame(mb=mb-(window.size/2000))
  Y.piece.data = rbind(Y.piece.data, tail(mb, n=1)+(window.size/2000))
  
  #predict
  piece.fun=rec_piecewise[[z]][["func"]]
  Y.piece.raw = predict(piece.fun, newdata=Y.piece.data)
  
  #make table
  Y.piece.raw[Y.piece.raw < 0] = 0 #no negative cM
  
  if(loop == 0){
    loc_scaf_cM <-   cbind.data.frame(  loc_scaf_loop, "cM" = Y.piece.raw[-c(1)])
    
  }else{
    loc_scaf_cM <- rbind.data.frame(loc_scaf_cM,
                                    cbind.data.frame(  loc_scaf_loop, "cM" = Y.piece.raw[-c(1)])
    )
  }
  loop =+ 1
}

loc_scaf_cM$chrom = factor(loc_scaf_cM$LG, levels=c('chrA1','chrA2','chrA4','chrX','chrA3'))

#mareyFit_plot <-
ggplot(loc_scaf_cM,aes(x=mb,y=cM)) + geom_point()  +  facet_grid(. ~ chrom,scales = "free") + theme_bw()

  windowedrec <- cbind.data.frame("LG"=NA,"winNum"=NA,"mb_s"=NA,"mb_e"=NA,"cmmb"=NA)

  for(z in unique(loc_scaf_cM$chrom)){
    loc_scaf_cM_loop <- loc_scaf_cM[loc_scaf_cM$LG == z,]
    mb <- seq(from=0,to=max(loc_scaf_cM_loop$mb),0.001)
  
    Y.piece.data = data.frame(mb=mb-(window.size/2000))
    Y.piece.data = rbind(Y.piece.data, tail(mb, n=1)+(window.size/2000))
  
    #predict
    piece.fun=rec_piecewise[[z]][["func"]]
    Y.piece.raw = predict(piece.fun, newdata=Y.piece.data)
  
    #make table
    Y.piece.raw[Y.piece.raw < 0] = 0 #no negative cM
    loc_scaf_cM_loop <-   cbind.data.frame(  mb, "cM" = Y.piece.raw[-c(1)])
  
    win = 0
    wl = 0
  
    while(wl < max(loc_scaf_cM_loop$mb)){
      wh = wl + wSize
    
      loc_scaf_cM_loop_w <- loc_scaf_cM_loop[loc_scaf_cM_loop$mb > wl & loc_scaf_cM_loop$mb <= wh,]

      if(all(is.na(loc_scaf_cM_loop_w$mb)) || all(is.na(loc_scaf_cM_loop_w$cM))){
        localRec <- NA
      
      } else{
        localRec <- lm(loc_scaf_cM_loop_w$cM ~ loc_scaf_cM_loop_w$mb  ) 
        localRec <- localRec$coefficients[2]
        if(localRec < 0 || is.na(localRec)){
          localRec <- NA
        }
      }
    
      windowedrec <- rbind.data.frame(windowedrec, 
                                    cbind.data.frame("LG"=z,"winNum"=win,"mb_s"=wl,"mb_e"=wh,"cmmb"=localRec)
      )
      win = win + 1
      wl = wl + wSize
      cat(z,wl,"\n")
    }
  
  }



#windowedreC <- windowedrec[complete.cases(windowedrec), ]
  windowedreC <- windowedrec
  
  windowedreC$LG <- gsub("(.{3})", "\\1 ", windowedreC$LG)
windowedreC <- separate(windowedreC, LG, c("chr","LG"), 
                 sep = " ", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")

windowedreC$chrom = factor(windowedreC$LG, levels=c('A1','A2','A4','X','A3'))
windowedreC$mb <- as.numeric(windowedreC$mb_s)
#windowedreC$mb <- round(windowedreC$mb,digits = 2)





####strcture cline (HZAR) by chrom####
constrcPop <- fread("GBS.pop.num")

windowedreC$recBin <- "low"
windowedreC$recBin[windowedreC$cmmb > 0.1] <- "high"

loc$mb_r <- round(loc$mb)

locHI_rec <- left_join(loc,windowedreC[,c(2,4,9)],by=c("LG"="LG","mb_r"="mb_s"))
locHI_rec <- locHI_rec[!is.na(locHI_rec$recBin),]

s012m <- as.matrix(s012)

clines_recLG <-foreach(chrom=unique(locHI_rec$LG),.combine=rbind) %do% {
  

  print(c("LG=",chrom) )
  local_rows  <- locHI_rec[locHI_rec$chrom %in% chrom ,]
  
  clines_rec_tmp <- 
  foreach(recBin=unique(local_rows$recBin),.combine=rbind) %do% {
  
    print(c("recBin=",recBin) )
    local_cols  <- local_rows$SNP[local_rows$recBin ==  recBin]
    
  s012m_sub <-  as.data.frame( s012m[,c(local_cols)])
  
  pop.data.matrix <- matrix(NA,nrow=45,ncol=ncol(s012m_sub))
  for(i in 1:nrow(pop.data.matrix)){
    pop.data.matrix[i,] <- colMeans(
      s012m_sub[
        which(constrcPop$V2==i),, ##
        drop=FALSE
        ],na.rm=TRUE
    )/2
  }


  pop.data.matrix[is.nan(pop.data.matrix)] <- NA

  my.run <- conStruct(spatial = FALSE, #non-spatial
                      K = as.numeric(2), 
                      freqs = pop.data.matrix,
                      geoDist = as.matrix(fields::rdist.earth(x1 =  pop_pos_gen[order(mixName),c(2,1)])), 
                      coords = as.matrix(pop_pos_gen[order(mixName),c(2,1)]),
                      prefix = paste("nsp_k2",chrom,recBin,sep="_"),
                      make.figs=FALSE)
  
 
  admix.props <- as.data.frame(my.run$chain_1$MAP$admix.proportions)
  names(admix.props) <- c("K1","K2")
  admix.props.df <- cbind.data.frame( admix.props,pop_pos_gen[order(pop_pos_gen$mixName),])
  
  ##start HZAR
  cnstrc_dist <-  diag(distm(cbind(min(admix.props.df$Longitude),admix.props.df$Latitude),cbind(admix.props.df$Longitude,admix.props.df$Latitude)))
  
  cnstrc_hzar_table <- cbind.data.frame("Locality_ID"=admix.props.df$Pop,"Locality"=admix.props.df$State,"dist"=cnstrc_dist,"strc1"=admix.props.df$K2,"strc2"=admix.props.df$K1,"Nsamples"=admix.props.df$n)
  
  cnstrc_cline <- hzar.doMolecularData1DPops(cnstrc_hzar_table$dist, cnstrc_hzar_table$strc1, cnstrc_hzar_table$Nsamples)

  
  ##start fitting
  chainLength=1e5;                       
  
  ## Make each model run off a separate seed
  mainSeed=
    list(A=c(596,528,124,978,544,99),
         B=c(528,124,978,544,99,596),
         C=c(124,978,544,99,596,528))
  
  ## Blank out space in memory to hold molecular analysis
  if(length(apropos("^cnstrc_cline$",ignore.case=FALSE)) == 0 ||
     !is.list(cnstrc_cline) ) cnstrc_cline <- list()
  ## We are doing just the one allele at one locus, but it is
  ## good to stay organized.
  cnstrc_cline$WGS <- list();
  ## Space to hold the observed data
  cnstrc_cline$WGS$obs <- list();
  cnstrc_cline$WGS$obs <- hzar.doMolecularData1DPops(cnstrc_hzar_table$dist, cnstrc_hzar_table$strc1, cnstrc_hzar_table$Nsamples)
  
  ## Space to hold the models to fit
  cnstrc_cline$WGS$models <- list();
  ## Space to hold the compiled fit requests
  cnstrc_cline$WGS$fitRs <- list();
  ## Space to hold the output data chains
  cnstrc_cline$WGS$runs <- list();
  ## Space to hold the analysed data
  cnstrc_cline$WGS$analysis <- list();
  
  cnstrc_cline.loadWGSmodel <- function(scaling,tails,
                                        id=paste(scaling,tails,sep="."))
    cnstrc_cline$WGS$models[[id]] <<- hzar.makeCline1DFreq(cnstrc_cline$WGS$obs, scaling, tails)
  
  cnstrc_cline.loadWGSmodel("fixed","none","modelI");
  cnstrc_cline.loadWGSmodel("fixed" ,"mirror","modelII");
  cnstrc_cline.loadWGSmodel("fixed","both","modelIII");
  
  cnstrc_cline.loadWGSmodel("free" ,"none","modelIV");
  cnstrc_cline.loadWGSmodel("free" ,"mirror","modelV");
 cnstrc_cline.loadWGSmodel("free" ,"both","modelVI");
  
  #print(cnstrc_cline$WGS$models)
  
  ## Compile each of the models to prepare for fitting
  cnstrc_cline$WGS$fitRs$init <- sapply(cnstrc_cline$WGS$models,
                                        hzar.first.fitRequest.old.ML,
                                        obsData=cnstrc_cline$WGS$obs,
                                        verbose=FALSE,
                                        simplify=FALSE)
  
  ## Run just one of the models for an initial chain
  cnstrc_cline$WGS$runs$init <- list()
  cnstrc_cline$WGS$runs$init$modelI <-  hzar.doFit(cnstrc_cline$WGS$fitRs$init$modelI)
  cnstrc_cline$WGS$runs$init$modelII <-  hzar.doFit(cnstrc_cline$WGS$fitRs$init$modelII)
  cnstrc_cline$WGS$runs$init$modelIII <-  hzar.doFit(cnstrc_cline$WGS$fitRs$init$modelIII)
  cnstrc_cline$WGS$runs$init$modelIV <-  hzar.doFit(cnstrc_cline$WGS$fitRs$init$modelIV)
  cnstrc_cline$WGS$runs$init$modelV <-  hzar.doFit(cnstrc_cline$WGS$fitRs$init$modelV)
  cnstrc_cline$WGS$runs$init$modelVI <-  hzar.doFit(cnstrc_cline$WGS$fitRs$init$modelVI)
  
  ## Compile a new set of fit requests using the initial chains 
  cnstrc_cline$WGS$fitRs$chains <-
    lapply(cnstrc_cline$WGS$runs$init,
           hzar.next.fitRequest)
  
  ## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
  cnstrc_cline$WGS$fitRs$chains <-
    hzar.multiFitRequest(cnstrc_cline$WGS$fitRs$chains,
                         each=3,
                         baseSeed=NULL)
  
  ## Just to be thorough, randomize the initial value for each fit
  #for(chain in c(1:9)){
  for(chain in c(1:18)){
    cnstrc_cline$WGS$fitRs$chains[[chain]]$modelParam$init["center"] <- runif(1,0,1809877)
    cnstrc_cline$WGS$fitRs$chains[[chain]]$modelParam$init["width"] <- runif(1,0,1809877)
  }
  
  #free freqs
  for(chain in c(10:18)){
    cnstrc_cline$WGS$fitRs$chains[[chain]]$modelParam$init["pMin"] <- runif(1,0,1)
    cnstrc_cline$WGS$fitRs$chains[[chain]]$modelParam$init["pMax"] <- runif(1,0,1)
  }
  
  #mirror tails
  for(chain in c(4:6,13:15)){
  #for(chain in c(4:6)){
    cnstrc_cline$WGS$fitRs$chains[[chain]]$modelParam$init["deltaM"] <- runif(1,0,1809877)
    cnstrc_cline$WGS$fitRs$chains[[chain]]$modelParam$init["tauM"] <- runif(1,0,1)
  }
  
  #both tails
  for(chain in c(7:9,16:18)){
  #for(chain in c(7:9)){
    cnstrc_cline$WGS$fitRs$chains[[chain]]$modelParam$init["deltaL"] <- runif(1,0,1809877)
    cnstrc_cline$WGS$fitRs$chains[[chain]]$modelParam$init["deltaR"] <- runif(1,0,1809877)
    
    cnstrc_cline$WGS$fitRs$chains[[chain]]$modelParam$init["tauL"] <- runif(1,0,1)
    cnstrc_cline$WGS$fitRs$chains[[chain]]$modelParam$init["tauR"] <- runif(1,0,1)
    
  }
  
  cnstrc_cline$WGS$runs$chains <-  hzar.doChain.multi(cnstrc_cline$WGS$fitRs$chains,
                                                      doPar=TRUE,
                                                      inOrder=FALSE,
                                                      count=3)
  

  
  ## Start aggregation of data for analysis
  ## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
  cnstrc_cline$WGS$analysis$initDGs <- list(  nullModel =  hzar.dataGroup.null(cnstrc_cline$WGS$obs))
  
  ## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
  cnstrc_cline$WGS$analysis$initDGs$modelI <-  hzar.dataGroup.add(cnstrc_cline$WGS$runs$init$modelI)
  cnstrc_cline$WGS$analysis$initDGs$modelII <- hzar.dataGroup.add(cnstrc_cline$WGS$runs$init$modelII)
  cnstrc_cline$WGS$analysis$initDGs$modelIII <- hzar.dataGroup.add(cnstrc_cline$WGS$runs$init$modelIII)
  cnstrc_cline$WGS$analysis$initDGs$modelIV <- hzar.dataGroup.add(cnstrc_cline$WGS$runs$init$modelIV)
  cnstrc_cline$WGS$analysis$initDGs$modelV <- hzar.dataGroup.add(cnstrc_cline$WGS$runs$init$modelV)
  cnstrc_cline$WGS$analysis$initDGs$modelVI <- hzar.dataGroup.add(cnstrc_cline$WGS$runs$init$modelVI)
  
  ## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created
  cnstrc_cline$WGS$analysis$oDG <-  hzar.make.obsDataGroup(cnstrc_cline$WGS$analysis$initDGs)
  cnstrc_cline$WGS$analysis$oDG <- hzar.copyModelLabels(cnstrc_cline$WGS$analysis$initDGs,
                                                        cnstrc_cline$WGS$analysis$oDG)
  
  ## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
  cnstrc_cline$WGS$analysis$oDG <-  hzar.make.obsDataGroup(lapply(cnstrc_cline$WGS$runs$chains,
                                                                  hzar.dataGroup.add),
                                                           cnstrc_cline$WGS$analysis$oDG);
  
  ## Do model selection based on the AICc scores
  cnstrc_cline$WGS$analysis$AICcTable <- hzar.AICc.hzar.obsDataGroup(cnstrc_cline$WGS$analysis$oDG)
  
  ## Print out the model with the minimum AICc score
  cnstrc_cline$WGS$analysis$model.name <-
    rownames(cnstrc_cline$WGS$analysis$AICcTable
    )[[ which.min(cnstrc_cline$WGS$analysis$AICcTable$AICc )]]
  
  ## Extract the hzar.dataGroup object for the selected model
  cnstrc_cline$WGS$analysis$model.selected <- cnstrc_cline$WGS$analysis$oDG$data.groups[[cnstrc_cline$WGS$analysis$model.name]]
  
  ## Look at the variation in parameters for the selected model
  
  ## Print the maximum likelihood cline for the selected model

  
    
  cline_tmp <- unlist(  c(
  "LG"=chrom,
  "recBin" = recBin,
  "bestModel"=  cnstrc_cline$WGS$analysis$model.name <-
              rownames(cnstrc_cline$WGS$analysis$AICcTable
              )[[ which.min(cnstrc_cline$WGS$analysis$AICcTable$AICc )]],
    "center"=hzar.get.ML.cline(cnstrc_cline$WGS$analysis$model.selected)$param.all$center,
    "width"=hzar.get.ML.cline(cnstrc_cline$WGS$analysis$model.selected)$param.all$width,
    
    hzar.getLLCutParam(cnstrc_cline$WGS$analysis$model.selected,   names(cnstrc_cline$WGS$analysis$model.selected$data.param))[1:4]
    ))
  
  cline_tmp


  }
  clines_rec_tmp
}


write.table( clines_recLG, file = "clines_recLG.txt", append = FALSE, quote = FALSE, sep = "\t", 
             eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
             col.names = TRUE, qmethod = c("escape", "double"))

clines_recLG_raw <- fread("clines_recLG.txt")

clines_recLG_raw$lon <- ((clines_recLG_raw$center  /1809877)*19.99)-96.86
clines_recLG_raw$lonLow <- ((clines_recLG_raw$center2LLLow  /1809877)*19.99)-96.86
clines_recLG_raw$lonHigh <- ((clines_recLG_raw$center2LLHigh  /1809877)*19.99)-96.86

clines_recLG_center <- clines_recLG_raw[,c(1,2,10,12,11)]
clines_recLG_center$var <- "center"
names(clines_recLG_center)[3:5] <- c("est","low","high")
clines_recLG_width <- clines_recLG_raw[,c(1,2,5,8:9)]
clines_recLG_width$var <- "width"
names(clines_recLG_width)[3:5] <- c("est","low","high")

clines_recLG <- rbind.data.frame(clines_recLG_center,clines_recLG_width)


ggplot(clines_recLG,aes(x=LG,y=est,color=recBin)) +
  geom_point(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin=low,ymax=high,width=0.2),position = position_dodge(width = 0.9)) +
  facet_grid(var~.,scales = "free") + theme_bw(base_size=18) + labs(x="",y="")

####rec rate and PC1####
#hybPCA <- PCA(s012[hybPCAcoord$flower == "F"], graph = FALSE) #missMDA assumes continuous variables
#hybPCA <- PCA(s012, graph = FALSE) #missMDA assumes continuous variables

hybPCAcoord <- cbind.data.frame(hybPCA$ind$coord,indsInfo[,-1])
#hybPCAcoord <- cbind.data.frame(hybPCA$ind$coord,indsInfo[hybPCAcoord$flower == "F",-1])
#ggplot(hybPCAcoord,aes(x=Dim.1,y=Dim.2,color=riverside)) + geom_point() +theme_bw()

hybPCAloc <- cbind.data.frame(loc,hybPCA$var$contrib)

rec_PC1_raw <- hybPCAloc[hybPCAloc$LG != "U",]

windowed_PC1 <- cbind.data.frame("LG"=NA,"winNum"=NA,"mb_s"=NA,"mb_e"=NA,"Sdim1"=NA)


for(z in unique(rec_PC1_raw$chrom)){
  rec_PC1_raw_loop <- rec_PC1_raw[rec_PC1_raw$LG == z,]
  mb <- seq(from=0,to=max(rec_PC1_raw_loop$mb),0.001)
  
  win = 0
  wl = 0
  
  while(wl < max(rec_PC1_raw_loop$mb)){
    wh = wl + wSize
    
    rec_PC1_raw_loop_w <- rec_PC1_raw_loop[rec_PC1_raw_loop$mb > wl & rec_PC1_raw_loop$mb <= wh,]
    
    if(all(is.na(rec_PC1_raw_loop_w$mb)) || all(is.na(rec_PC1_raw_loop_w$Dim.1))){
      localPC1 <- NA
      
    } else{
      localPC1 <- sum(rec_PC1_raw_loop_w$Dim.1) / length(rec_PC1_raw_loop_w$Dim.1)
    }
    
    windowed_PC1 <- rbind.data.frame(windowed_PC1, 
                                     cbind.data.frame("LG"=z,"winNum"=win,"mb_s"=wl,"mb_e"=wh,"Sdim1"=localPC1)
    )
    win = win + 1
    wl = wl + wSize
    cat(z,wl,"\n")
  }
}

rec_PC1  <- left_join(windowedreC,  windowed_PC1,by=c("chrom"="LG","mb"="mb_s"))
rec_PC1_comp <- rec_PC1[complete.cases(rec_PC1),]
rec_PC1_comp$Chrom = factor(rec_PC1_comp$chrom, levels=c('A1','A2','A4','X','A3'))

# rec_PC1_plot <- 
ggplot(rec_PC1_comp,aes(x=cmmb,y=Sdim1)) + theme_bw(base_size=15) +
  geom_point(alpha=0.2) + stat_smooth(method="lm",color="black") + facet_grid(.~Chrom,scales = "free") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Recombination rate (cM/mb)",y="Mean PC1 Effect Size") +
  theme(strip.background =element_rect(fill="white")) 

PC1_rec_lm <-   lm(data=rec_PC1_comp,formula = Sdim1 ~ cmmb*chrom )
PC1_rec_lmsm <- summary(PC1_rec_lm)  

write.table(PC1_rec_lmsm$coefficients, file = "pc1_rec_lm", append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = TRUE,    col.names = TRUE, qmethod = c("escape", "double"))

ggplot(rec_PC1_comp,aes(x=mb_s,y=Sdim1)) + theme_bw(base_size=15) +
  geom_point(alpha=0.2) + 
  #stat_smooth(method="lm",color="black") + 
  facet_grid(.~Chrom,scales = "free") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Recombination rate (cM/mb)",y="Mean PC1 Effect Size") +
  theme(strip.background =element_rect(fill="white")) 

rec_PC1_comp$recBin <- "high"
rec_PC1_comp$recBin[rec_PC1_comp$cmmb < 0.1] <- "low"
ggplot(rec_PC1_comp,aes(x=mb,y=cmmb)) + facet_grid(. ~ LG,scale="free",space="free") + geom_point(aes(color=recBin)) +theme_bw()

#####bgc####
#hybrid index
BGCHI <- cbind.data.frame( fread('GBS.mis60.maf01.noY.bgc40k.hi.out'),
                          "inds"= fread('GBS.admx.inds',header=FALSE))

strc_Pop_k1 <- strc_Pop[strc_Pop$variable == "k1" & strc_Pop$k == "2",]

strcBGC <- left_join(BGCHI,strc_Pop_k1,by=c("inds.V1"="inds"))

ggplot(strcBGC,aes(x=median,y=value)) + 
  geom_point(size=3,aes(color=Longitude)) +
  # geom_text(aes(label=momID)) +
  labs(x="bgc HI",y="STRUCTURE HI") + geom_abline(slope=1) +
  theme_bw()  + theme_bw(base_size = 18) + xlim(0,1) + ylim(0,1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_colour_gradient2(low = "darkblue", high = "darkgreen", mid = "lightgray", midpoint = -91, na.value = NA) 


#strbhc_lm <-  lm(formula = value ~ median * Longitude, data = strcBGC) 
#summary(strbhc_lm)
cor.test(strcBGC$value, strcBGC$median, method=c("pearson"))

#cline analysis
#
loc_Yin <- fread("GBS.mis60.maf01.bgc.012.chrom.pos")
names(loc_Yin) <- c("scaffold","sPos","LG","lPos")

loc_Yin$LG[loc_Yin$LG == "L.10"] <- "X"
loc_Yin$LG[loc_Yin$LG == "L.7"] <- "A1"
loc_Yin$LG[loc_Yin$LG == "L.8"] <- "A2"
loc_Yin$LG[loc_Yin$LG == "L.5"] <- "A3"
loc_Yin$LG[loc_Yin$LG == "L.3"] <- "A4"
loc_Yin$LG[is.na(loc_Yin$LG)] <- "U"

maxA3pos <- 175014771
loc_Yin$lPos[loc_Yin$LG == "A3"] <- maxA3pos - loc_Yin$lPos[loc_Yin$LG == "A3"] 

loc_Yin$mb <- loc_Yin$lPos / 1000000
loc_Yin$SNP <- seq(1,length(loc_Yin$lPos),1)
loc_Yin$chrom = factor(loc_Yin$LG, levels=c('A1','A2','A4','X','A3','U'))

alpha <-   fread('GBS.mis60.maf01.bgc40k.alpha.out')[,-1] 
#alpha <-   fread('GBS.mis60.maf01.noY.bgc40k.alpha.out')[,-1] 
names(alpha) <- c("var_mean","var_median","var_CI_LB","var_CI_UB")

gamma <- fread('GBS.mis60.maf01.bgc40k.gamma-quantile.out')[,-1] 
#gamma <- fread('GBS.mis60.maf01.noY.bgc40k.gamma-quantile.out')[,-1] 
names(gamma) <- c("Q_mean","Q_median","Q_LB","Q_UB")
alpha_bound <- cbind.data.frame(loc_Yin,alpha,gamma,"var"="alpha")

beta <-   fread('GBS.mis60.maf01.bgc40k.beta.out')[,-1] 
#beta <-   fread('GBS.mis60.maf01.noY.bgc40k.beta.out')[,-1] 
names(beta) <- c("var_mean","var_median","var_CI_LB","var_CI_UB")

zeta <- fread('GBS.mis60.maf01.bgc40k.zeta-quantile.out')[,-1] 
#zeta <- fread('GBS.mis60.maf01.noY.bgc40k.zeta-quantile.out')[,-1] 
names(zeta) <- c("Q_mean","Q_median","Q_LB","Q_UB")
zeta_bound <- cbind.data.frame(loc_Yin,beta,zeta,"var"="beta")

bgc_bound <- rbind.data.frame(alpha_bound,zeta_bound)
bgc_bound$chrom = factor(bgc_bound$LG, levels=c('A1','A2','A4','X','A3','U'))

bgc_bound$excess[ bgc_bound$var_CI_UB < 0] <- T
bgc_bound$excess[bgc_bound$var_CI_LB > 0 ] <- T

bgc_bound$outlier[bgc_bound$Q_median < 0.05 ] <- T
bgc_bound$outlier[ bgc_bound$Q_median > 0.95] <- T

ggplot(data=bgc_bound,aes(x=var_median)) + 
  geom_vline(xintercept=0, color = "gray", size=1,alpha=0.25) +  
  
  geom_histogram()+
  facet_grid( . ~ var,scales = "free",space="free") + 
  theme_bw()  + theme_bw(base_size = 18) + 
  labs(x="",y="Counts") +
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values=c( 'black', 'gray'))


ggplot(data=bgc_bound,aes(x=var_median)) + 
  geom_vline(xintercept=0, color = "gray", size=1,alpha=0.25) +  
  
  geom_histogram()+
  facet_grid( chrom ~ var,scales = "free",space="free") + 
  theme_bw()  + theme_bw(base_size = 18) + 
  labs(x="",y="Counts") +
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values=c( 'black', 'gray'))

#library(e1071)
#skewness(bgc_bound$var_median[bgc_bound$var == "alpha"])  
#t.test(bgc_bound$var_median[bgc_bound$var == "alpha"])
#skewness(bgc_bound$var_median[bgc_bound$var == "beta"])  

#library(mclust) #map function conflicts with tidyverse^
#density_model <- densityMclust(bgc_bound$var_median[bgc_bound$var == "alpha"])
#summary(density_model)
#m_model <- Mclust(bgc_bound$var_median[bgc_bound$var == "alpha"])
#summary(m_model)

#par(mfrow = c(2,1))
#plot(density_model, what = "density", data = bgc_bound$var_median[bgc_bound$var == "alpha"], breaks = 15)
#plot(m_model, what="classification")
#par(mfrow = c(1,1))
#detach("package:mclust", unload=TRUE)


#bgc_bound <- bgc_bound[bgc_bound$chrom != "U",]

#bgc_plot <- 
  ggplot() + 
    geom_line(data=bgc_bound,color="black",alpha=0.1,aes(x=mb,y=var_median)) +
    
  geom_point(data=bgc_bound,alpha=0.1,color="white",aes(x=mb,y=var_median)) +
  geom_point(data=bgc_bound[!is.na(bgc_bound$excess),],shape=4,aes(x=mb,y=var_median),color="gray") +
  geom_point(data=bgc_bound[!is.na(bgc_bound$outlier),],color="black",shape=8,aes(x=mb,y=var_median)) +
    
  geom_hline(yintercept=0, color = "gray", size=1,alpha=0.25) +  
  facet_grid( var ~ chrom,scales = "free",space="free_x") + 
  theme_bw()  + theme_bw(base_size = 18) + 
  labs(x="Position (mb)",y=" ") +
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values=c( 'black', 'gray'))
  

  bgc_bound_beta <- bgc_bound[bgc_bound$var == "beta" & bgc_bound$LG != "U",]
 # bgc_bound_beta$linkage <- "A"
  bgc_bound_beta$LG[bgc_bound_beta$chrom == "A3" & bgc_bound_beta$mb > 80 ] <- "A3PAR"
  bgc_bound_beta$LG[bgc_bound_beta$chrom == "X" & bgc_bound_beta$mb < 100 ] <- "XPAR"
  
  beta_linkage <-   lm(data=bgc_bound_beta,formula = var_median ~ LG )
  summary(beta_linkage)
  

  ##find the beta outliers
  
  bgc_bound[!is.na(bgc_bound$outlier) &  bgc_bound$var =="beta" ,] # 2032  2034  5382 13296
  bgc_bound$SNP[!is.na(bgc_bound$outlier) &  bgc_bound$var =="beta" ] # 2032  2034  5382 13296
  
  bgc_bound_beta <- bgc_bound[bgc_bound$var == "beta" & bgc_bound$LG != "U",]
  title_b <- expression(paste(beta))
  
  s012_yin <- fread('GBS.mis60.maf01.bgc.012',na.strings = "-1")
  s012_yin <- s012_yin[,-1]
#9666 10468 13620 (+1) 9667 10469 13621
  
 # ggplot( cbind("V10469"= s012_yin$V13621,indsInfo),aes(x=Longitude,y=flower,color=as.factor(V13621))) + geom_jitter(size=2) + theme_bw() +
  #  scale_color_manual(values=c( '#ffd700', '#9d02d7'),na.value="lightgrey") + 
 #   labs(y="Sex",color="Genotype")

  
  V9667 <-  cbind("genotype"= s012_yin$V9667,indsInfo)
  V9667 <- V9667[V9667$flower != "Inc" & !is.na(V9667$genotype),]
  V10469 <-  cbind("genotype"= s012_yin$V10469,indsInfo)
  V10469 <- V10469[V10469$flower != "Inc" & !is.na(V10469$genotype),]
  V13621 <-  cbind("genotype"= s012_yin$V13621,indsInfo)
  V13621 <- V13621[V13621$flower != "Inc" & !is.na(V13621$genotype),]

  beta_outliers_genotyped <- 
  rbind.data.frame(
    cbind.data.frame(V9667,"SNP"=9666,"SNPpos"="A3 3mb"),
    cbind.data.frame(V13621,"SNP"=13620,"SNPpos"="A3 41mb"),
    cbind.data.frame(V10469,"SNP"=10469,"SNPpos"="X 134mb")
    
  )
  
  
pdf(paste("fig5_beta.",today,".pdf",sep=''),width=5,height=4)


plot_grid(  
  
  ggplot() + 
    geom_line(data=bgc_bound_beta,color="black",alpha=0.1,size=0.25,aes(x=mb,y=var_median)) +
    
    geom_point(data=bgc_bound_beta,alpha=0.1,color="white",aes(x=mb,y=var_median)) +
  #  geom_point(data=bgc_bound_beta[!is.na(bgc_bound_beta$excess) ,],shape=4,aes(x=mb,y=var_median),color="gray") +
    geom_point(data=bgc_bound_beta[!is.na(bgc_bound_beta$outlier) ,],size=0.5,color="black",shape=8,aes(x=mb,y=var_median)) +
    
    geom_hline(yintercept=0, color = "gray", size=1,alpha=0.25) +  
    facet_grid( . ~ chrom,scales = "free",space="free_x") + 
    theme_bw()  + theme_bw(base_size = 8) + 
    labs(x="Position (mb)",y=title_b) +
    theme(strip.background =element_rect(fill="white"))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_color_manual(values=c( 'black', 'gray'))
  
  ,
  ggplot( beta_outliers_genotyped,aes(x=Longitude,y=as.factor(genotype),shape=flower,color=as.factor(genotype))) + 
    geom_vline(xintercept = -88.55, linetype="dashed",alpha=0.5) +
    geom_jitter(alpha=0.7,) + theme_bw(base_size = 8) +
    facet_grid(SNPpos~.,space="free",scale="free")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(x="Longitude",y="Genotype",shape="Sex",color="Genotype") +  scale_shape_manual(values=c(16,1))
  ,  labels = c('A', 'B'), label_size = 12,ncol = 1, align = 'vh',axis='tbrl',rel_heights = c(1,2))

dev.off()



####Substructure in XX/XY####

XY_PCA <- PCA(s012[indsInfo$riverside == "W",], graph = FALSE) #missMDA assumes continuous variables

XY_PCAcoord <- cbind.data.frame(XY_PCA$ind$coord,indsInfo[indsInfo$riverside == "W",-1])

##longitude##
ggplot(XY_PCAcoord,aes(y=Dim.1, x=Dim.2)) + 
  geom_point(size=3,alpha=1) + 
  theme_bw(base_size = 18) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  scale_shape_manual(values=c(3,5,4))

cor.test(XY_PCAcoord$Dim.1, XY_PCAcoord$Longitude, method=c("pearson"))


hybPCAcoord_2 <- hybPCAcoord
names(hybPCAcoord_2)[1:5] <- c("dim1_all","dim2_all","dim3_all","dim4_all","dim5_all")

XY_all_PC <- left_join(XY_PCAcoord,hybPCAcoord_2)

ggplot(XY_all_PC,aes(x=Dim.1,y=dim1_all)) + geom_point()
#ggplot(XY_all_PC,aes(x=Dim.2,y=dim1_all)) + geom_point()

cor.test(XY_all_PC$Dim.1, -XY_all_PC$dim1_all, method=c("pearson"))


##landscape prediction##
worldClim <- fread('~/Google Drive/Research/Data2/spider/worldClimBio.txt',header=FALSE)

ind_latlon <- indsInfo[indsInfo$riverside == "W",c(9:8)]

coordinates(ind_latlon) <- ~Longitude+Latitude
climate <- getData('worldclim', var='bio', res=2.5)

predictors <- stack(climate)
presvals <- raster::extract(predictors, ind_latlon)

bioClimbyPop <- cbind.data.frame(indsInfo[indsInfo$riverside == "W",],presvals,XY_PCAcoord)

M<-cor(bioClimbyPop[,c(8:9,14,18:41)])
corrplot(M, type="upper", order="hclust")


cor.test(-bioClimbyPop$Dim.1, bioClimbyPop$elevation, method=c("pearson"))
cor.test(bioClimbyPop$Dim.1, bioClimbyPop$Longitude, method=c("pearson"))


XYPCA_elev_lm_collec <-  lm(formula = Dim.1 ~ Longitude + Collec , data = bioClimbyPop) 
summary(XYPCA_elev_lm_collec)



####SNP associations####
XYPCAloc <- cbind.data.frame(loc,XY_PCA$var$contrib)

XYPCAloc$chrom = factor(XYPCAloc$LG, levels=c('A1','A2','A4','X','A3','U'))

  ggplot(XYPCAloc[XYPCAloc$LG != "U" ,]) + 
  geom_point(aes(y=Dim.4,x=mb,alpha=Dim.4),color="#ffd700") + 
  geom_point(aes(y=Dim.3,x=mb,alpha=Dim.3),color="#fa8775") + 
  geom_point(aes(y=Dim.2,x=mb,alpha=Dim.2),color="#cd34b5") + 
  geom_point(aes(y=Dim.1,x=mb,alpha=Dim.1),color="#0000ff") + 
    
  #  stat_smooth(span = 0.25) + 
  facet_grid( . ~ chrom,scales = "free",space="free_x") + theme_bw(base_size = 18) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Position (mb)",y="PC1 eff.") +
  theme(strip.background =element_rect(fill="white"))  +
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank())+ theme(axis.title.x = element_blank())





####LD A2####
#pos_conversion <- fread('GBS.mis60.maf01.012.chrom.pos')

#ld <- fread('GBS.A2.XY.plink.ld')
ld <- fread('GBS.XY.v2.ld.within')
names(ld) <- c("CHR_A","BP_A","SNP_A","CHR_B","BP_B","SNP_B","R2")
#ld <- separate(ld, SNP_A , c("chrom_1","pos_1"), sep = ":", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")
#ld <- separate(ld, SNP_B , c("chrom_2","pos_2"), sep = ":", remove = TRUE, convert = FALSE, extra = "merge", fill = "left")

#ld <- ld[ld$chrom_2  ==  ld$chrom_1  ,]
ld_pos_chrom <- ld

ld_pos_chrom$CHR_A[ld_pos_chrom$CHR_A == "L.10"] <- "X"
ld_pos_chrom$CHR_A[ld_pos_chrom$CHR_A == "L.7"] <- "A1"
ld_pos_chrom$CHR_A[ld_pos_chrom$CHR_A == "L.8"] <- "A2"
ld_pos_chrom$CHR_A[ld_pos_chrom$CHR_A == "L.5"] <- "A3"
ld_pos_chrom$CHR_A[ld_pos_chrom$CHR_A == "L.3"] <- "A4"

#ld_pos_chrom$pos_1 <- as.numeric(ld_pos_chrom$pos_1)
#ld_pos_chrom$pos_2 <- as.numeric(ld_pos_chrom$pos_2)

maxA3pos <- 175014771
ld_pos_chrom$BP_A[ld_pos_chrom$CHR_A == "A3"] <- maxA3pos - ld_pos_chrom$BP_A[ld_pos_chrom$CHR_A == "A3"] 
ld_pos_chrom$BP_B[ld_pos_chrom$CHR_A == "A3"] <- maxA3pos - ld_pos_chrom$BP_B[ld_pos_chrom$CHR_A == "A3"] 

ld_pos_chrom$logR2 <- log(ld_pos_chrom$R2)

ggplot(ld_pos_chrom,aes(x=R2,fill=CHR_A))+ geom_histogram(position="stack")

ggplot(ld_pos_chrom,aes(x=logR2,fill=CHR_A))+ geom_histogram(position="stack") + geom_vline(xintercept = -7)


ld_pos_chrom_cut <- ld_pos_chrom[ld_pos_chrom$R2 > 0.5]

ggplot(ld_pos_chrom_cut,aes(x=BP_A/1000000,y=BP_B/1000000,color=R2)) + geom_point(size=0.01,shape=15) + 
  facet_grid(. ~ CHR_A) +
  theme_bw() +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_gradient2(low = "white",mid="yellow",high="blue")

####Fig S4: EEMS####
##locations export##
HZ <- indsInfo[indsInfo$state %in% c("MS","LA","AR","AL"),]
write.table(HZ$name, file = "HZEEMS.indv", append = FALSE, quote = FALSE, sep = "\t", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = FALSE, qmethod = c("escape", "double"))


HZ.coord <- HZ[,c(8:9)]
write.table(format(HZ.coord, digits=4), file = "HZEEMS.coord", append = FALSE, quote = FALSE, sep = "\t", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = FALSE, qmethod = c("escape", "double"))

ggplot() + 
  geom_path(data = states, aes(x = long, y = lat, group = group),color="lightgray") +
  geom_path(data = rivers, aes(x = long, y = lat, group = group),color="lightblue",size=2,alpha=0.75) +
  
  geom_point(data= HZ.coord, aes(x=Longitude,y=Latitude),alpha=0.75) +
  
  xlim(-95,-85)+ ylim(30,34)

#start EEMS import##

eems_results <- file.path(c(
  "~/Google Drive/Research/Data2/hybridGenome/GBShyb-EEMS-nDemes800-chain1"
  ,"~/Google Drive/Research/Data2/hybridGenome/GBShyb-EEMS-nDemes800-chain2"
  ,"~/Google Drive/Research/Data2/hybridGenome/GBShyb-EEMS-nDemes800-chain3"
  ,"~/Google Drive/Research/Data2/hybridGenome/GBShyb-EEMS-nDemes800-chain4"
  ,"~/Google Drive/Research/Data2/hybridGenome/GBShyb-EEMS-nDemes800-chain5"
  ,"~/Google Drive/Research/Data2/hybridGenome/GBShyb-EEMS-nDemes800-chain6"
))

projection_none <- "+proj=longlat +datum=WGS84"
projection_mercator <- "+proj=merc +datum=WGS84"

map_world <- getMap()
map_NA <- map_world[which(map_world@data$country == "United States"), ]
name_figures <- file.path(path.expand("~/Google Drive/Research/Data2/hybridGenome/eems/"), "eems_800d")

eems.plots(mcmcpath = eems_results,
           plotpath = paste0(name_figures, "-demes-and-edges"),
           longlat = F,
           add.grid=TRUE,
           col.grid="gray90",
           lwd.grid = 2,
           plot.width=10,
           add.outline = FALSE,
           #col.outline = "blue",
           #lwd.outline = 5,
           add.demes = TRUE,
           col.demes = "red",
           pch.demes = 5,
           min.cex.demes = 0.5,
           max.cex.demes = 1.5,
           out.png=F,
           projection.in = projection_none,
           projection.out = projection_mercator,
           add.map = TRUE,
           col.map = "black",
           lwd.map = 5,
           m.plot.xy = { plot(map_NA, col = NA, add = TRUE) },
           q.plot.xy = { plot(map_NA, col = NA, add = TRUE) }
)


####fertility analysis####
babies <- fread('hyb_offspringcounts.txt')

baby_var <- cbind.data.frame(
  aggregate(babies[, c(1,3)], list(babies$pop), mean),
  aggregate(babies[, c(1,3)], list(babies$pop), FUN=std.error)
)
baby_var <- baby_var[,c(1,3,6)]
names(baby_var) <- c("pop","baby_mean","baby_se")

baby_var_pop <- data.frame(sqldf('select pop_pos.*, baby_var.* from pop_pos
                                 left join baby_var on baby_var.pop = pop_pos.Pop  '))

babies$coast <- factor(babies$coast, levels=c('W','E'))


baby_aov <- aov(offspring_count ~ coast * pop,data=babies)
summary(baby_aov)


babies_strc_pop <- left_join(strc_Pop,babies,by=c("momID"="motherID"))

babies_strc_pop <- babies_strc_pop[babies_strc_pop$variable == "k3" & babies_strc_pop$k == "3",]

baby_HI_mod <-  lm(formula = offspring_count ~ value , data = babies_strc_pop) 
summary(baby_HI_mod)

ggplot(babies_strc_pop,aes(x=value)) + geom_histogram()

babies_strc_pop$hyb[babies_strc_pop$value < 0.05] <- "XY"
babies_strc_pop$hyb[babies_strc_pop$value > 0.05 & babies_strc_pop$value < 0.95] <- "Hybrid"
babies_strc_pop$hyb[babies_strc_pop$value > 0.95] <- "XYY"

babies_strc_pop$hyb <- factor(babies_strc_pop$hyb, levels=c('XY','Hybrid',"XYY"))

babies_strc_pop <- babies_strc_pop[complete.cases(babies_strc_pop), ]

baby_str_aov <- aov(offspring_count ~ hyb * pop,data=babies_strc_pop)
summary(baby_str_aov)


#baby FIT
babyFIT <-  fread('plink.het')
names(babyFIT)[3] <- "OHOM"
babies_strc_pop_FIT <- left_join(babies_strc_pop,babyFIT,by=c("inds"="FID"))


baby_str_fit_aov <- aov(FIT ~ hyb * pop,data=babies_strc_pop_FIT)
summary(baby_str_fit_aov)

baby_str_fit_aov <- aov(OHOM ~ hyb ,data=babies_strc_pop_FIT)
summary(baby_str_fit_aov)






