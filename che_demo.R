## Working example of using Convex Hull Ensembles and range bagging with the {che} package
## Matt Hill 2019

# this library is on my github
library (che)
# this data is the in the {che} package - will implement neater soon....
load("data/insect_dist.Rdata")
load("data/terres.Rdata")


library (raster)
library (maptools)
library (rgeos)

bioclimall <- raster::getData('worldclim', var='bio', res=2.5)


e <- extent (10, 40, -35, -15)
SAfrica <- crop (bioclimall, e)
#convert to a stack
SAfrica <- stack(unstack(SAfrica))




dist <- insect_dist[insect_dist$Species =="h_destructor" & insect_dist$Range == "Native",]

sp_df <- cbind(as.data.frame(dist[,c("Longitude", "Latitude")]),raster::extract(SAfrica, as.data.frame(dist[,c("Longitude", "Latitude")])))
plot (SAfrica[[1]])
points(as.data.frame(dist[,c("Longitude", "Latitude")]), pch=20)

#load shapefiles of continent/country and some biogeographical zones (Olsen et al. biomes / Koppen-Geiger)

#from maptools package
data("wrld_simpl")


# biomes available here
temp <- tempfile()
download.file("https://c402277.ssl.cf1.rackcdn.com/publications/15/files/original/official_teow.zip",temp)
unzip(temp)
unlink(temp)
terres <- readShapePoly("official/wwf_terr_ecos.shp")


#setwd for project (output directory for raster dump etc)
dir.create("output")
wd = "./output/"; setwd(wd)

clim <- cbind(as.data.frame(dist[,c("Longitude", "Latitude")]),raster::extract(SAfrica, as.data.frame(dist[,c("Longitude", "Latitude")])))
#clim <- dplyr::rename(clim, "Longitude"="Long", "Latitude"="Lat")


# Uses the background builder script to set up sampling environment
nat_ras <- background_builder(clim, terres, wrld_simpl, ref_rast = bioclimall[[1]])


## create a stack of predictors masked to projection region
predict.nat <- crop (bioclimall, nat_ras)
predict.nat <- mask (predict.nat, nat_ras)

predict.to <- SAfrica

###########
natraster <- rasterize (clim[,1:2], predict.nat)
natraster <- reclassify (natraster, c(0, Inf, 1))

bioclimPres <- mask (predict.nat, natraster)
natptsall <- as.data.frame(rasterToPoints(bioclimPres, xy=TRUE))
natptsall <- natptsall[complete.cases(natptsall),]


#bioclimAbs <- mask (predict.nat, natrasterAbs)
## keep all the presence points in the background (ignore P/A at this point)
backall <- as.data.frame(rasterToPoints(predict.nat, xy=TRUE))
backall <- backall[complete.cases(backall),]

if (length(backall[,1]) > 50000){
  backall <- dplyr::sample_n(backall, 50000, replace=FALSE)
}


pbg.env <- rbind(natptsall, backall)
pbg.env <- pbg.env[,3:ncol(pbg.env)]
pbg.which <- c(rep(1, nrow(natptsall)), rep(0, nrow(backall)))
spp.data <- data.frame(cbind(pbg.which, pbg.env))
PA <- pbg.which
PA2 <- as.factor (PA)

## background weighting
cell_size<-area(nat_ras, na.rm=TRUE, weights=FALSE)
cell_size<-cell_size[!is.na(cell_size)]
raster_area<-length(cell_size)*median(cell_size)
back.area <- raster_area
w <- rep(1.e-6, nrow(spp.data))
w[spp.data$pbg.which ==0] <- back.area / sum(spp.data$pbg.which == 0)


## Parallel processing of GAMs

varsStack <- stack()
convStack <- stack()
AUC.result <- list()


spp <- "h_destuctor"

parallel = TRUE

if (parallel==TRUE){

require(parallel)

n.cores <- detectCores() - 1
cl <- makeCluster(n.cores, type="FORK", methods=F)
clusterExport(cl=cl, varlist=c('predict.to', 'convStack', 'PA', 'PA2', 'w', 'natptsall', 'pbg.which', 'backall', 'AUC.result'))
clusterEvalQ(cl, {
  require(methods)
  require(AUC)
  require(dismo)
  require(mgcv)})

# Run models on the cluster (cl), for variables (j) using the runGams function...
test <- parLapply(cl, 12:19, runCHE)
stopCluster(cl)
} else{

  for (k in 12:19){
    runCHE(k)
    ## add list
}

}

  ###########################################
## reconcile outputs

##test[[i]][[1]] = convexhulls
##test[[i]][[2]] = GAMs
##test[[i]][[3]] = AUC values

## just reset these all first
AUC.result <- list()
convStack <- stack()
varStack <- stack ()

for (i in 1:8){
  convStack <- stack(convStack, test[[i]][[1]])
  varsStack <- stack(varsStack, test[[i]][[2]])

  temp.AUC <- unlist(test[[i]][[3]])
  AUC.result[[i]] <- temp.AUC

}


AUC.result <- unlist(AUC.result)


quantile_thresh <- (quantile (AUC.result, 0.25))
keep_AUC <- AUC.result[AUC.result > quantile_thresh]
keep_AUC <- tibble::rownames_to_column(as.data.frame(keep_AUC))



q_CHE <- sum.stack(convStack,keep_AUC[,1])

plot (q_CHE)
points(dist[,c("Longitude", "Latitude")], pch=20)


## Range bagging

devtools::install_version("geometry", version = "0.3-6", repos = "http://cran.us.r-project.org")

source("./Rcode/rb.R")

temp<- rb(x = natptsall[,3:21], d=2)

test <- as.data.frame(getValues(SAfrica))
temp2 <- rb.test(temp, test)

pred <-SAfrica[[1]]
pred[] <- temp2
points(as.data.frame(dist[,4:3]), pch=20)

library (rasterVis)
library (viridis)

pred <- mask(pred, SAfrica[[1]])

rasterVis::gplot(pred) + geom_tile(aes(fill = value)) +
  scale_fill_viridis(option="viridis", begin=0.0, direction=1, na.value = NA,
                     breaks=seq(0, 1, 0.25), limits=c(0,1))+
  coord_equal()+
  scale_x_continuous(expand = c(0,0), limits= c(10, 40)) +
  scale_y_continuous(expand = c(0,0), limits= c(-35, -15)) +
  geom_path (data=worldshp, aes(x=long, y=lat, group=group), alpha=0.5, size=0.1)+
  ggtitle("Range Bagging Prediction")+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_minimal()+
  theme(axis.text.y=element_text(size =12),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.x=element_text(size =12),
        legend.title = element_text(size=14),
        legend.text=element_text(size=12))


rasterVis::gplot(q_CHE) + geom_tile(aes(fill = value)) +
  scale_fill_viridis(option="viridis", begin=0.0, direction=1, na.value = NA,
                     breaks=seq(0, 1, 0.25), limits=c(0,1))+
  coord_equal()+
  scale_x_continuous(expand = c(0,0), limits= c(10, 40)) +
  scale_y_continuous(expand = c(0,0), limits= c(-35, -15)) +
  geom_path (data=worldshp, aes(x=long, y=lat, group=group), alpha=0.5, size=0.1)+
  ggtitle("Convex Hull Ensemble Prediction")+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_minimal()+
  theme(axis.text.y=element_text(size =12),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.x=element_text(size =12),
        legend.title = element_text(size=14),
        legend.text=element_text(size=12))

## tmap

library(tmap)
library (classInt)


c_dist <- dist[,4:3]
coordinates(c_dist) <- ~Long+Lat

tm_shape(q_CHE)+
  tm_raster(palette=viridis(10), n=10,
            style="pretty", title="Suitability")+
  tm_shape(wrld_simpl)+
  tm_borders(col="grey60", lwd=2)+
  tm_shape(c_dist)+
  tm_dots(size=0.5, col="red")+
  tm_layout(title="Convex Hull Ensembles", title.position = c("left", "BOTTOM"),
            legend.position=c("right", "BOTTOM"))
save_tmap(filename="che_map.png")

tm_shape(pred)+
  tm_raster(palette=viridis(10), n=10,
            style="pretty", title="Suitability")+
  tm_shape(wrld_simpl)+
  tm_borders(col="grey60", lwd=2)+
  tm_shape(c_dist)+
  tm_dots(size=0.5, col="red")+
  tm_layout(title="Range Bagging", title.position = c("left", "BOTTOM"),
            legend.position=c("right", "BOTTOM"))
save_tmap(filename="rb_map.png")

