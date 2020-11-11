## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(che)
library (raster)
library (maptools)
library (rgeos)

# simple world outline from maptools package
data("wrld_simpl")

## -----------------------------------------------------------------------------
load("../data/insect_dist.Rdata")

## -----------------------------------------------------------------------------
bioclimall <- raster::getData('worldclim', var='bio', res=2.5)

## -----------------------------------------------------------------------------
e <- extent (10, 40, -35, -15)
SAfrica <- crop (bioclimall, e)
#convert to a stack
SAfrica <- stack(unstack(SAfrica))

## ----fig.width=5, fig.height=5------------------------------------------------
dist <- insect_dist[insect_dist$Species =="h_destructor" & insect_dist$Range == "Native",]

sp_df <- cbind(as.data.frame(dist[,c("Longitude", "Latitude")]),raster::extract(SAfrica, as.data.frame(dist[,c("Longitude", "Latitude")])))
plot (SAfrica[[1]])
points(as.data.frame(dist[,c("Longitude", "Latitude")]), pch=20)

## -----------------------------------------------------------------------------

if (!file.exists("official/wwf_terr_ecos.shp")){
  temp <- tempfile()
  download.file("https://c402277.ssl.cf1.rackcdn.com/publications/15/files/original/official_teow.zip",temp)
  unzip(temp)
  unlink(temp)
}
terres <- readShapePoly("official/wwf_terr_ecos.shp")

## ----fig.width= 5, fig.height=5-----------------------------------------------
clim <- cbind(as.data.frame(dist[,c("Longitude", "Latitude")]),raster::extract(SAfrica, as.data.frame(dist[,c("Longitude", "Latitude")])))

nat_ras <- background_builder(distrib=clim, 
                              terres=terres,
                              wrld_simpl= wrld_simpl, 
                              ref_rast = bioclimall[[1]])

plot (nat_ras)

## -----------------------------------------------------------------------------

model.data <- prepare_data(rasters=bioclimall, mask_raster=nat_ras, predict.to=SAfrica, 
             spp.dist=clim[,1:2], w.val=1.e-6, bg.sample = 50000)

## -----------------------------------------------------------------------------

model.out <- che_model(spp = "h_destuctor", parallel = TRUE, model.data=model.data)



## ----fig.width= 5, fig.height=5-----------------------------------------------

q_CHE <- che_out(model.out, 0.25)


## -----------------------------------------------------------------------------
hist(q_CHE$AUC.scores)
abline(v=q_CHE$AUC.thresh, col="red")

## -----------------------------------------------------------------------------

library(tmap)
library (classInt)
library (viridis)


c_dist <- dist[,2:3]
coordinates(c_dist) <- ~Longitude+Latitude

tm_shape(q_CHE$predicted)+
  tm_raster(palette=viridis(10), n=10,
            style="pretty", title="Suitability")+
  tm_shape(wrld_simpl)+
  tm_borders(col="grey60", lwd=2)+
  tm_shape(c_dist)+
  tm_dots(size=0.3, col="red")+
  tm_layout(legend.position=c("right", "BOTTOM"))

