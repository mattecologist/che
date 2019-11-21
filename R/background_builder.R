#' Background Builder for Species Distribution Modelling
#'
#' Before this function is used you'll need to load wrld_simply
#' @param distrib Distribution data with "Longitude" and "Latitude" colnames
#' @param terres The terrestrial biomes layer from Olson et al.
#' @param wrld_simpl The simple country/ continent boundaries from maptools
#' @param ref_rast Reference raster for grid cell processing
#' @param no_biome Flag to turn off the biome part and just return continent level backgrounds
#' @return A \code{raster} file that can be used for sampling background
#' @export

background_builder <- function(distrib=distrib,
                               terres=terres,
                               wrld_simpl = wrld_simpl,
                               ref_rast = ref_rast,
                               no_biome = FALSE){

  ## First part, make the XY data into coordinates (spatial points)
  coordinates(distrib) <- ~Longitude+Latitude
  distrib@proj4string <- CRS("+proj=longlat +datum=WGS84")

  # # #determine unique CONTINENTS
  wrld_simpl@proj4string <- CRS("+proj=longlat +datum=WGS84")

  continents <- na.exclude (over (distrib, wrld_simpl))

  continents <- unique(continents$SUBREGION)

  continentspoly <- wrld_simpl[wrld_simpl$SUBREGION %in% continents, ]
  continentspoly <- gUnaryUnion(continentspoly)
  projection(continentspoly) <- CRS('+proj=longlat')


  if (no_biome==FALSE){

    ## Extract which biomes the distribution points fall on, create a polygon with this.
    terres@proj4string <- CRS("+proj=longlat +datum=WGS84")
    biomes <- over (distrib, terres)
    biomes <- unique(biomes$BIOME); biomes <- biomes[which(!is.na(biomes))]
    biomespoly <- terres[terres$BIOME %in% biomes, ]

    #break down polygons
    biomespoly <- gUnaryUnion(biomespoly)
    projection(biomespoly) <- CRS('+proj=longlat')

    distpoly <- gIntersection (biomespoly, continentspoly, byid=TRUE)
  } else{
    distpoly <- continentspoly}


  dist_ras <- rasterize (distpoly, ref_rast)
  dist_ras <- trim(dist_ras)

  return (dist_ras)
}

#' Background weighting function
#'
#' @param bg_ras The background raster to sample across
#' @param weight The value for background cells
#' @param spp.data The presence data - must include a column of 0/1 for presence background samples
#' @return who knows
#' @export

bg_weighting <- function(bg_ras=bg_ras,
                         weight=1.e-6,
                         spp.data=spp.data){
  cell_size<-area(bg_ras, na.rm=TRUE, weights=FALSE)
  cell_size<-cell_size[!is.na(cell_size)]
  raster_area<-length(cell_size)*median(cell_size)
  back.area <- raster_area
  w <- rep(weight, nrow(spp.data))
  w[spp.data$pbg.which ==0] <- back.area / sum(spp.data$pbg.which == 0)
  return (w)
}


