#' Background Builder for Species Distribution Modelling
#'
#'
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

