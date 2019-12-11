#' Convex Hull Ensembles
#'
#' Core function for 
#' @param j Counter for bioclim raster position
#' @return A \code{list} comprising two stacks of rasters
#' @export


# modified function from original code
runCHE <- function(j){

  for (i in 1:11){ cat(i,'\n')

    VarX <- paste0("bio", i)
    VarY <- paste0("bio", j)

    # if (i < 10) {VarX <- paste0("bio0", i)}
    # if (j < 10) {VarY <- paste0("bio0", j)}

    predict.to.spp <- stack(predict.to[[VarX]], predict.to[[VarY]])

    natptsR <- cbind(natptsall[,c("x", "y")], natptsall[VarX], natptsall[VarY])
    backR <- cbind (backall[,c("x", "y")], backall[VarX], backall[VarY])

    pbg.envR <- rbind(natptsR, backR)
    pbg.envR <- pbg.envR[,3:ncol(pbg.envR)]
    spp.dataR <- cbind (pbg.which, pbg.envR)

    #### model
    # mgcv package
    mod.formula <- as.formula(paste('pbg.which/w ~ s(', VarX, ", bs='ts') + s(", VarY, ", bs='ts')"))
    temp <- mgcv::gam(mod.formula, data=spp.dataR, family='poisson', method='REML', weights=w)
    pred <- predict(temp, type='response')

    ## model evalution section
    roc.test <- roc(pred, PA2)
    auc.test <- auc(roc.test)

    ## predict to invasive range
    pred2 <- predict(predict.to, temp,type="response")

    idx <- paste(spp, i, j, sep="_")
    AUC.result[[idx]] <- auc.test
    names (pred2) <- idx
    varsStack <- stack(varsStack, pred2)
    train <- natptsR[,c(3,4)]

    ## create convex hull in environmental space
    ch <- convHull(train)
    predict.to.spp.pts <- rasterToPoints(predict.to.spp)
    ## this needs to be here for some predictor / data set combinations (e.g. l_humile and bio03)
    predict.to.spp.pts <- predict.to.spp.pts[complete.cases(predict.to.spp.pts),]
    # create 1s for inside hull, 0s for outside
    inside.hull <- predict(ch, predict.to.spp.pts[,3:4])
    # add the environmental data (from geographic space) to the hull
    occur.pred <- cbind(predict.to.spp.pts[,1:2], inside.hull)

    temp <- rasterize(x = occur.pred[,1:2], y = predict.to.spp[[1]], field = occur.pred[,3])
    names (temp) <- paste(spp, i, j, sep="_")
    convStack <- stack(convStack, temp)
  }


  test <- list (convStack, varsStack, AUC.result)
}

#' Stacking function for che outputs
#'
#' @param stackoflayers Raster layers of che results
#' @param keep_layers Vector of layer names to keep
#' @return A single layer
#' @export

sum_stack <- function(stackoflayers, keep_layers=keep_layers){
  x <- stack()
  for (i in keep_layers){
    x <- stack(x, stackoflayers[[i]])
  }
  x <- sum(x)/nlayers(x)
  return (x)
}

#' Prepare data for use in the CHE model
#'
#' @param rasters Stack of rasters that are bioclim variables and cover extent of study (global fine)
#' @param mask_raster Background mask, usualyl created using \code{background_builder}
#' @param predict.to Raster stack where the model will be predicted/ projected to
#' @param spp.dist X & Y species points
#' @param w.val Weighting value for background
#' @param bg.sample Size of background samples (default 50000)
#' @return A list of data objects to be used with che model.
#' @export

prepare_data <- function(rasters = rasters, 
                         mask_raster = mask, 
                         predict.to = predict.to, 
                         spp.dist=spp.dist, 
                         w.val=1.e-6,
                         bg.sample = 50000){
  
  ## create a stack of predictors masked to training region
  train.stack <- crop (rasters, mask_raster)
  train.stack <- mask (train.stack, mask_raster)
  
  ## Rasterize the points to the same
  train.dist <- rasterize (spp.dist, train.stack)
  train.dist <- reclassify (train.dist, c(0, Inf, 1))
  
  #make a stack of layers masked to the rasterized distribution
  bioclimPres <- mask (train.stack, train.dist)
  
  #turn this into a data.frame
  natptsall <- as.data.frame(rasterToPoints(bioclimPres, xy=TRUE))
  natptsall <- natptsall[complete.cases(natptsall),]
  
  ## keep all the presence points in the background (ignore P/A at this point)
  backall <- as.data.frame(rasterToPoints(train.stack, xy=TRUE))
  backall <- backall[complete.cases(backall),]
  
  if (length(backall[,1]) > bg.sample){
    backall <- dplyr::sample_n(backall, bg.sample, replace=FALSE)
  }
  
  # Create data.frames used for the modelling
  pbg.env <- rbind(natptsall, backall)
  pbg.env <- pbg.env[,3:ncol(pbg.env)]
  pbg.which <- c(rep(1, nrow(natptsall)), rep(0, nrow(backall)))
  spp.data <- data.frame(cbind(pbg.which, pbg.env))
  PA <- pbg.which
  PA2 <- as.factor (PA)
  
  ## background weighting
  cell_size<-area(mask_raster, na.rm=TRUE, weights=FALSE)
  cell_size<-cell_size[!is.na(cell_size)]
  raster_area<-length(cell_size)*median(cell_size)
  back.area <- raster_area
  w <- rep(1.e-6, nrow(spp.data))
  w[spp.data$pbg.which ==0] <- back.area / sum(spp.data$pbg.which == 0)
  
  out_list <- list(predict.to, PA, PA2, w, natptsall, pbg.which, backall)
  names(out_list) <- c("predict.to", "PA", "PA2", "w", "natptsall", "pbg.which","backall")
  
  return(out_list)
}

#' Perform convex hulls models
#'
#' Wrapper function for the CHE modelling method
#' @param spp Name of species
#' @param parallel T/F flag for using parallel - must be TRUE for now
#' @param model.data A model.data list of formatted data from \code{prepare_data} function.
#' @return Modelled outputs
#' @export
#' 
che_model <- function(spp = "spp", 
                      parallel = TRUE,
                      model.data = model.data){
  
  #model.data <- lapply (model.data, list)
  #create lists, unpack model data
  varsStack <- stack()
  convStack <- stack()
  AUC.result <- list()
  
  

  predict.to <- model.data$predict.to
   PA <- model.data$PA
   PA2 <- model.data$PA2
   pbg.which <- model.data$pbg.which
  
   w <- model.data$w
   natptsall <- model.data$natptsall
   backall <- model.data$backall
  
  ### run models using runCHE function
  if (parallel==TRUE){
    
    require(parallel)
    
    n.cores <- detectCores() - 1
    cl <- makeCluster(n.cores, type="FORK", methods=F)
    clusterExport(cl=cl, varlist=c('predict.to', 'convStack', 'PA', 'PA2', 'w', 'natptsall', 'pbg.which', 'backall', 'AUC.result', 'spp', 'varsStack'),
                  envir=environment())
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
   return(test)
}







