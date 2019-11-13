#' Convex Hull Ensembles
#'
#'
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
