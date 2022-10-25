# You need the following libraries installed to run this:
# - raster (depends on sp, possibly on rgdal)
# - gbm

remotes::install_github("SEEG-Oxford/seegSDM")


library('raster')
library('rgdal')
library(devtools)
library(seegSDM)
library(dismo)
library(gbm)
#source('brt.functions.R')


make.raster.stack <- function(raster.fnames) {
    
    # This function takes:
    #   - A list of raster file names
    #   - The name of the globcover raster file
    #   - The globcover channels to use as predictors
    # It returns a list containing:
    #   - A raster stack with all the predictors in it
    #   - The names of the rasters.

    raster.names <- c()

    # Load the rasters from files into memory as a raster stack.
    rasters <- c()
    for (n in raster.fnames) {
        print(paste('Loading raster file', n))
        # Load the raster.
        r <- raster::raster(n)
        # Record the name of the raster.
        raster.names <- c(raster.names, sub('.asc','',n))
        # Append the raster to the list of rasters.
        rasters <- c(rasters, r)
    }
    # Make a raster stack from the list of rasters. A raster stack is a special object 
    # from the 'raster' library (see the documentation in the docs folder) that efficiently
    # supports operations like extraction over all of the rasters.
    sk <- raster::stack(rasters)
    

    
    list(stack=sk, names=raster.names)
}

extract.as.data.frame <- function(raster.stack, raster.names, data.locations, data.response) {

    # This function takes:
    #     - A raster stack
    #     - The names of the rasters
    #     - The locations of the datapoints
    # It returns the extraction of the raster stack at the data locations, as an R data frame.

    print('Extracting raster data to data locations')
    # Do the extractions. The result will be an R matrix.
    d <- raster::extract(raster.stack,data.locations)
    
    # Convert the matrix to an R data frame.
    extraction.list <- list(response=data.response)
    for (i in 1:ncol(d)) {
        extraction.list[[raster.names[i]]] = d[1:nrow(d),i]
    }
    
    as.data.frame(extraction.list)
}

vector.to.raster <- function(stack, vector){

    # This function takes a raster stack and a vector of new values,
    # and returns a new raster with the new values in it.
    # It's used by gbm.object.to.raster.

    unstack <- raster::unstack(stack)
    r <- raster::raster(unstack[[1]])
    r@data@inmemory <- TRUE
    where.na <- is.na(raster::getValues(unstack[[1]]))
    vector[where.na] = NA
    raster::setValues(r, vector)
}

gbm.object.to.raster <- function(gbm.object, raster.stack, raster.names) {
    
    # This function takes:
    #     - A GBM object (the output of the BRT code)
    #     - A raster stack
    #     - The names of the rasters
    # It returns:
    #     - A new raster of the prediction map.

    x <- raster::getValues(raster.stack)
    n.pixels <- dim(x)[1]
    data.list <- list()
    for (i in 1:length(raster.names)){data.list[[raster.names[i]]] = x[1:n.pixels,i]}
    df <- as.data.frame(data.list)

    # Make the prediction in chunks to keep memory usage down.
    prediction <- c()
    n.chunks <- 100
    chunk.size <- floor(n.pixels/n.chunks)
    for (i in 1:n.chunks){
        print(paste('Generating section ',i,' of ',n.chunks,' of prediction map.'))
        new.prediction <- predict.gbm(gbm.object, df[((i-1)*chunk.size+1):(i*chunk.size),], n.trees=Inf)
        prediction <- c(prediction, new.prediction)
    }
    if (chunk.size*n.chunks<n.pixels){
        new.prediction <- predict.gbm(gbm.object, df[(chunk.size*n.chunks+1):n.pixels,], n.trees=Inf)
        prediction <- c(prediction, new.prediction)
    }
    vector.to.raster(raster.stack, prediction)
}

#############################################################################################################
#############################################################################################################


# What predictors are we using?

raster.fnames <- c('rasters/prec57mx_10k_CLEAN_africa.asc','rasters/mod_dem_10k_CLEAN_africa.asc','rasters/wd0107a0_10k_CLEAN_africa.asc','rasters/wd0114a0_10k_CLEAN_africa.asc')

# Read in the data locations and the presence/absence values.

absence.locations <- as.matrix(read.csv('arab_abs.csv'))
presence.locations <- as.matrix(read.csv('arab_pres.csv'))
data.locations <- rbind(presence.locations, absence.locations)
presence <- c(rep(1,nrow(presence.locations)), rep(0,nrow(absence.locations)))

# Load the rasters and extract them to the data locations.
stack <- make.raster.stack(raster.fnames)
extraction <- extract.as.data.frame(stack$stack, stack$names, data.locations, presence)

# Fit the BRT model
gbm.obj <- gbm.step(extraction, 2:ncol(extraction), 1)

# Make a map.
result.raster <- gbm.object.to.raster(gbm.obj, stack$stack, stack$names)
raster::image(result.raster)


#############################################################################################################

#Statistics

# For explanations of everything below, see docs/JAE_1390_sm_Online_Tutorial.doc.

# visualise how the variables were used.
gbm.plot(gbm.obj)

# ROC (area under the curve)
roc(presence, gbm.obj$fitted>.5)

# Kappa and correlation (note the mean values of each)
gbm.obj$cv.statistics

# Deviance
deviance <- calc.deviance(presence, gbm.obj$fitted)
print(deviance)

# more detail on variable use in the model.
gbm.obj$contributions