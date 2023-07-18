library(sf)
library(units)
library(vegan)

# Prep subsample data -----------------------------------------------------

# subsampled dataset for results reported in main text
subsmplFl <- 'Data/DS5_assemblage_data_subsampled.csv'
subsmpl <- read.csv(subsmplFl)
bins <- unique(subsmpl$bin)

# remove the 2 subsample iterations where all species were singletons
subsmpl <- subsmpl[!is.na(subsmpl$mst_med),]

# Prep global data --------------------------------------------------------

globlOccFl <- 'Data/DS1_unique_taxonLocationTime_occurrences.csv'
globlOcc <- read.csv(globlOccFl)

# calculate study variables over all occurrences in each time bin

# square-root minimum spanning tree length for each species
calcMst <- function(x, dat, nameCol, coordCols, prj){
  rows <- which(dat[, nameCol] == x)
  if (length(rows) < 2) { 
    return(NA) 
  } else {
    coordSf <- st_as_sf(dat[rows,], coords = coordCols, crs = prj)
    gcdists <- st_distance(coordSf) |> set_units(value = 'km')
    tr <- spantree(drop_units(gcdists))
    sum(tr$dist) |> sqrt()
  }
} 

# return study variables in one time bin
calcVars <- function(bin, nameCol, coordCols, prj){
  binBool <- globlOcc$stage == bin
  binDat <- globlOcc[binBool, ]
  # total species count
  binSpp <- unique(binDat[, nameCol])
  rich <- length(binSpp)
  
  # mean occupancy without singletons
  occSpp <- table(binDat[, nameCol]) |> as.numeric()
  occMean <- occSpp[occSpp > 1] |> mean()
  
  # median square-root MST
  mstSpp <- sapply(binSpp, calcMst, dat = binDat, nameCol = nameCol, 
                   coordCols = coordCols, prj = prj)
  mstMed <- median(mstSpp, na.rm = TRUE)
  
  cbind('species_count' = rich,
             'oc_wo.singl_mean' = occMean, 
             'mst_med' = mstMed)
}

# iterate over bins
globlList <- sapply(bins, calcVars, nameCol = 'unique_name', 
                coordCols = c('centroid_long', 'centroid_lat'),
                prj = 'ESRI:54010', # Eckhert VI
                simplify = FALSE # retain column names in output
                ) 
globlMeta <- do.call(rbind, globlList) |> data.frame()
globlMeta$bin <- names(globlList)
