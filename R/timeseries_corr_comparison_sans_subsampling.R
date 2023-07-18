library(divvy) # uniqify and sdsumry functions
# for stationarity and cross-correlation tests on time series
library(fUnitRoots) 
library(TSA)
# table and plot outputs
library(stringr)
library(xtable)
options(xtable.timestamp = "")

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

coordCols <- c('centroid_long', 'centroid_lat')
prj <- 'ESRI:54010' # Eckhert VI

# calculate study variables over all occurrences in each time bin
globlOcc <- uniqify(globlOcc, taxVar = 'unique_name', xy = coordCols)

globlList <- sapply(bins, function(b){
  binBool <- globlOcc$stage == b
  sdsumry(globlOcc[binBool,], taxVar = 'unique_name', 
          xy = coordCols, crs = prj)
  }, simplify = FALSE)
globlMeta <- do.call(rbind, globlList) |> data.frame()

# mean proportional occupancy without singletons, all taxa in one time bin
getPropOcc <- function(bin, nameCol, coordCols, prj){
  binBool <- globlOcc$stage == bin
  binDat <- globlOcc[binBool, ]
  nSite <- unique(binDat$cell_number) |> length()
  occSpp <- table(binDat[, nameCol]) |> as.numeric()
  occProp <- occSpp[occSpp > 1] / nSite
  mean(occProp)
}
# iterate over bins
globlMeta$propOcc <- sapply(bins, getPropOcc, nameCol = 'unique_name', 
                            coordCols = coordCols, prj = prj
                            ) 
globlMeta$bin <- names(globlList)

# Subsample correlations table --------------------------------------------

# determine the number of successful bootstrap replications
n <- sapply(bins, function(t) sum(subsmpl$bin==t))
nruns <- min(n)

resp_vars <- c('oc_wo.singl_mean', 'species_count') # 'mst_med'
pred_vars <- c('species_count','mst_sample') # 'pool_n'

cors <- matrix(nrow = length(resp_vars), ncol = length(pred_vars),
               dimnames = list(resp_vars, pred_vars))

# extract the i-th subsample from every time step
get_seq <- function(t, iter, cols, dat){
  series_rows <- which(dat$bin == t)[iter]
  dat[series_rows, cols]
} 

# Non-parametric correlation for a single time series at a time
# This function will run in parallel across bootstrap replicates,
# within a loop across time bins
get_cor <- function(pred, r_var, dat, iter){ 
  cols <- c(pred, r_var)
  series_list <- lapply(bins, get_seq, iter=iter, cols=cols, dat=dat)
  series <- matrix(unlist(series_list), ncol=length(cols), 
                   dimnames=list(bins, cols), byrow=TRUE)
  
  # fit AR1 model to the predictor, non-stationary series
  AR <- arima(series[,pred], order=c(1, 0, 0)) 
  series[,pred] <- as.numeric(AR$residuals) 
  
  # return only the point estimate of the coefficient
  # error bars come from bootstrapping, so no need to save significance values here
  cor(series[,c(r_var,pred)], method='kend')[1,2] 
} 

for (k in 1:length(pred_vars)){
  pred <- pred_vars[k]
  
  for (i in 1:length(resp_vars)){ 
    r_var <- resp_vars[i]
    
    # calculate across all iterations
    cors_i <- sapply(1:nruns, get_cor, pred = pred, r_var = r_var, dat = subsmpl)
    
    # 95% CI for 2-sided tests
    cors_avg <- mean(cors_i, na.rm = TRUE)
    cors_CI <- quantile(cors_i, probs=c(0.025, 0.975), na.rm=TRUE)
    
    # round (preserve trailing zeroes) and combine into 1 string
    roundSpcl <- function(x){
      sprintf("%.2f", round(x, digits=2))
    }
    rounded <- sapply(c(cors_avg, cors_CI), roundSpcl) 
    cors[i,k] <- paste(rounded[1], '_(', paste(rounded[2:3], collapse=', '), ')', sep='')
  } # end loop through range metrics
} # loop through each predictor variable

# format and save output as LaTeX table

sp_count_cors <- unlist(str_split(cors[,'species_count'], '_'))
sp_count_m <- matrix(sp_count_cors, ncol=2, byrow=TRUE)
agg_cors <- unlist(str_split(cors[,'mst_sample'], '_'))
agg_m <- matrix(agg_cors, ncol=2, byrow=TRUE)
fancy_cors <- cbind(resp_vars, sp_count_m, agg_m)

tabl <- xtable(fancy_cors, type = "latex", align=c('l','l','r','c','r','c') )  
caption(tabl) <- "Correlation between subsampled time series"
addtorow <- list()
addtorow$pos <- list(0, 0, 0)
addtorow$command <- c("& \\multicolumn{4}{c}{Correlation with predictor time series} \\\\ \\cline{2-5} \n",
                      "& \\multicolumn{2}{c}{Species count} & \\multicolumn{2}{c}{Sampling aggregation} \\\\\n",
                      "Response series & Mean & 95\\% & Mean & 95\\% \\\\\n")
day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%Y-%m-%d') 
name_rich <- paste0('tau_range_vs_richness_subsampled_', day, '.tex')
noDot <- function(str){ gsub("_", " ", str, fixed = TRUE) }
print(tabl, add.to.row = addtorow, file = name_rich, 
      caption.placement = 'top', booktabs = TRUE,
      include.colnames = FALSE, include.rownames = FALSE, 
      sanitize.text.function = noDot)
