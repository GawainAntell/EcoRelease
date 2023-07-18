library(divvy) # uniqify and sdsumry functions
library(fUnitRoots) # for timeseries stationarity tests
library(NSM3) # for bootstrapped confidence intervals
# table and plot outputs
library(stringr)
library(xtable)
options(xtable.timestamp = "")
library(ggplot2)
library(cowplot)

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

# Global correlations -----------------------------------------------------

# * Time series analysis --------------------------------------------------

# (1) richness vs. proportional occupancy

# check stationarity of predictor series (and fit AR1 model to achieve it)

globlRich <- globlMeta$nTax
adfTest(globlRich) 
acf(globlRich) 
rich_AR <- arima(globlRich, order = c(1, 0, 0)) # fit AR1 model
rich_e <- as.numeric(rich_AR$residuals)
acf(rich_e)
adfTest(rich_e) 
rich_AR$coef
#       ar1   intercept 
# 0.6368807 420.7497075 

cor(rich_e, globlMeta$propOcc, method = 'kend')
# [1] -0.4203789

kendall.ci(rich_e, globlMeta$propOcc, bootstrap = TRUE, B = 10000) # 95% CI
# 1 - alpha = 0.95 two-sided CI for tau:
# -0.553, -0.276 

# (2) site number vs. richness

globlArea <- globlMeta$nLoc
adfTest(globlArea)
area_AR <- arima(globlArea, order = c(1, 0, 0))
area_e <- as.numeric(area_AR$residuals)
acf(area_e)
adfTest(area_e)
area_AR$coef
#       ar1   intercept 
# 0.4006415 102.0221630

cor(area_e, globlRich, method = 'kend')
# [1] 0.4144468

kendall.ci(area_e, globlRich, bootstrap = TRUE, B = 10000)
# 1 - alpha = 0.95 two-sided CI for tau:
# 0.259, 0.554

# (3) site aggregation

globlAgg <- globlMeta$minSpanTree
adfTest(globlAgg)
agg_AR <- arima(globlAgg, order = c(1, 0, 0))
agg_e <- as.numeric(agg_AR$residuals)
acf(agg_e)
adfTest(agg_e)
agg_AR$coef
#          ar1    intercept 
# 3.736689e-01 6.346312e+04

# site aggregation vs. richness

cor(agg_e, globlRich, method = 'kend')
# [1] 0.4441599

kendall.ci(agg_e, globlRich, bootstrap = TRUE, B = 10000)
# 1 - alpha = 0.95 two-sided CI for tau:
# 0.313, 0.561 

# cross-correlation functions (use package TSA)
# pw1 <- prewhiten(rich, span, x.model=rich_AR) 
# pw2 <- prewhiten(rich, occ, x.model=rich_AR) 

# * Scatterplots ----------------------------------------------------------

artArea <- max(globlMeta$nLoc)
arty <- max(globlMeta$nTax)
p1 <- ggplot(globlMeta, aes(x = nLoc, y = nTax)) +
  theme_bw() +
  geom_point() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_x_continuous(labels = rep('', 5)) +
  annotate('text', x = artArea, y = arty, label = 'Artinskian', hjust = 1.2)

artMst <- max(globlMeta$minSpanTree)
p2 <- ggplot(globlMeta, aes(x = minSpanTree, y = nTax)) +
  theme_bw() +
  geom_point() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_x_continuous(labels = rep('', 5)) +
  scale_y_continuous(labels = rep('   ', 4)) +
  annotate('text', x = artMst, y = arty, label = 'Artinskian', hjust = 1.2)

gzhelArea <- min(globlMeta$nLoc)
gzhely <- max(globlMeta$propOcc)
p3 <- ggplot(globlMeta, aes(x = nLoc, y = propOcc)) +
  theme_bw() +
  geom_point() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank() 
  ) +
  annotate('text', x = gzhelArea, y = gzhely, label = 'Gzhelian', hjust = -0.2)

gzhelMst <- globlMeta$minSpanTree[globlMeta$bin == 'Gzhelian']
p4 <- ggplot(globlMeta, aes(x = minSpanTree, y = propOcc)) +
  theme_bw() +
  geom_point() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_x_continuous(labels = c(30, 50, 70, 90, '')) +
  scale_y_continuous(labels = rep('', 4)) +
  annotate('text', x = gzhelMst, y = gzhely, label = 'Gzhelian', hjust = -0.2)

quadPlot <- plot_grid(p1, p2, p3, p4,
                      labels = c('AUTO'),
                      hjust = -1.25, vjust = 1.6,
                      ncol = 2,
                      align = 'hv'
                      )

# add x-axis titles
ttlX <- ggdraw() + 
  draw_label('Site count', # fontface = 'bold',
    x = 0.205, hjust = 0
  ) + 
  draw_label('Site dispersion (10^3 km)', # fontface = 'bold',
             x = 0.59, hjust = 0
  )
plotaddx <- plot_grid(
  quadPlot, ttlX,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(1, 0.04)
)

# add y-axis titles
ttlY <- ggdraw() +
  draw_label('Proportional occupancy', # fontface = 'bold',
             y = 0.29, angle = 90
             ) +
  draw_label('Species count', # fontface = 'bold',
             y = 0.78, angle = 90)

plotaddxy <- plot_grid(
  ttlY, plotaddx,
  ncol = 2, rel_widths = c(0.04, 1)
)

day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%Y-%m-%d')
panelsNm <- paste0('scatterplots-global-bias_', day, '.pdf')
pdf(panelsNm, width = 5.9, height = 5.75) # Palaeobiology 2-col width = 15cm
plotaddxy
dev.off()
