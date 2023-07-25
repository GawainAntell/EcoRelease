library(divvy) # uniqify and sdsumry functions
library(fUnitRoots) # for timeseries stationarity tests
library(NSM3) # for bootstrapped confidence intervals
# table and plot outputs
library(stringr)
library(xtable)
options(xtable.timestamp = "")
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(velociraptr) # geologic timescale for time series plot

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

# adapted from R/correlate_Tseries_GSA.R

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

# adapted from R/correlate_Tseries_GSA.R

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

# (2) site number

globlArea <- globlMeta$nLoc
adfTest(globlArea)
area_AR <- arima(globlArea, order = c(1, 0, 0))
area_e <- as.numeric(area_AR$residuals)
acf(area_e)
adfTest(area_e)
area_AR$coef
#       ar1   intercept 
# 0.4006415 102.0221630

# site number vs. richness
cor(area_e, globlRich, method = 'kend')
# [1] 0.4144468

kendall.ci(area_e, globlRich, bootstrap = TRUE, B = 10000)
# 1 - alpha = 0.95 two-sided CI for tau:
# 0.259, 0.554

# site number vs. proportional occupancy
cor(area_e, globlMeta$propOcc, method = 'kend')
# [1] -0.6682028

kendall.ci(area_e, globlMeta$propOcc, bootstrap = TRUE, B = 10000)
# 1 - alpha = 0.95 two-sided CI for tau:
# -0.769, -0.546 

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

# site aggregation vs. proportional occupancy

cor(agg_e, globlMeta$propOcc, method = 'kend')
# [1] -0.6507937

kendall.ci(agg_e, globlMeta$propOcc, bootstrap = TRUE, B = 10000)
# 1 - alpha = 0.95 two-sided CI for tau:
# -0.743, -0.544 

# * Sampling scatterplots -------------------------------------------------

artArea <- max(globlMeta$nLoc)
arty <- max(globlMeta$nTax)
p1 <- ggplot(globlMeta, aes(x = nLoc, y = nTax)) +
  theme_bw() +
  geom_point() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_x_continuous(limits = c(0, 385),  expand = c(0, 0),
                     labels = rep('', 5)) +
  scale_y_continuous(limits = c(0, 1500), expand = c(0, 0), 
                     labels = c(0, 500, 1000, '')) +
  annotate('text', x = artArea, y = arty, label = 'Artinskian', hjust = 1.2)

artMst <- max(globlMeta$minSpanTree)
p2 <- ggplot(globlMeta, aes(x = minSpanTree, y = nTax)) +
  theme_bw() +
  geom_point() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_x_continuous(labels = rep('', 5)) +
  scale_y_continuous(limits = c(0, 1500), expand = c(0, 0),
                     labels = rep('   ', 4)) +
  annotate('text', x = artMst, y = arty, label = 'Artinskian', hjust = 1.2)

gzhelArea <- min(globlMeta$nLoc)
gzhely <- max(globlMeta$propOcc)
p3 <- ggplot(globlMeta, aes(x = nLoc, y = propOcc)) +
  theme_bw() +
  geom_point() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank() 
  ) +
  scale_x_continuous(limits = c(0, 385),  expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.15), expand = c(0, 0)) +
  annotate('text', x = gzhelArea, y = gzhely, label = 'Gzhelian', hjust = -0.2)

gzhelMst <- globlMeta$minSpanTree[globlMeta$bin == 'Gzhelian']
p4 <- ggplot(globlMeta, aes(x = minSpanTree, y = propOcc)) +
  theme_bw() +
  geom_point() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_x_continuous(labels = c(30, 50, 70, 90, '')) +
  scale_y_continuous(limits = c(0, 0.15), expand = c(0, 0),
                     labels = rep('', 4)) +
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
    x = 0.205, vjust = 0.2, hjust = 0.2
  ) + 
  draw_label('Site dispersion (10^3 km)', # fontface = 'bold',
             x = 0.59, vjust = 0.2, hjust = 0
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

# * Range-vs-richness plots -----------------------------------------------

# PANEL 1: global range size vs. species count scatterplot

gzhelPO <- max(globlMeta$propOcc)
gzhelRich <- globlMeta$nTax[globlMeta$bin == 'Gzhelian']

pGlobl <- ggplot(data = globlMeta, aes(x = nTax, y = propOcc)) +
  theme_bw() +
  geom_point() +
  scale_x_continuous(limits = c(0, 1500), expand = c(0, 0), position = 'top',
                     labels = c(0, 500, 1000, '')) +
  scale_y_continuous(name = 'Proportion of global cells',
                     limits = c(0, 0.15), expand = c(0, 0)) +
  theme(axis.title.x = element_blank()) +
  annotate('text', x = gzhelRich, y = gzhelPO, label = 'Gzhelian', hjust = -0.2)

# PANEL 2: subsampled range size vs. species count scatterplot

# error bar function
barPos <- function(t, v, dat){
  rows <- which(dat[,'bin'] == t)
  ci <- quantile(dat[rows, v], probs = c(0.25,0.75), na.rm = TRUE)
  avg <- mean(dat[rows, v], na.rm = TRUE)
  out <- c(ci, avg)
}

# summarize mean and interquartile range for each variable
vars <- c('species_count', 'oc_wo.singl_mean') # , 'mst_med'
plotDatL <- lapply(vars, function(v){
  wide <- sapply(bins, barPos, v = v, dat = subsmpl)
  t(wide)
})
plotDatMat <- do.call(cbind, plotDatL)
plotDat <- data.frame(plotDatMat)
colnames(plotDat) <- sapply(vars, paste, c('lwr','upr','avg'), sep = '_')

lw <- 0.5 # interquanrtile range line width
hirnOcc <- max(plotDat$oc_wo.singl_mean_avg)
hirnRich <- plotDat$species_count_avg[ rownames(plotDat) == 'Hirnantian' ]

pSubsmpl <- ggplot(data = plotDat, aes(x = species_count_avg, y = oc_wo.singl_mean_avg)) + 
  theme_bw() +
  geom_errorbarh(aes(xmin = species_count_lwr, xmax = species_count_upr),
                linewidth = lw, colour = 'grey') +
  geom_errorbar(aes(ymin = oc_wo.singl_mean_lwr, ymax = oc_wo.singl_mean_upr),
                 linewidth = lw, colour = 'grey') +
  geom_point() +
  scale_x_continuous(limits = c(0, 225), expand = c(0, 0), # don't truncate 
                     position = 'top') + 
  scale_y_continuous(name = 'Count in subsample cells',
                     labels = c(2, 3, 4, '    5')) +
  theme(axis.title.x = element_blank()) +
  annotate('text', x = hirnRich, y = hirnOcc, label = 'Hirnantian', hjust = -0.2)
# NB: axes are flipped here compared to Fig 1 of Antell et al. (2020)
# That's because 2020 figure plotted y-axis as shared between panels A and B.
# Here all axis scales differ (different B subplot) and need to match axes above

# PANEL 3: time series of species count

# get geologic period ages for timescale axis
periods <- velociraptr::downloadTime('international periods')
firstStep <- which(periods$name == 'Ordovician')
periods <- periods[ 1:firstStep, ]
periods$abbrev[1] <- '' # not evough room to plot 'Q'
periods$duration <- periods$b_age - periods$t_age
periods <- periods[order(periods$b_age, decreasing = TRUE),]
periods$b_age <- - periods$b_age
periods$t_age <- - periods$t_age

# attach midpoint of stage/bin age for subsample and global data

stages <- velociraptr::downloadTime('international ages')
stages <- stages[order(stages$b_age, decreasing=TRUE), ]

# lump certain stages together
stages2omit <- c('Stage 2','Stage 3','Stage 4','Stage 5',
                 'Drumian','Guzhangian','Paibian','Jiangshanian',
                 'Stage 10',
                 'Floian','Darriwilian',
                 'Katian', # otherwise no seed cells for Sandbian
                 'Aeronian', # otherwise no Rhuddanian or Aeronian seed cells
                 'Homerian','Ludfordian',
                 'Pragian', 
                 'Eifelian',
                 'Bashkirian',
                 'Kasimovian', 
                 'Sakmarian','Kungurian', # no Artinskian species records or Sakmarian/Artinskian grain size
                 'Olenekian', # otherwise only 1 abundance datum for Olenekian
                 'Sinemurian', # Hettangian is too poorly sampled
                 'Bajocian', # otherwise no Aalenian seed cells
                 'Hauterivian','Barremian', # no seed cells for Haut., Barremian or Valanginian alone
                 'Santonian', # otherwise nothing survives from Coniacian
                 'Thanetian',
                 'Bartonian', # otherwise no environmental data for Bartonian
                 'Aquitanian', # otherwise no seeds here or in Chattian
                 'Serravallian', # otherwise no seed cells for Langhian
                 'Messinian', # otherwise no seed cells for Messinian
                 'Calabrian','Middle Pleistocene','Late Pleistocene') # otherwise weird extinction rates
stages <- stages[!(stages$name %in% stages2omit),] 

# save median age for each record

# most recent time bin (Pleistocene) should count as running up to present
age2adjust <- which(stages$name == 'Gelasian') +1
stages$b_age[age2adjust] <- 0

stages$time_mid <- NA
getMid <- function(row){
  startAges <- stages$b_age[c(row, row+1)]
  mean(startAges)
}
stages$time_mid[1:(nrow(stages)-1)] <- sapply(1:(nrow(stages)-1), getMid)

# remove Cambrian and Holocene/latest Pleistocene stages outside study interval
# Note: updates to velociraptr's retrieved timescale since 2020 study mean
# there are additional stages to remove here compared to the
# 'read_data_for_bootstrapping' R script
firstRow <- which(stages$name == 'Tremadocian')
lastRow <-  which(stages$name == 'Gelasian')
stages <- stages[firstRow:lastRow, ]

plotDat$t <- -stages$time_mid[ stages$name %in% rownames(plotDat) ]
globlMeta$t <- -stages$time_mid[ stages$name %in% globlMeta$bin ]

max_x <- min(periods$b_age)
min_y <- 0 
max_y <- max(globlMeta$nTax) * 1.05
mainLw <- 0.7 # width for time series line
attach(periods)

empty <- ggplot() + 
  theme_bw() + 
  labs(title = '') + # create space for geologic scale to go
  scale_x_continuous(expand = c(0, 0), limits = c(max_x, 0),
                     breaks = seq(0, max_x, by = -100),
                     labels = paste( -seq(0, max_x, by = -100) )
                     ) +
  scale_y_continuous(name = 'Species count',
                     expand = c(0, 0), limits = c(min_y, max_y) 
                     ) +
  theme(axis.title.x = element_blank())

tseries <- 
  empty +
  geom_line(data = plotDat, aes(x = t, y = species_count_avg), 
            lwd = mainLw) +
  geom_linerange(data = plotDat, aes(x = t, ymin = species_count_lwr, ymax = species_count_upr), 
                 linewidth = lw) +
  geom_line(data = globlMeta, aes(x = t, y = nTax), 
            lwd = mainLw, linetype = 'twodash') +
  
  # add opaque rectangular bands for alternate geo periods
  # note that these lines must come after some other object (e.g. line) is plotted
  annotate('rect', xmin = t_age[1],  xmax = b_age[1],  ymin = -Inf, ymax = Inf, alpha = 0.2) +
  annotate('rect', xmin = t_age[3],  xmax = b_age[3],  ymin = -Inf, ymax = Inf, alpha = 0.2) +
  annotate('rect', xmin = t_age[5],  xmax = b_age[5],  ymin = -Inf, ymax = Inf, alpha = 0.2) +
  annotate('rect', xmin = t_age[7],  xmax = b_age[7],  ymin = -Inf, ymax = Inf, alpha = 0.2) +
  annotate('rect', xmin = t_age[9],  xmax = b_age[9],  ymin = -Inf, ymax = Inf, alpha = 0.2) +
  annotate('rect', xmin = t_age[11], xmax = b_age[11], ymin = -Inf, ymax = Inf, alpha = 0.2) 

# add geologic time series scale bar

pg <- ggplotGrob(tseries) #convert to grob object
j.plot <- unique(gtable::gtable_filter(pg, "panel", trim = FALSE)$layout$l) 

# create empty gtable
gt <- gtable::gtable(widths = unit(periods$duration,'null'), heights = unit(1, "null"))

# fill gtable with individual table grobs for each geologic time bin
for(i in 1:nrow(periods)){
  tt <- tableGrob(d = periods$abbrev[i],
                  cols = NULL, rows = NULL,
                  heights = unit(1, "null"),
                  widths = unit(periods$duration[i],'null')
                  # theme = ttheme_minimal(core = list(bg_params = list(fill=NULL)))
  )
  gt <- gtable::gtable_add_grob(x = gt, grobs = tt, t = 1, l = i)
}

pg_w_axis <- gtable::gtable_add_grob(x = pg, grobs = gt, t = 3, l = j.plot) 

# COMBINE PANELS

# align the bottom-row plot (time series) with the left-most plot of the top row 
# because y-axis labels have different number of digits
pAlignLft <- align_plots(pGlobl, pg_w_axis, align = 'v', axis = 'l')

plotTop <- plot_grid(
  pAlignLft[[1]], pSubsmpl,
  ncol = 2, rel_widths = c(1,1),
  labels = c('AUTO'),
  hjust = -1.3, vjust = 1.6
)

# add x-axis titles
ttlTop <- ggdraw() +
  draw_label('Species count (regional or global)'  
             ) 
ttlBase <- ggdraw() +
  draw_label('Time (millions of years ago)',
             vjust = 0 
  )

plotTriad <- plot_grid(
  ttlTop, plotTop, pAlignLft[[2]], ttlBase,
  labels = c('', '', 'C', ''),
  hjust = -1.3, vjust = 1.6,
  ncol = 1, 
  rel_heights = c(0.06, 0.7, 0.4, 0.05)
)

# add y-axis title for top row
ttlLft <- ggdraw() +
  draw_label('Mean species\' occupancy', # fontface = 'bold',
             y = 0.65, angle = 90, 
             vjust = 0.5 # more positive = shifts right
  ) 

plotFull <- plot_grid(
  ttlLft, plotTriad,
  ncol = 2, rel_widths = c(0.04, 1)
)

day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%Y-%m-%d')
panelsNm <- paste0('scatterplots-tseries-richness-rangesize_', day, '.pdf')
pdf(panelsNm, width = 5.9, height = 5) # Palaeobiology 2-col width = 15cm
plotFull
dev.off()

