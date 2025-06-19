source('code/functions.R')

l2.orig <- read_sf('data/L2.dbf')
obs <- raster::stack("data/sbw_defol_stack.grd")

## Sample parameters
set.seed(100)
theta.lst <- post_parsamp(1)
theta.lst$est.prob <- 1
# $altitude
# [1] 0.6
# 
# $temp.min.to
# [1] 14.34946
# 
# $temp.max.to
# [1] 27.41412
# 
# $altitude.disp
# [1] 44.77669
# 
# $temp.min.disp
# [1] 15.11697

## Assume that for each year, there's one date on which all flights occur
## Let's assume it's July 15th every year

rs.all <- read.csv('data/res_simul_all.csv', header = T)
rs.all.cut1 <- subset(rs.all, Year > 12)
rs.all.cut2 <- subset(rs.all.cut1, Month < 8 | Day <= 6)
rm(rs.all, rs.all.cut1)
rs.all.cut2 <- data.table(rs.all.cut2)

all.dates <- unique(rs.all.cut2[,YMD])

rs.all.0715 <- rs.all %>%
  subset(Month == 7 & Day == 15 &
           Year > 12 &
         PBL == theta.lst$altitude) %>%
  data.table()

## Calculate endpoints based on selected date and parameters
endpoints.all <- post_eps(rs.all.0715, theta.lst$temp.min.to, theta.lst$temp.max.to, 
                          theta.lst$altitude.disp, theta.lst$temp.min.disp) #%>%
#   st_as_sf(coords = c('Lon', 'Lat'), crs = crs(l2.orig)) %>%
#   st_transform(crs = CRS('+init=epsg:6623'))

years <- unique(endpoints.all$Year)
est.prob <- 1
year.lst <- lapply(years, function(x) {
  endpoints.y <- subset(endpoints.all, Year == x)
  n.success <- round(nrow(endpoints.y)*est.prob)
  endpoints.samp <- endpoints.y[sample(1:nrow(endpoints.y), n.success),]
  ends.all <- as.data.frame(endpoints.samp[,c('Lon', 'Lat')])
  cord.dec <- SpatialPoints(ends.all, proj4string = CRS("+proj=longlat"))
  cord.utm <- as.data.frame(spTransform(cord.dec, CRS('+init=epsg:6623')))
  
  endpoints.samp[,c('x.coord', 'y.coord')] <- cord.utm
  
  ends.df <- endpoints.samp %>% st_as_sf(coords = c('x.coord', 'y.coord'))
  
  end.rast <- rasterize(ends.df, rast, rep(1, nrow(ends.df)), max, na.rm = TRUE)
  ## 35 percent of cells in original L2 data are unobserved so we'll replicate this
  values(end.rast) <- sapply(values(end.rast), function(x) {ifelse(runif(1) <= 0.35, NA,
                                                                   ifelse(is.na(x), 0, x))})
  
  return(end.rast)
})
names(year.lst) <- years
end.stack <- stack(year.lst)


##### Run ABC SMC #####

ptm.full.start <- proc.time()
ptm <- proc.time()
set.seed(101)
results.orig <- sample.n(2500, quant1 = 0.25, quant2 = 0.25, rast = year.lst) 
results <- results.orig
results$wvec <- 1/nrow(results)
results$index <- 1
print(proc.time() - ptm)
ptm <- proc.time()
results <- subset(results, l2hit > l2.ep & accuracy > acc.ep)
res.lst <- list(results)
for (i in 1:5) {
  print(proc.time() - ptm)
  ptm <- proc.time()
  res <- res.lst[[i]]
  set.seed(101 + i)
  results <- post_next.sample(res, index = i, q1 = 0.25, q2 = 0.25, rst = year.lst)
  results$index <- i + 1
  res.lst <- append(res.lst, list(results))
}

ptm.full.end <- proc.time()
ptm.full <- ptm.full.end - ptm.full.start

comp.df <- bind_rows(res.lst)
write.csv(comp.df, 'code/output/simulated_data_results.csv', row.names = FALSE)



