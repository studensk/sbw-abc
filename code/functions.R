dfapply <- function(...) {
  l <- lapply(...)
  df <- as.data.frame(do.call('rbind', l))
  return(df)
}

##### MODEL SET-UP #####
library(lubridate)
library(parallel)
#library(tictoc)
library(raster)
#library(rgdal)
library(sf)
library(ClustGeo)
library(dplyr)
library(caret)
library(data.table)

## Read in and clean observation data; EPSG 6623 is the mapping for measurements in meters in Quebec
# Allows us to calculate distances in meters
l2 <- read_sf('data/L2.dbf') %>%
  st_transform(crs = CRS('+init=epsg:6623'))

#x.mat <- l2[,grep('X201', names(l2))]
x.mat <- l2[,grep('^201', names(l2))]
x.mat <- x.mat[,2:ncol(x.mat)] 
imm.mat <- l2[,grep('Imm_201', names(l2))]
all.mat <-  x.mat*imm.mat
rast.df2 <- all.mat %>%
  st_as_sf() %>%
  st_drop_geometry() %>%
  mutate(geometry = l2$geometry) %>%
  st_as_sf()


ext.ind <- st_bbox(rast.df2)
ext.ind.10k <- sapply(1:length(ext.ind), function(x) {
  if (substr(names(ext.ind)[x], 2, 4) == 'min') {
    floor(ext.ind[x]/10000)*10000
  }
  else {ceiling(ext.ind[x]/10000)*10000}
})
rast <- raster(xmn = ext.ind.10k['xmin'], xmx = ext.ind.10k['xmax'],
               ymn = ext.ind.10k['ymin'], ymx = ext.ind.10k['ymax'],
               resolution = c(10000, 10000)) 

rlst <- lapply(2013:2018, function(x) {
  colnm <- paste0('X', x)
  rasterize(rast.df2, rast, field = colnm, fun = 'sum')
})

names(rlst) <- 2013:2018
rstack <- stack(rlst)

## Read in start point rasters
obs <- raster::stack("data/sbw_defol_stack.grd")

## Define years, boundary layer levels and initial parameter boundaries
years <- 2007:2017
pbls <- c(0.4, 0.6, 0.8, 1, 1.2)
lower <- c(0.4, 13, 27, 0, 13)
upper <- c(1.2, 17, 31, 100, 17)
param.names <- c('altitude', 'temp.min.to', 'temp.max.to', 'altitude.disp', 'temp.min.disp')



params <- c('altitude', 'temp.min.to', 'temp.max.to',
            'altitude.disp', 'temp.min.disp', 'est.prob')

sample.discrete <- function(x, prob) {
  prob <- prob/sum(prob)
  cprob <- cumsum(prob)
  randos <- runif(x)
  values <- sapply(randos, function(r) {
    min(which(cprob > r))
  })
  return(values)
}

tkernel.sample <- function(theta, indx) {
  ind <- 1.5
  r.alt <- (round(rnorm(1, theta$altitude, 0.2/ind)/0.2 - 1) + 1)*0.2
  while(!(r.alt %in% seq(0.4, 1.2, by = 0.2))) {
    r.alt <- (round(rnorm(1, theta$altitude, 0.2/ind)/0.2 - 1) + 1)*0.2
  }
  theta$altitude <- r.alt
  theta$temp.min.to <- rnorm(1, theta$temp.min.to, 1/ind)
  theta$temp.max.to <- rnorm(1, theta$temp.max.to, 1/ind)
  r.ad <- rnorm(1, theta$altitude.disp, 16/ind)
  while(r.ad < 0) {
    r.ad <- rnorm(1, theta$altitude.disp, 16/ind)
  }
  theta$altitude.disp <- r.ad
  r.ep <- rnorm(1, theta$est.prob, 0.1/ind)
  while(r.ep < 0 | r.ep > 1) {
    r.ep <- rnorm(1, theta$est.prob, 0.1/ind)
  }
  theta$temp.min.disp <- rnorm(1, theta$temp.min.disp, 1/ind)
  theta$est.prob <- r.ep
  return(theta)
}

tkernel.density <- function(theta.obs, theta.samp, ind) {
  density <- theta.samp[params]
  density$altitude <- dnorm(theta.samp$altitude,
                            theta.obs$altitude, 0.2/ind)
  density$temp.min.to <- dnorm(theta.samp$temp.min.to, 
                               theta.obs$temp.min.to, 1/ind)
  density$temp.max.to <- dnorm(theta.samp$temp.max.to, 
                               theta.obs$temp.max.to, 1/ind)
  density$altitude.disp <- dnorm(theta.samp$altitude.disp, 
                                 theta.obs$altitude.disp, 16/ind)
  density$temp.min.disp <- dnorm(theta.samp$temp.min.disp, 
                                 theta.obs$temp.min.disp, 1/ind)
  density$est.prob <- dnorm(theta.samp$est.prob, 
                            theta.obs$est.prob, 0.1/ind)
  return(density)
}

##### Model functions #####
post_eps <- function(df.orig,tmin.to,tmax.to,altitude.disp,temp.min.disp) {
  ## Filters out the trajectories that did not start because the T was
  ## outside (tmin.to,tmax.to)
  ids <- df.orig[AgeTraj==0
  ][,start.temp:=AIR_TEMP-273.15
  ][between(start.temp, tmin.to, tmax.to)]$ID2
  df <- df.orig[ID2 %in% ids] 
  cols <- c('Year','Lat','Lon','YMD','ID2', 'AgeTraj')
  # This function finds the first occurence (i.e. minimum AgeTraj) when the
  # conditions are met for an early landing
  if (dim(df)[1] > 0) {
    res.prelim <- df[Elev < altitude.disp | 
                       AIR_TEMP - 273.15 < temp.min.disp | 
                       AgeTraj == 9, 
                     lapply(.SD, min), 
                     by = ID2, .SDcols = 'AgeTraj'
    ][order(ID2)]
    #res <- res.prelim[df, on = .(ID2, AgeTraj), ..cols]
    res <- df[res.prelim, on = .(ID2, AgeTraj), ..cols][AgeTraj > 0]
    
  } else {
    res <- data.frame(Year=integer(),Lat=double(),Lon=double(),YMD=double(),ID2=double())
  }
  res$Year <- res$Year + 2000
  return(as.data.frame(res))
}

## Function to sample parameters given bounds (output is a list of parameter values)
post_parsamp <- function(n) {
  s1 <- round(seq(0.4, 1.2, by = 0.2), 1)
  altitude <- sample(s1, n, replace = TRUE)
  temp.min.to <- rnorm(n, 15, 1)
  temp.max.to <- rnorm(n, 29, 1)
  altitude.disp <- rgamma(n, shape = 7, scale = 7)
  temp.min.disp <- rnorm(n, 15, 1)
  est.prob <- runif(n, 0, 1)
  d <- data.frame(altitude, temp.min.to, temp.max.to,
                  altitude.disp, temp.min.disp, 
                  est.prob)#,temp.max.disp)
  d <- as.list(d)
  return(d)
}

parsamp.density <- function(parsamp) {
  lst <- list(altitude = 1/5,
              temp.min.to = dnorm(parsamp$temp.min.to, 15, 1),
              temp.max.to = dnorm(parsamp$temp.max.to, 29, 1),
              altitude.disp = dgamma(parsamp$altitude.disp,
                                     shape = 7, scale = 7),
              temp.min.disp = dnorm(parsamp$temp.min.disp, 15, 1),
              est.prob = 1)
  return(lst)
}

post_theta.sample <- function(origin, theta, rasters = rlst) {
  altitude <- theta$altitude
  temp.min.to <- theta$temp.min.to
  temp.max.to <- theta$temp.max.to
  altitude.disp <- theta$altitude.disp
  temp.min.disp <- theta$temp.min.disp
  #est.prob <- theta$est.prob
  est.prob <- 1
  
  s.origin <- origin[PBL == altitude]
  
  endpoints.all <- post_eps(s.origin, temp.min.to, temp.max.to, altitude.disp, temp.min.disp)
  n.success <- round(nrow(endpoints.all)*est.prob)
  if (n.success == 0) {
    return(c('accuracy' = 0, 'l2hit' = 0, 'n.ends' = 0))}
  year <- unique(endpoints.all$Year)
  
  endpoints.samp <- endpoints.all[sample(1:nrow(endpoints.all), n.success),]
  ends.all <- as.data.frame(endpoints.samp[,c('Lon', 'Lat')])
  cord.dec <- SpatialPoints(ends.all, proj4string = CRS("+proj=longlat"))
  cord.utm <- as.data.frame(spTransform(cord.dec, CRS('+init=epsg:6623')))
  
  endpoints.samp[,c('x.coord', 'y.coord')] <- cord.utm
  
  ends.df <- endpoints.samp %>% st_as_sf(coords = c('x.coord', 'y.coord'))
  
  end.rast <- rasterize(ends.df, rast, rep(1, nrow(ends.df)), max, na.rm = TRUE)
  values(end.rast) <- sapply(values(end.rast), function(x) {ifelse(is.na(x), 0, x)})
  
  obs <- rasters[[as.character(year)]]
  tot <- sum(values(obs), na.rm = TRUE)
  obs.bin <- obs
  values(obs.bin) <- sapply(values(obs.bin), function(x) {
    ifelse(x == 0, 0, x/x)
  })
  df <- data.frame('obs' = values(obs.bin),
                   'pred' = values(end.rast))
  df.nao <- na.omit(df)
  acc <- length(which(df.nao$obs == df.nao$pred))/nrow(df.nao)
  
  htwt <- sum(na.omit(values(obs))[which(df.nao$pred == 1)])
  htwt.p <- htwt/tot
  return(c('accuracy' = acc, 'l2hit' = htwt.p, 'n.ends' = n.success))
}

sample.n <- function(n = NULL, data = NULL, quant1 = 0.25, quant2 = 0.25, rast = rlst) {
  if (!is.null(data)) {
    n <- nrow(data)
    ags <- aggregate(data = data, cbind(l2.ep, acc.ep) ~ year, unique)
    
    l2.ep <- ags$l2.ep
    names(l2.ep) <- as.character(ags$year)
    
    acc.ep <- ags$acc.ep
    names(acc.ep) <- as.character(ags$year)
    
    rep.lst <- lapply(1:n, function(i) {
      c.l2 <- 0
      c.acc <- 0
      y <- as.character(data$year[i])
      while(c.l2 < l2.ep[[y]] | c.acc < acc.ep[[y]]) {
        s.date <- sample(all.dates, 1)
        #th.row <- results[sample.discrete(1, results$wvec),params]
        th.row <- data[sample.discrete(1, data$wvec),params]
        tdf <- tkernel.sample(th.row, index)
        #ps <- append(th.row, list('date' = s.date))
        ps <- append(tdf, list('date' = s.date))
        dat <- rs.all.cut2[YMD == s.date]
        s.year <- as.character(year(s.date))
        pts <- as.list(post_theta.sample(dat, ps, rasters = rast))
        t.res <- append(ps, pts)
        c.l2 <- t.res$l2hit
        c.acc <- t.res$accuracy
        y <- as.character(year(t.res$date))
      }
      return(t.res)
    })
  }
  else {
    theta.df <- bind_cols(post_parsamp(n))
    rep.lst <- lapply(1:n, function(i) {
      s.date <- sample(all.dates, 1)
      ps <- append(theta.df[i,], list('date' = s.date))
      dat <- rs.all.cut2[YMD == s.date]
      pts <- as.list(post_theta.sample(dat, ps, rasters = rast))
      t.res <- append(ps, pts)
      return(t.res)
    })
  }
  results <- bind_rows(rep.lst)
  results$year <- year(results$date)
  res.strat <- dfapply(unique(results$year), function(y) {
    sub <- subset(results, year == y)
    l2hit <- sub$l2hit
    accuracy <- sub$accuracy
    
    l2.ep <- quantile(l2hit, quant1)
    acc.ep <- quantile(accuracy, quant2)
    
    sub$l2.ep <- l2.ep
    sub$acc.ep <- acc.ep
    
    return(sub)
  })
  return(res.strat)
}

post_next.sample <- function(results, q1 = 0.25, q2 = 0.25, index, rst = rlst) {
  nt.df <- sample.n(data = results, quant1 = q1, quant2 = q2, rast = rst)
  nt.df$wvec <- sapply(1:nrow(nt.df), function(i) {
    tsamp <- nt.df[i,]
    denom.vec <- sapply(1:nrow(nt.df), function(j) {
      tobs <- nt.df[j,]
      kd.lst <- tkernel.density(tobs, tsamp, index)
      kd <- prod(unlist(kd.lst))
      w <- results$wvec[j]
      return(w*kd)
    })
    denom <- sum(denom.vec)
    prior.lst <- parsamp.density(tsamp)
    prior <- prod(unlist(prior.lst))
    wnew <- prior/denom
    return(wnew)
  }) 
  return(nt.df)
}
