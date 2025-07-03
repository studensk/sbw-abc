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
#library(ClustGeo)
library(dplyr)
#library(caret)
library(data.table)
library(tictoc)
library(arrow)

## Read in and clean observation data; EPSG 6623 is the mapping for measurements in meters in Quebec
# Allows us to calculate distances in meters
l2 <- read_sf('data/L2.dbf') %>%
  st_transform(crs = '+init=epsg:6623')

#x.mat <- l2[,grep('X201', names(l2))]
x.mat <- l2[,grep('^201', names(l2))]
x.mat <- x.mat[,2:ncol(x.mat)] 
imm.mat <- l2[,grep('Imm_201', names(l2))]
all.mat <-  x.mat*imm.mat
rast.df <- all.mat %>%
  st_as_sf() %>%
  st_drop_geometry() %>%
  mutate(geometry = l2$geometry) %>%
  st_as_sf(crs = st_crs(l2)) #%>%
  # st_buffer(dist = 5000)


## Read in start point rasters
obs <- raster::stack("data/sbw_defol_stack.grd")

## Define years, boundary layer levels and initial parameter boundaries
years <- 2007:2017
pbls <- c(0.4, 0.6, 0.8, 1, 1.2)
lower <- c(0.4, 13, 27, 0, 13)
upper <- c(1.2, 17, 31, 100, 17)
param.names <- c('altitude', 'temp.min.to', 'temp.max.to', 
                 'altitude.disp', 'temp.min.disp')



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

# tkernel.density <- function(theta.obs, theta.samp, ind) {
#   density <- theta.samp[params]
#   density$altitude <- dnorm(theta.samp$altitude,
#                             theta.obs$altitude, 0.2/ind)
#   density$temp.min.to <- dnorm(theta.samp$temp.min.to, 
#                                theta.obs$temp.min.to, 1/ind)
#   density$temp.max.to <- dnorm(theta.samp$temp.max.to, 
#                                theta.obs$temp.max.to, 1/ind)
#   density$altitude.disp <- dnorm(theta.samp$altitude.disp, 
#                                  theta.obs$altitude.disp, 16/ind)
#   density$temp.min.disp <- dnorm(theta.samp$temp.min.disp, 
#                                  theta.obs$temp.min.disp, 1/ind)
#   density$est.prob <- dnorm(theta.samp$est.prob, 
#                             theta.obs$est.prob, 0.1/ind)
#   return(density)
# }

tkernel.density <- function(theta.obs, theta.samp, ind) {
  sd.vec <- c('altitude' = 0.2,
              'temp.min.to' = 1,
              'temp.max.to' = 1,
              'altitude.disp' = 16,
              'temp.min.disp' = 1,
              'est.prob' = 1)
  density <- dnorm(theta.obs, theta.samp, sd.vec/ind)
  return(density)
}

##### Model functions #####
post_eps <- function(df.orig,tmin.to,tmax.to,altitude.disp,temp.min.disp) {
  ## Filters out the trajectories that did not start because the T was
  ## outside (tmin.to,tmax.to)
  
  cols <- c('Year','Lat','Lon','YMD','ID2', 'AgeTraj')
  
  df1 <- df.orig |> 
    mutate(start.temp = AIR_TEMP - 273.15) |> 
    filter(AgeTraj == 0 &
             between(start.temp, tmin.to, tmax.to)) |>
    select(ID2) |>
    left_join(df.orig) |>
    select(all_of(c(cols, 'Elev', 'AIR_TEMP'))) |>
    arrange(ID2, AgeTraj) |>
    mutate(sub.elev = Elev < altitude.disp,
           sub.temp  = AIR_TEMP - 273.15 < temp.min.disp)
  
  df2 <- df1 |>
    group_by(ID2) |>
    collect() |>
    summarize(at1 = min(c(which(sub.elev) - 1, 9)),
              at2 = min(c(which(sub.temp) - 1, 9)),
              AgeTraj = min(at1, at2)) |>
    select(ID2, AgeTraj) |>
    arrow_table(schema = schema(select(df1, c(ID2, AgeTraj))))
  
  res <- df1 |>
    select(all_of(cols)) |>
    right_join(df2) |>
    filter(AgeTraj > 0) |>
    mutate(Year = Year + 2000)
  
  # else {
  #   res <- data.frame(Year=integer(),Lat=double(),Lon=double(),YMD=double(),
  #                     ID2=double())
  # }
  return(res)
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

post_theta.sample <- function(origin, theta, rast.df2) {
  altitude <- theta$altitude
  temp.min.to <- theta$temp.min.to
  temp.max.to <- theta$temp.max.to
  altitude.disp <- theta$altitude.disp
  temp.min.disp <- theta$temp.min.disp
  endpoints.all <- post_eps(origin, temp.min.to, temp.max.to,
                            altitude.disp, temp.min.disp)
  
  year <- year(theta$date)
  # y.rast <- na.omit(rast.df2[,as.character(year)])
  # values <- y.rast[[as.character(year)]]
  values <- rast.df2[[as.character(year)]]
  values.bin <- values/pmax(values, 0.001)
  
  endpoints.df <- endpoints.all |>
    collect() |>
    st_as_sf(coords = c('Lon', 'Lat'), crs = "+proj=longlat") |>
    st_transform(crs = st_crs(rast.df2))
  
  intersection <- st_intersects(rast.df2, endpoints.df)
  ints <- sapply(intersection, function(x) {ifelse(length(x) == 0, 0, 1)})
  
  acc <- length(which(values.bin == ints))/length(values.bin)
  l2hit <- sum(values*ints)/sum(values)
  n.ends <- nrow(endpoints.df)
  
  return(c('accuracy' = acc, 'l2hit' = l2hit, 'n.ends' = n.ends))
}

sample.n <- function(n = NULL, data = NULL, 
                     quant1 = 0.25, quant2 = 0.25, 
                     clusters = 50) {
  if (!is.null(data)) {
    n <- nrow(data)
    ags <- aggregate(data = data, cbind(l2.ep, acc.ep) ~ year, unique)
    
    l2.ep <- ags$l2.ep
    names(l2.ep) <- as.character(ags$year)
    
    acc.ep <- ags$acc.ep
    names(acc.ep) <- as.character(ags$year)
    
    cl <- makeCluster(clusters)
    clusterEvalQ(cl, {
      library(sf)
      library(tidyverse)
      library(arrow)
    })
    clusterExport(cl, c('all.dates', 'post_theta.sample', 
                        'post_eps', 'rast.df', 'sample.discrete',
                        'tkernel.sample', 'data', 
                        'l2.ep', 'acc.ep', 'n', 'params'),
                  envir = environment())
    rep.lst <- parLapply(cl, 1:n, function(i) {
      
      c.l2 <- 0
      c.acc <- 0
      y <- as.character(data$year[i])
      while(c.l2 < l2.ep[[y]] | c.acc < acc.ep[[y]]) {
        date.orig <- data$date[i]
        #s.date <- sample(all.dates, 1)
        m <- 8
        d <- 6
        while (m == 8 & d == 6) {
          #s.date <- sample(all.dates, 1)
          off <- round(rnorm(1, 0, 1))
          s.date <- as.Date(date.orig) + off
          #ymd <- as.Date(s.date)
          m <- month(s.date)
          d <- day(s.date)
        }
        th.row <- data[sample.discrete(1, data$wvec),params]
        tdf <- tkernel.sample(th.row, index)
        ps <- append(tdf, list('date' = s.date))
        dat <- open_dataset('data/archive/res_simul_cut.parquet') |>
          filter(YMD == s.date & PBL == ps$altitude)
        rast2 <- rast.df |>
          dplyr::select(all_of(y)) |>
          na.omit() |>
          st_buffer(dist = 5000)
        #dat <- rs.all |> filter(YMD == s.date)
        #s.year <- as.character(year(s.date))
        pts <- as.list(post_theta.sample(dat, ps, rast.df2 = rast2))
        t.res <- append(ps, pts)
        c.l2 <- t.res$l2hit
        c.acc <- t.res$accuracy
        y <- as.character(year(t.res$date))
      }
      return(t.res)
    })
    stopCluster(cl)
  }
  else {
    # theta.df <- bind_cols(post_parsamp(n))
    # s.dates <- sample(all.dates, n)
    cl <- makeCluster(clusters)
    clusterEvalQ(cl, {
      library(sf)
      library(tidyverse)
      library(arrow)
    })
    clusterExport(cl, c('post_parsamp', 'all.dates', 'post_theta.sample', 
                        'post_eps', 'rast.df'))
    rep.lst <- parLapply(cl, 1:n, function(x) {
      m <- 8
      d <- 6
      while (m == 8 & d == 6) {
        s.date <- sample(all.dates, 1)
        ymd <- as.Date(s.date)
        m <- month(ymd)
        d <- day(ymd)
      }
      theta.df <- bind_cols(post_parsamp(1)) %>%
        mutate('date' = s.date)
      #write.csv(theta.df, paste0('code/output/theta/', x, '.csv'))
      theta.df <- select(theta.df, -date)
      yr <- year(s.date)
      rast2 <- rast.df |>
        dplyr::select(all_of(as.character(yr))) |>
        na.omit() |>
        st_buffer(dist = 5000)
      ps <- append(theta.df, list('date' = s.date))
      dat <- open_dataset('data/archive/res_simul_cut.parquet') |>
        filter(YMD == s.date & PBL == ps$altitude)
      pts <- as.list(post_theta.sample(dat, ps, rast.df2 = rast2))
      t.res <- append(ps, pts)
      #write.csv(as.data.frame(t.res), paste0('code/output/results/', x, '.csv'))
      return(t.res)
    })
    stopCluster(cl)
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

sd.vec <- c('altitude' = 0.2,
            'temp.min.to' = 1,
            'temp.max.to' = 1,
            'altitude.disp' = 16,
            'temp.min.disp' = 1,
            'est.prob' = 1)

post_next.sample <- function(results, q1 = 0.25, q2 = 0.25, index) {
  nt.df <- sample.n(data = results, quant1 = q1, quant2 = q2)
  # nt.df$wvec <- sapply(1:nrow(nt.df), function(i) {
  #   tsamp <- nt.df[i,]
  #   denom.vec <- sapply(1:nrow(nt.df), function(j) {
  #     tobs <- nt.df[j,]
  #     kd.lst <- tkernel.density(tobs, tsamp, index)
  #     kd <- prod(unlist(kd.lst))
  #     w <- results$wvec[j]
  #     return(w*kd)
  #   })
  #   denom <- sum(denom.vec)
  #   prior.lst <- parsamp.density(tsamp)
  #   prior <- prod(unlist(prior.lst))
  #   wnew <- prior/denom
  #   return(wnew)
  # }) 
  ind.vec <- 1:nrow(nt.df)
  eg.nt.df <- expand.grid('ind.obs' = ind.vec, 'ind.samp' = ind.vec)
  
  prior.prod <- sapply(ind.vec, function(ind) {
    tsamp <- nt.df[ind,]
    prior.lst <- parsamp.density(tsamp)
    prior <- prod(unlist(prior.lst))
    return(prior)
  })
  
  combo.df <- expand.grid('prior' = prior.prod, 
                          'wv' = results$wvec) %>%
    cbind(eg.nt.df)
  s <- sapply(params, function(p) {
    vec.orig <- nt.df[,p]
    eg <- expand.grid(vec.orig, vec.orig)
    k <- dnorm(eg[,2], eg[,1], sd.vec[p]/index)
    return(k)
  })
  mult.df <- combo.df |>
    mutate(k = apply(s, 1, prod),
           denom.vec = wv*k) |>
    group_by(ind.obs, prior) |>
    summarize(denom.val = sum(denom.vec)) |>
    mutate(wvec = prior/denom.val) |>
    ungroup() |>
    select(wvec)
  nt.df.new <- cbind(nt.df, mult.df)
  return(nt.df.new)
}


