library(MASS)
library(raster)
library(sf)
library(tictoc)
library(arrow)
library(tidyverse)

#source('code/functions.R')

st_switch <- function(st.obj, ll = TRUE) {
  coords <- st_coordinates(st.obj) |>
    as.data.frame()
  new.df <- st_drop_geometry(st.obj) |>
    cbind(coords) |>
    rename(Lon = X, Lat = Y)
  return(new.df)
}

## Create arrow for trajectory (HYSPLIT) dataset
rs.all <- open_dataset('data/res_simul_all_ws.parquet')

l2.orig <- read_sf('data/L2.dbf') 
ll.crs <- st_crs(l2.orig)
dec.crs <- '+init=epsg:6623'

## Read in L2 data (end points)
l2 <- l2.orig |>
  st_transform(crs = dec.crs)
x.mat <- l2[,grep('^201', names(l2))]
x.mat <- x.mat[,2:ncol(x.mat)] 
imm.mat <- l2[,grep('Imm_201', names(l2))]
all.mat <-  x.mat*imm.mat
l2.df <- all.mat %>%
  st_as_sf() %>%
  st_drop_geometry() %>%
  mutate(geometry = l2$geometry) %>%
  st_as_sf(crs = st_crs(l2)) 
all.years <- 2013:2018
l2.lst <- lapply(all.years, function(year) {
  sub <- l2.df[,as.character(year)] |>
    na.omit() |>
    as.data.frame() 
  sub.st <- sub |>
    st_as_sf() |>
    st_coordinates() |>
    bind_cols(sub[,1]) |>
    as.data.frame() |>
    mutate(year = year)
  names(sub.st) <- c('x', 'y', 'l2count', 'year')
  return(sub.st)
})
#names(l2.lst) <- as.character(all.years)
l2.geo.df <- bind_rows(l2.lst)

d.lims <- st_bbox(l2.df)[c(1, 3, 2, 4)]

## Read in defoliation data (start points)
obs <- raster::stack("data/sbw_defol_stack.grd") |>
  rasterToPoints()
obs.sf <- st_as_sf(as.data.frame(obs), 
                   coords = c('x', 'y'), crs = dec.crs) |>
  st_transform(crs = crs(l2.df))
obs.df <- obs.sf[,7:ncol(obs.sf)] |>
  pivot_longer(contains('sbw'), 
               names_to = 'year', values_to = 'defoliation') |>
  mutate(year = substr(year, 4, 7))

## theta is a list of parameter values; proposes parameter values for next
##  iteration of MCMC

theta.start <- list('c.temp' = 22,
                    'r.temp' = 14,
                    'altitude' = 0.6,
                    'h.wind' = 25,
                    'd.altitude' = 50,
                    'd.tmin' = 15,
                    'alpha' = 10,
                    'l.lambda' = 2,
                    's.date' = as.Date('2013-07-01'))

prior.density <- function(theta) {
  p <- with(theta, {
    #altitude.p = log(1/5)
    c.temp.p = dnorm(c.temp, 22, 0.7, log = TRUE)
    r.temp.p = dnorm(r.temp, 14, 1.5, log = TRUE)
    d.altitude.p = dgamma(d.altitude, shape = 7, scale = 7, log = TRUE)
    d.tmin.p = dnorm(d.tmin, 15, 1, log = TRUE)
    h.wind.p = dlnorm(h.wind, log(10), 1, log = TRUE)
    l.lambda.p = dlnorm(l.lambda, log(1.5), 0.5, log = TRUE)
    alpha.p = dlnorm(alpha, log(5), 0.75, log = TRUE)
    #s.date.p = log(1/42)
    sum(c.temp.p, r.temp.p, d.altitude.p, d.tmin.p, h.wind.p, 
        l.lambda.p, alpha.p)
  }) 
  return(p)
}

proposal <- function(theta) {
  theta.new <- with(theta, {

    c.temp.new <- rnorm(1, c.temp, 1)
    r.temp.new <- rnorm(1, r.temp, 1)

    alt.add <- round(rnorm(1, 0, 0.75))*0.2
    altitude.new <- altitude + alt.add
    while(!between(altitude.new, 0.4, 1.2)) {
      alt.add <- round(rnorm(1, 0, 0.75))*0.2
      altitude.new <- altitude + alt.add
    }

    windspeed.new <- rlnorm(1, log(h.wind), 0.05)
    d.altitude.new <- rlnorm(1, log(d.altitude), 0.03)
    d.tmin.new <- rnorm(1, d.tmin, 1)
    
    l.lambda.new <- rlnorm(1, log(l.lambda), 0.05)
    alpha.new <- rlnorm(1, log(alpha), 0.5)

    date.start <- as.Date(paste0(year(s.date), '-06-15'))
    date.end <- as.Date(paste0(year(s.date), '-07-18'))
    date.window <- seq(date.start, date.end, by = 'day')

    date.add <- round(rnorm(1, 0, 2))
    date.new <- s.date + date.add
    while(!between(date.new, date.start, date.end)) {
      date.add <- round(rnorm(1, 0, 2))
      date.new <- date.new + date.add
    }
    list('c.temp' = c.temp.new,
         'r.temp' = r.temp.new,
         'altitude' = altitude.new,
         'h.wind' = windspeed.new,
         'd.altitude' = d.altitude.new,
         'd.tmin' = d.tmin.new,
         'alpha' = alpha.new,
         'l.lambda' = l.lambda.new,
         's.date' = date.new)
  })
  return(theta.new)
}

proposal.density <- function(theta.new, theta.old) {
  names(theta.old) <- paste0(names(theta.old), '.old')
  names(theta.new) <- paste0(names(theta.new), '.new')
  
  theta.lst <- append(as.list(theta.old), as.list(theta.new))
  
  prior <- with(theta.lst, {
    
    c.temp.d <- dnorm(c.temp.new, c.temp.old, 1, log = TRUE)
    r.temp.d <- dnorm(r.temp.new, r.temp.old, 1, log = TRUE)
    
    altitude.d <- dnorm(altitude.new, altitude.old, 0.75, log = TRUE)
    
    windspeed.d <- dlnorm(h.wind.new, log(h.wind.old), 0.05, log = TRUE)
    d.altitude.d <- dlnorm(d.altitude.new, 
                           log(d.altitude.old), 0.03, log = TRUE)
    d.tmin.d <- dnorm(d.tmin.new, d.tmin.old, 1, log = TRUE)
    
    l.lambda.d <- dlnorm(l.lambda.new, 
                         log(l.lambda.old), 0.05, log = TRUE)
    alpha.d <- dlnorm(alpha.new, log(alpha.old), 0.5, log = TRUE)
    
    s.date.diff <- as.numeric(s.date.new - s.date.old)
    
    s.date.d <- dnorm(s.date.diff, 0, 2, log = TRUE)
    sum(c.temp.d, r.temp.d, altitude.d, windspeed.d, d.altitude.d, d.tmin.d,
         l.lambda.d, alpha.d, s.date.d)
  })
  return(prior)
}

get_starts <- function(df, id = FALSE) {
  if (id) {
    vars <- c('Lon', 'Lat', 'ID2')
  }
  else {vars <- c('Lon', 'Lat')}
  df |>
    filter(AgeTraj == 0) |>
    select(all_of(vars)) |>
    unique()
}

post_eps_ws <- function(df.orig, theta) {
  res <- with(theta, {
    cols <- c('Year','Lat','Lon','YMD','ID2', 'AgeTraj')
    
    df1 <- df.orig |> 
      mutate(start.temp = AIR_TEMP - 273.15,
             tmin = c.temp - r.temp/2,
             tmax = c.temp + r.temp/2) |> 
      filter(AgeTraj == 0 &
               between(start.temp, tmin, tmax) &
               h.windspeed >= h.wind) |>
      dplyr::select(ID2) |>
      left_join(df.orig) |>
      dplyr::select(all_of(c(cols, 'Elev', 'AIR_TEMP'))) |>
      arrange(ID2, AgeTraj) |>
      mutate(sub.elev = Elev < d.altitude,
             sub.temp  = AIR_TEMP - 273.15 < d.tmin)
    
    df2 <- df1 |>
      group_by(ID2) |>
      collect() |>
      summarize(at1 = min(c(which(sub.elev) - 1, 9)),
                at2 = min(c(which(sub.temp) - 1, 9)),
                AgeTraj = min(at1, at2)) |>
      dplyr::select(ID2, AgeTraj) |>
      arrow_table(schema = schema(select(df1, c(ID2, AgeTraj))))
    
    df1 |>
      dplyr::select(all_of(cols)) |>
      right_join(df2) |>
      filter(AgeTraj > 0) |>
      mutate(Year = Year + 2000)
  })
  return(res)
}

mcmc.iter <- function(theta, traj.data = rs.all) {
  ns.rows <- 1
  theta <- proposal(theta)
  s.date.orig <- theta$s.date
  year <- year(s.date.orig)
  s.date <- s.date.orig
  ends.full <- list()
  short.year <- year - 2000
  
  traj.data.iter <- traj.data |>
    filter(Year == short.year)
  
  starts <- traj.data.iter |>
    get_starts() |>
    collect()
  
  id.starts <- traj.data.iter |>
    filter(Year == short.year) |>
    get_starts(id = TRUE) |>
    collect()
  
  while(ns.rows > 0 & s.date <= (s.date.orig + 19)) {
    #tic()
    #print(s.date)
    rs.date <- traj.data |>
      filter(YMD == s.date)
    ends.df <- post_eps_ws(rs.date, theta) |>
      collect()
    ends.full <- append(ends.full, list(ends.df))
    ids <- unique(ends.df$ID2)
    starts.iter <- id.starts |>
      filter(ID2 %in% ids)
    starts.iter.ll <- starts.iter |>
      select(Lon, Lat) |>
      unique()
    no.starts <- setdiff(starts, starts.iter.ll)
    #mg <- merge(no.starts, traj.data.iter)
    mg <- left_join(no.starts, id.starts,
                    by = join_by(Lon, Lat))
    ns.ids <- unique(mg$ID2)
    ns.rows <- nrow(no.starts)
    #print(ns.rows)
    s.date <- s.date + days(1)
    starts.iter <- starts.iter |>
      filter(ID2 %in% ns.ids)
    #print(nrow(traj.data.iter))
    #toc()
  }
  ends.full.df <- bind_rows(ends.full) |>
    st_as_sf(coords = c('Lon', 'Lat'), crs = ll.crs) |>
    st_transform(l2.df, crs = dec.crs) |>
    st_coordinates() |>
    as.data.frame()
  return(list('endpoints' = ends.full.df, 'theta' = as.data.frame(theta)))
}

## obs and pred are both lists (or dfs) with elements x and y
kds2d <- function(obs, pred, h = 50000) {
  names(obs) <- tolower(names(obs))
  names(pred) <- tolower(names(pred))
  nobs <- length(obs$x)
  npred <- length(pred$x)
  
  #h <- with(obs, c(bandwidth.nrd(x), bandwidth.nrd(y))/4)
  
  ax <- outer(obs$x, pred$x, "-")/h
  ay <- outer(obs$y, pred$y, "-")/h
  
  k <- dnorm(ax)*dnorm(ay)
  az <- apply(k, 2, sum)/(nobs*h^2)
  
  # az <- c()
  # for (i in 1:npred) {
  #   z <- sum(dnorm(ax[,i])*dnorm(ay[,i]))/(nobs*h^2)
  #   az <- c(az, z)
  # } 
  result <- pred
  result$z <- (az/sum(az))*npred
  return(result)
}

ll.fun <- function(iter, l2 = l2.geo.df) {
  iter.e <- iter$endpoints
  iter.t <- iter$theta
  
  l2.y <- l2 |>
    subset(year == year(iter.t$s.date))
  
  s.d <- kds2d(obs = iter.e, pred = l2.y)$z
  
  s.lambda <- with(iter.t, l.lambda + alpha*s.d)
  ll <- dpois(round(l2.y$l2count), s.lambda, log = TRUE)
  return(sum(ll))
}

