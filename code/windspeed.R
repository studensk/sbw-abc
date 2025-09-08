library(sf)
library(tidyverse)
library(tictoc)
library(arrow) 
library(parallel)

l2.orig <- read_sf('data/L2.dbf') 
ll.crs <- st_crs(l2.orig)
dec.crs <- '+init=epsg:6623'
rs.all <- open_dataset('data/archive/res_simul_cut.parquet')

all.dates <- rs.all |>
  select(YMD) |>
  unique() |>
  collect() |>
  unlist()

tic()
cl <- makeCluster(70)
clusterExport(cl, c('ll.crs', 'dec.crs', 'all.dates'))
clusterEvalQ(cl, {
  library(sf)
  library(tidyverse)
  library(arrow) 
})
toc()
tic()
rs.lst <- parLapply(cl, all.dates, function(dt) {
  rs.all <- open_dataset('data/archive/res_simul_cut.parquet')
  rs.samp <- rs.all |>
    filter(YMD %in% dt) |> 
    collect()
  
  dec.geo <- rs.samp |>
    st_as_sf(coords = c('Lon', 'Lat'), crs = ll.crs) |>
    st_transform(crs = dec.crs) |>
    st_coordinates() |>
    as.data.frame() |>
    cbind(rs.samp)
  
  rs.sf <- dec.geo |>
    group_by(ID2) |>
    mutate(x.diff = c(diff(X), NA),
           y.diff = c(diff(Y), NA),
           e.diff = c(diff(Elev), NA)) |>
    ungroup() |>
    mutate(t.windspeed = sqrt(x.diff^2 + y.diff^2 + e.diff^2)/1000,
           h.windspeed = sqrt(x.diff^2 + y.diff^2)/1000,
           v.windspeed = e.diff/1000) |>
    select(-c(x.diff, y.diff, e.diff, X, Y))
  return(rs.sf)
})
stopCluster(cl)
toc()

rs.df <- bind_rows(rs.lst)
write_parquet(rs.df, 'data/res_simul_all_ws.parquet')

