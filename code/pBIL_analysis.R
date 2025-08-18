source("code/pBIL.R")

run.lst <- lapply(2013:2018, function(yr) {
  path <- paste0('code/output/')
  name <- paste0('long_run_', yr, '.csv')
  data <- read.csv(paste0(path, name))
  data$year <- yr
  return(data)
})
run.df <- bind_rows(run.lst)