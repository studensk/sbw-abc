source('code/pBIL.R')
library(ramcmc)
library(MASS)
library(raster)
library(sf)
library(tictoc)
library(arrow)
library(tidyverse)

scale.default <- c('c.temp' = 1,
                   'r.temp' = 1, 
                   'altitude' = 0.75, 
                   'h.wind' = 0.75,
                   'd.altitude' = 0.04,
                   'd.tmin' = 1,
                   'alpha' = 0.5,
                   'l.lambda' = 0.05,
                   #'s.date' = 2)
                   's.date' = 4)

proposal <- function(theta, iter.block = c('climatic', 'date', 'count'),
                     scale = scale.default) {
  date <- theta$s.date
  y <- year(date)
  date.min <- as.numeric(as.Date(paste0(y, '-06-15')))
  date.max <- as.numeric(as.Date(paste0(y, '-07-18')))
  scale.df <- data.frame('variable' = c('c.temp', 'r.temp', 'altitude', 
                                        'h.wind', 'd.altitude', 'd.tmin', 
                                        'alpha', 'l.lambda', 's.date'),
                         'scale' = scale,
                         'log.trans' = c(FALSE, FALSE, FALSE, 
                                   TRUE, TRUE, FALSE,
                                   TRUE, TRUE, FALSE),
                         'block' = c(rep('climatic', 6), rep('count', 2), 
                                     'date')) |>
    mutate(theta = as.numeric(theta),
           theta.trans = log.trans*log(theta) + (1-log.trans)*theta,
           norm.add = rnorm(length(scale), 0, scale),
           theta.new = round(log.trans*pmin(1e+100, 
                                            exp(theta.trans + norm.add)) +
             (1 - log.trans)*(theta.trans + norm.add), 3),
           theta.new = ifelse(block == 'date', round(theta.new, 0),
                              ifelse(variable == 'altitude',
                                     round(theta.new/0.2)*0.2, theta.new)),
           theta.new = ifelse(block %in% iter.block, theta.new, theta))
  
  
  theta.new.lst <- as.list(scale.df$theta.new)
  names(theta.new.lst) <- scale.df$variable
  
  date.new <- theta.new.lst$s.date
  alt.new <- theta.new.lst$altitude

  while(!between(date.new, date.min, date.max)) {
    date.new <- with(subset(scale.df, variable == 's.date'), {
      add <- rnorm(1, 0, scale)
      new <- round(date.new + add)
      return(new)
    })
  }
  
  while(!between(alt.new, 0.2, 1.2)) {
    alt.new <- with(subset(scale.df, variable == 'altitude'), {
      add <- rnorm(1, 0, scale)
      new <- round((alt.new + add)/0.2)*0.2
      return(new)
    })
  }
  
  theta.new.lst$s.date <- as.Date(date.new)
  theta.new.lst$altitude <- alt.new
  
  if (!is.null(block))
  
  return(theta.new.lst)
}

accept.reject <- function(theta.prev, ll.prev, theta, ll) {
  rterm.num <- ll + prior.density(theta) + 
    proposal.density(theta.prev, theta)
  rterm.den <- ll.prev + prior.density(theta.prev) +
    proposal.density(theta, theta.prev)
  rterm.full <- pmin(1, exp(rterm.num - rterm.den))
  if (is.na(rterm.full)) {rterm.full <- 0}
  u <- runif(1)
  if (u < rterm.full) {
    theta.prev <- theta
    ll.prev <- ll
    theta$ll <- ll
    theta$result <- 'accept'
  }
  else {
    theta$ll <- ll
    theta$result <- 'reject'
  }
  return(theta)
}

theta.start.iter <- mcmc.iter(theta.start, new = FALSE)
theta.start.df <- theta.start.iter$theta
theta.start.df$ll <- ll.fun(theta.start.iter)
theta.start.df$result <- 'start'
  

one.step <- function(theta.df, block = c('climatic', 'date', 'count')) {
  ll <- theta.df$ll
  theta.sub <- select(theta.df, -c(ll, result))
  iter <- mcmc.iter(theta.sub, iter.block = block)
  iter.ll <- ll.fun(iter)
  theta <- accept.reject(theta.df, ll, iter$theta, iter.ll)
  theta.df <- bind_rows(theta.df, as.data.frame(theta))
  return(theta.df)
}

three.steps <- function(theta.df) {
  blocks <- c('date', 'climatic', 'count')
  step <- 1
  for (b in blocks) {
    print(b)
    # theta <- theta.df[step,]
    # ll <- theta$ll
    # theta.sub <- select(theta, -c(ll, result))
    # iter <- mcmc.iter(theta.sub, iter.block = b)
    # iter.ll <- ll.fun(iter)
    # theta <- accept.reject(theta, ll, iter$theta, iter.ll)
    # theta.df <- bind_rows(theta.df, as.data.frame(theta))
    theta.df <- one.step(theta.df, block = b)
    if (theta$result == 'accept') {step <- step + 1}
  }
  return(theta.df[2:4,])
}

theta.start.df$result <- 'accept'
theta.df <- theta.start.df
theta.full <- theta.df

# out.file <- 'code/output/mh_gibbs_block_run.csv'
# write.csv(theta.df, out.file, row.names = FALSE)

for (i in 1:1920) {
  print(i)
  theta.new <- three.steps(theta.df)
  theta.full <- bind_rows(theta.full, theta.new)
  theta.df <- theta.full[max(which(theta.full$result == 'accept')),]
}


out.file <- 'code/output/mh_gibbs_block_run.csv'
write.csv(theta.df, out.file, row.names = FALSE)
theta$s.date <- as.character(theta$s.date)
theta.line <- paste0(unlist(theta), collapse = ',')
theta$s.date <- as.Date(theta$s.date)
write_lines(theta.line, out.file, append = TRUE)