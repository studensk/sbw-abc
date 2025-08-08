library(parallel)
tic()
cl <- makeCluster(33)
clusterEvalQ(cl, {
  source('code/pBIL.R')
})
toc()
tic()
date.lst <- parLapply(cl, 0:32, function(d) {
  s.date <- as.Date('2013-06-15') + days(d)
  theta.start$s.date <- s.date
  iter <- mcmc.iter(theta.start)
  theta.prev <- iter$theta
  theta.prev$s.date <- s.date
  theta.df <- theta.prev
  ll.prev <- ll.fun(iter)
  theta.df$ll <- ll.prev
  theta.df$result <- 'start'
  times <- c()
  n <- 20
  for (i in 1:n) {
    start.time <- proc.time()
    iter <- mcmc.iter(theta.prev)
    theta <- iter$theta
    theta$s.date <- s.date
    ll <- ll.fun(iter)
    rterm.num <- ll + prior.density(theta) + proposal.density(theta.prev, theta)
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
    end.time <- proc.time() - start.time
    theta.df <- bind_rows(theta.df, theta)
    write.csv(theta.df, paste0('code/output/pbil', d, '.csv'), 
              row.names = FALSE)
    times <- c(times, end.time[3])
  }
  #theta.df$iteration <- 1:(n+1)
  return(theta.df)
})
stopCluster(cl)
toc()

ar.df <- date.df |>
  filter(result != 'start') |>
  count(s.date, result)

ggplot(data = ar.df) +
  geom_bar(aes(x = s.date, y = n, fill = result), stat = 'identity') +
  theme_minimal()

### No parallel ###
source('code/pBIL.R')
iter <- mcmc.iter(theta.start)
theta.prev <- iter$theta
theta.df <- theta.prev
ll.prev <- ll.fun(iter)
theta.df$ll <- ll.prev
theta.df$result <- 'start'
times <- c()
n <- 2000
for (i in 1:n) {
  start.time <- proc.time()
  iter <- mcmc.iter(theta.prev)
  theta <- iter$theta
  ll <- ll.fun(iter)
  rterm.num <- ll + prior.density(theta) + proposal.density(theta.prev, theta)
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
  end.time <- proc.time() - start.time
  theta.df <- bind_rows(theta.df, theta)
  if (i %% 10 == 0) {
    write.csv(theta.df, paste0('code/output/aug7_8_overnight.csv'), 
              row.names = FALSE)
  }
  times <- c(times, end.time[3])
}