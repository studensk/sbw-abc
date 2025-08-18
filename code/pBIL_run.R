library(parallel)

### Parallel years ###
cl <- makeCluster(6)
clusterEvalQ(cl, {
  source('code/pBIL.R')
})
mcmc.lst <- parLapply(cl, 2013:2018, function(y) {
  iter <- mcmc.iter(theta.start)
  theta.prev <- iter$theta
  year(theta.prev$s.date) <- y
  theta.df <- theta.prev
  ll.prev <- ll.fun(iter)
  theta.df$ll <- ll.prev
  theta.df$result <- 'start'
  theta.df$time <- 0
  out.file <- paste0('code/output/long_run_', y, '.csv')
  write.csv(theta.df, out.file, row.names = FALSE)
  n <- 50000
  #n <- 3
  for (i in 1:n) {
    start.time <- proc.time()
    iter <- mcmc.iter(theta.prev)
    theta <- iter$theta
    ll <- ll.fun(iter)
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
    end.time <- proc.time() - start.time
    theta$time <- end.time[3]
    theta$s.date <- as.character(theta$s.date)
    theta.line <- paste0(unlist(theta), collapse = ',')
    theta$s.date <- as.Date(theta$s.date)
    write_lines(theta.line, out.file, append = TRUE)
    theta.df <- bind_rows(theta.df, theta)
  }
  theta.df$year <- y
  return(theta.df)
})
stopCluster(cl)
theta.full <- bind_rows(mcmc.lst)