source('code/functions.R')

##### Run ABC SMC #####

tic()
results <- sample.n(2500) 
results$wvec <- 1/nrow(results)
results$index <- 1
toc()
results <- subset(results, l2hit > l2.ep & accuracy > acc.ep)
write.csv(results, 'code/output/results1.csv', row.names = FALSE)
res.lst <- list(results)
for (i in 1:5) {
  tic()
  res <- res.lst[[i]]
  results <- post_next.sample(res, index = i, q1 = 0.25, q2 = 0.25)
  new.ind <- i + 1
  results$index <- new.ind
  write.csv(results, paste0('code/output/results', new.ind, '.csv'), row.names = FALSE)
  res.lst <- append(res.lst, list(results))
  toc()
}

comp.df <- bind_rows(res.lst)
write.csv(comp.df, 'code/output/new2_results_test.csv', row.names = FALSE)
