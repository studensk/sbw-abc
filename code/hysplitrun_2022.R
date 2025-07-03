source('code/functions.R')

## Read Hysplit trajectory data
## This file is the parquet version of rs.all.cut2
rs.all <- open_dataset('data/archive/res_simul_cut.parquet')

all.dates <- rs.all |>
  select(YMD) |>
  unique() |>
  collect() |>
  unlist(use.names = FALSE)
rm(rs.all)

##### Run ABC SMC #####

tic()
results <- sample.n(2500) 
results$wvec <- 1/nrow(results)
results$index <- 1
toc()
results <- subset(results, l2hit > l2.ep & accuracy > acc.ep)
res.lst <- list(results)
for (i in 1:5) {
  print(proc.time() - ptm)
  ptm <- proc.time()
  res <- res.lst[[i]]
  results <- post_next.sample(res, index = i, q1 = 0.25, q2 = 0.25)
  results$index <- i + 1
  res.lst <- append(res.lst, list(results))
}

ptm.full.end <- proc.time()
ptm.full <- ptm.full.end - ptm.full.start

comp.df <- bind_rows(res.lst)
write.csv(comp.df, 'code/output/new2_results_test.csv', row.names = FALSE)
