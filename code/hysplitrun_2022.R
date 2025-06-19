source('code/functions.R')

## Read Hysplit trajectory data
rs.all <- read.csv('data/res_simul_all.csv', header = T)
rs.all.cut1 <- subset(rs.all, Year > 12)
rs.all.cut2 <- subset(rs.all.cut1, Month < 8 | Day <= 6)
rm(rs.all, rs.all.cut1)
rs.all.cut2 <- data.table(rs.all.cut2)

all.dates <- unique(rs.all.cut2[,YMD])

##### Run ABC SMC #####

ptm.full.start <- proc.time()
ptm <- proc.time()
results <- sample.n(2500) 
results$wvec <- 1/nrow(results)
results$index <- 1
print(proc.time() - ptm)
ptm <- proc.time()
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
write.csv(comp.df, 'code/output/new_results_test.csv', row.names = FALSE)
