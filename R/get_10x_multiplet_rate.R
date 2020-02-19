library(assertthat)

args = commandArgs(trailingOnly=TRUE)
vals <- as.numeric(args)

x = assert_that(all(args[c(1:length(args))] > 0), msg="args must be greater than 0")

df <- read.table('~/important_data/10x_multiplet_rate.tsv', header=T)

y <- lm(multiplet_rate_pct ~ num_cells_loaded, data=df)
z <- lm(num_cells_recovered ~ num_cells_loaded, data=df)

pred <- data.frame(multiplet_rate=predict(y, new=data.frame(num_cells_loaded=vals)),
                   recovered_cells=predict(z, new=data.frame(num_cells_loaded=vals)))
rownames(pred) <- vals
print(pred)

