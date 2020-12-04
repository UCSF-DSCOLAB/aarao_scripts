suppressPackageStartupMessages({
  if (!'package:assertthat' %in% search()) library(assertthat)
})

#https://kb.10xgenomics.com/hc/en-us/articles/360001378811
df <- data.frame(multiplet_rate_pct=c(0.40, 0.80, 1.60, 2.30, 3.10, 3.90, 4.60, 5.40, 6.10, 6.90, 7.60), 
                 num_cells_loaded=c(870, 1700, 3500, 5300, 7000, 8700, 10500, 12200, 14000, 15700, 17400),
                 num_cells_recovered=c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000))


MULTIPLET_MODEL <- lm(multiplet_rate_pct ~ num_cells_loaded, data=df)
RECOVERY_MODEL <- lm(num_cells_recovered ~ num_cells_loaded, data=df)


get_10x_multiplet_rate <- function(x, precision=3, fraction=TRUE) {
    assert_that(is.numeric(x) && length(x) == 1)
    pred <- predict(MULTIPLET_MODEL, new=data.frame(num_cells_loaded=x))
    if (fraction){
        round(unname(pred)/100, precision)
    } else {
        round(unname(pred), precision)
    }    
}

get_10x_recovery_rate <- function(x) {
    assert_that(is.numeric(x) && length(x) == 1)
    pred <- predict(RECOVERY_MODEL, new=data.frame(num_cells_loaded=x))
    unname(pred)
}

# https://stackoverflow.com/a/47932989/3542991
# Equivalent to `if __name__ == '__main__':` in python
if (sys.nframe() == 0L) {
    args = commandArgs(trailingOnly=TRUE)
    if (length(args) == 0 | "-h" %in% args | "--help" %in% args){
        print('USAGE: Rscript get_10x_multiplet_rate.R cells_loaded1 cells_loaded2 cells ... cell_loadedn')
        print('where cells_loaded1,2..n are the number of cells loaded in the well.')
        quit("no")
    }

    vals <- sort(as.numeric(args))

    suppmsg = assert_that(all(args[c(1:length(vals))] > 0), msg="args must be greater than 0")

    results = data.frame(row.names=vals,
                         multiplet_rate=sapply(vals, function(x){ get_10x_multiplet_rate(x, fraction=FALSE) }),
                         recovery=sapply(vals, function(x){ get_10x_recovery_rate(x) }))
    results
}