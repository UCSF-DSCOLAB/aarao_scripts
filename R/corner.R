corner <- function(x, loc='topleft', nrow=5, ncol=5){
    if (!loc %in% c('topleft', 'topright', 'bottomleft', 'bottomright')) {
        print(paste0('ERROR: loc must be one of (', 
                     paste(c('topleft', 'topright', 'bottomleft', 'bottomright'),
                           collapse=', '),
                     ')'))
        return()
    }
    if (nrow < 0 || nrow > dim(x)[1]) {
        print(paste0('WARNING: nrow must be inthe range [1,', dim(x)[1], '], Using nrow=', dim(x)[1]))
        nrow=dim(x)[1]
    }
    if (ncol < 0 || ncol > dim(x)[2]) {
        print(paste0('WARNING: ncol must be inthe range [1,', dim(x)[2], '], Using ncol=', dim(x)[2]))
        ncol=dim(x)[2]
    }

    if (loc == 'topleft') {
        rows <- c(1:nrow)
        cols <- c(1:ncol)
    } else if (loc == 'topright') {
        rows <- c(1:nrow)
        cols <- c((dim(x)[2]-ncol+1):dim(x)[2])
    } else if (loc == 'bottomleft') {
        rows <- c((dim(x)[1]-nrow+1):dim(x)[1])
        cols <- c(1:ncol)
    } else {
        rows <- c((dim(x)[1]-nrow+1):dim(x)[1])
        cols <- c((dim(x)[2]-ncol+1):dim(x)[2])
    }
    x[rows, cols]
}