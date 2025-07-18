chrnum <- function(...) {
  args <- list(...)
  x <- args[[1]]
  
  if (length(args) == 1 || args[[2]] %in% c('hg19', 'hg38')) {
    a <- character(0)
    for (i in seq_along(x)) {
      if (x[i] <= 22) {
        a <- c(a, as.character(x[i]))
      } else if (x[i] == 23) {
        a <- c(a, 'X')
      } else if (x[i] == 24) {
        a <- c(a, 'Y')
      }
    }
  } else if (args[[2]] == 'mm10') {
    a <- character(0)
    for (i in seq_along(x)) {
      if (x[i] <= 19) {
        a <- c(a, as.character(x[i]))
      } else if (x[i] == 20) {
        a <- c(a, 'X')
      } else if (x[i] == 21) {
        a <- c(a, 'Y')
      } else if (x[i] == 22) {
        a <- c(a, 'M')
      }
    }
  } else if (args[[2]] == 'dm6') {
    chr_names <- c('2L', '2R', '3L', '3R', '4', 'X', 'Y')
    a <- character(0)
    for (i in seq_along(x)) {
      a <- c(a, chr_names[x[i]])
    }
  }
  
  return(a)
}
