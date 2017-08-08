# This script checks which variables can be repeated...
listing_duplicated_columns <- function(x) {

  # Fixing factors to strings
  x <- data.frame(lapply(x, function(y)
    if (is.factor(y)) as.character(y) else y
    ), stringsAsFactors = FALSE)

  k <- ncol(x)
  vnames <- colnames(x)
  ans <- NULL
  for (i in 1L:ncol(x))
    for (j in i:ncol(x)) {
      if (i == j) next
      if (all(x[,i] == x[,j], na.rm = TRUE))
        ans <- c(ans, list(c(vnames[c(i,j)])))
    }

  do.call(rbind, ans)
}

