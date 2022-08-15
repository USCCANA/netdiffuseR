x   <- readLines("DESCRIPTION")
ids <- which(grepl("^(Imports:|Suggests:|VignetteBuilder:)", x))
x   <- x[setdiff((ids[1] + 1):(ids[3] - 1), ids[2])]

x   <- gsub("^\\s+|(?<=[[:alnum:]])(\\s|,).*", "", x, perl = TRUE)
x   <- setdiff(x, c("MASS", "methods", "grDevices", "graphics", "stats", "utils"))

xfin <- x
for (i in x) {
  if (require(i))
    xfin <- setdiff(xfin, i)
}

writeLines(xfin, con = "pkgs.txt")
