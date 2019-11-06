# read_1 -----------------------------------------------------------------------
read_1 <- function(file)
{
  utils::read.csv2(
    file, stringsAsFactors = FALSE
  )
}

# read_2 -----------------------------------------------------------------------
read_2 <- function(file, ...)
{
  utils::read.table(
    file, header = TRUE, sep = ";", dec = ".", stringsAsFactors = FALSE, ...
  )
}
