# Finding the longest common starting substrings in a set of strings
# Source: https://stackoverflow.com/questions/28273716/r-implementation-for-finding-the-longest-common-starting-substrings-in-a-set-of

lcPrefix <- function (x, ignore.case = FALSE) 
{
  x <- as.character(x)
  if (ignore.case) 
    x <- toupper(x)
  nc <- nchar(x, type = "char")
  for (i in 1:min(nc)) {
    ss <- substr(x, 1, i)
    if (any(ss != ss[1])) {
      return(substr(x[1], 1, i - 1))
    }
  }
  substr(x[1], 1, i)
}
