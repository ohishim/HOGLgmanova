# internal functions

norm2 <- function(x){crossprod(x) %>% sqrt}

mprod3 <- function(A, B, C, size=FALSE, n=NULL, m=NULL, k=NULL, p=NULL){

  if(!size){n <- nrow(A); m <- ncol(A); k <- ncol(B); p <- ncol(C)}

  if(n*m*k + n*k*p < m*k*p + n*m*p)
  {
    return((A %*% B) %*% C)
  } else
  {
    return(A %*% (B %*% C))
  }
}

source("inst/defVari.R")
source("inst/subfunctions.R")

usethis::use_data(
  defVari, norm2, mprod3,
  lamCand, GMANOVA.HOGL.fit.MM, GMANOVA.HOGL.fit.BCDM, cw.update,
  internal=T, overwrite = T
)
