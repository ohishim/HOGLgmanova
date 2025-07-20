#' @title Generation of simulation dataset
#' @description \code{genData} this function generates simulation dataset (v0.1.0)
#'
#' @importFrom expm sqrtm
#' @importFrom ggplot2 aes element_blank element_text geom_line ggplot
#'    scale_color_discrete theme xlab
#' @importFrom magrittr %>% divide_by inset set_colnames
#' @importFrom mvnfast rmvn
#' @importFrom purrr exec map map_dbl
#' @importFrom stats runif
#' @importFrom tidyr all_of pivot_longer
#'
#' @param n the sample size
#' @param p the number of time points
#' @param k the number of explanatory variables satisfying `k >= 5`
#' @param q a dimension of polynomial basis satisfying `q >= 6`
#'    which means a polynomial of degree `q-1`
#' @param rho a parameter of autocorrelation matrix
#' @param SNR SNR
#'
#' @return a list with the following elements:
#' \item{Y}{a matrix of response variables}
#'
#' \item{A}{a matrix of explanatory variables}
#'
#' \item{THETA}{a matrix of true coefficients}
#'
#' \item{X}{a polynomial basis matrix}
#'
#' \item{MU}{the expectation of `Y`}
#'
#' \item{time.points}{a vector of time points}
#'
#' \item{SIGMA}{true covariance matrix}
#'
#' \item{degrees}{true degrees of polynomials}
#'
#' \item{trend}{a figure of varying coefficients}
#'
#' @export
#' @examples
#' #genData(n, p, k, q)

genData <- function(n, p, k, q, rho=0.5, SNR=1){

  k. <- 5

  if(k < k.){stop('`k` must be `k >= 5`')}
  if(q < 6){stop('`q` must be `q >= 6`')}

  time.points <- seq(-1, 1, length=p)

  THETA <- list(
    c(0.5, -3),
    c(-2, 1, 4),
    2*c(1, -2, -1, 3),
    3*c(0.5, -1, -4, 1, 4),
    c(-0.5, -1, 1, 15, -1, -12)
  ) %>% map(~{
    d <- length(.x)
    return(numeric(q) %>% inset(1:d, .x) %>% rev)
  }) %>% exec(rbind, !!!.) %>% rbind(matrix(0, k-k., q))

  A0 <- runif(n*k, -1, 1) %>% matrix(n, k)

  OMEGAp <- map(1:p, ~{rho^abs(.x - (1:p))}) %>% exec(rbind, !!!.)
  OMEGAk <- map(1:k, ~{rho^abs(.x - (1:k))}) %>% exec(rbind, !!!.)

  Rp <- diag(sqrt(1:p))
  Rk <- diag(sqrt(1:k))

  SIGMA <- Rp %*% OMEGAp %*% Rp
  PSI <- Rk %*% OMEGAk %*% Rk

  PSIsq <- expm::sqrtm(PSI)

  A <- A0 %*% PSIsq
  X <- map(time.points, ~{.x^(rev(0:(q-1)))}) %>% exec(rbind, !!!.)

  if(is.null(SNR))
  {
    SNR <- mprod3(PSIsq, THETA, t(X)) %>% crossprod %>% diag %>%
      divide_by(3*diag(SIGMA))
    snr <- mean(SNR)
  } else
  {
    .snr <- mprod3(PSIsq, THETA, t(X)) %>% crossprod %>% diag %>%
      divide_by(3*diag(SIGMA)) %>% mean

    THETA <- THETA*sqrt(SNR/.snr)

    SNR <- mprod3(PSIsq, THETA, t(X)) %>% crossprod %>% diag %>%
      divide_by(3*diag(SIGMA))
    snr <- mean(SNR)
  }

  MU <- A %*% THETA %*% t(X)
  Y <- MU + mvnfast::rmvn(n, numeric(p), SIGMA)

  tp <- seq(min(time.points), max(time.points), length=1000)

  trend <- map(1:k., function(l){
    tp %>% map_dbl(~sum(
      THETA[l,]*(.x^((q-1):0))
    ))
  }) %>% exec(cbind, !!!.) %>% set_colnames(paste0("V", 1:k.)) %>%
    cbind(time.points=tp) %>% as.data.frame %>% pivot_longer(
      all_of(1:k.), names_to = "variable"
    ) %>% ggplot() +
    geom_line(aes(x=time.points, y=.data$value, col=.data$variable), linewidth=1.5) +
    xlab("time points") +
    scale_color_discrete(
      labels = c(
        expression(beta[1]), expression(beta[2]), expression(beta[3]),
        expression(beta[4]), expression(beta[5])
      )
    ) + theme(
      axis.title.y = element_blank(), legend.title = element_blank(),
      text = element_text(size=20)
    )

  return(list(
    Y = Y,
    A = A,
    THETA = THETA,
    X = X,
    MU = MU,
    time.points = time.points,
    SIGMA = SIGMA,
    degrees = 1:k.,
    trend = trend
  ))
}
