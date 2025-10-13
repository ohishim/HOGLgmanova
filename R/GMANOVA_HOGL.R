#' @title Variable and basis function selections in GMANOVA model
#' @description \code{GMANOVA.HOGL} With hierarchical overlapping group Lasso,
#'    this package selects explanatory variables and basis functions in GMANOVA model
#'    (v0.1.0)
#'
#' @importFrom expm sqrtm
#' @importFrom dplyr left_join
#' @importFrom magrittr %>% multiply_by set_names set_colnames set_rownames
#' @importFrom pracma gramSchmidt
#' @importFrom purrr exec map map_dbl
#'
#' @param Y a matrix of response variables
#' @param A a matrix of explanatory variables without intercept
#' @param time.points a vector of time points
#' @param q a dimension of polynomial basis which means a polynomial of degree `q-1`
#' @param alpha a penalty of EGCV criterion; the default (`NULL`) is `log(n*p)`;
#'    this can have one or some positive real numbers or one or some of
#'    `"log-n"`, `"log-np"`, `"sqrt-n"`, and `"sqrt-np"`
#' @param delN the number of candidates of tuning parameter `delta`
#' @param lamN the number of candidates of tuning parameter `lambda`
#' @param lamR the ratio of the minimum `lambda` candidate to the maximum
#' @param ONB if `TRUE` (default), basis matrix is transformed to orthonormal
#'    basis matrix
#' @param output.all if `FALSE`, only optimized estimation result is output
#' @param tol a tolerance for algorithm convergence
#'
#' @return a list with the following elements:
#' \item{intercepts}{a vector of intercepts}
#'
#' \item{coefficients}{a matrix of coefficients}
#'
#' \item{fitted.values}{a matrix of fitted values}
#'
#' \item{polynomial.degrees}{estimated degrees of polynomials}
#'
#' \item{variable.selection}{variable selection results}
#'
#' \item{summary}{summary of estimation results}
#'
#' \item{trend}{figure(s) of varying coefficients}
#'
#' \item{log}{results for all candidates when `output.all=TRUE`}
#'
#' @export
#' @examples
#' #GMANOVA.HOGL(Y, A, time.points, q)

GMANOVA.HOGL <- function(
    Y, A, time.points, q, alpha=NULL, delN=10, lamN=100, lamR=NULL,
    ONB=TRUE, output.all=FALSE, tol=1e-5
){

  t1 <- proc.time()[3] %>% set_names(NULL)

  n <- nrow(Y); p <- length(time.points); k <- ncol(A)

  if(n <= p){stop('`n` must be `n > p`')}
  if(n <= k){stop('`n` must be `n > k`')}
  if(p < q){stop('`p` must be `p >= q`')}

  Vari <- defVari(Y, A, time.points, n, p, k, q, alpha, ONB, norm2)
  for(i in 1:length(Vari)){assign(names(Vari)[i], Vari[[i]])}

  if(ONB)
  {
    GMANOVA.HOGL.fit <- GMANOVA.HOGL.fit.BCDM
  } else
  {
    GMANOVA.HOGL.fit <- GMANOVA.HOGL.fit.MM
  }

  Delta <- seq(0, 1, length=delN+1)[-1]

  ##############################################################################
  ###   delta
  ##############################################################################

  alN <- length(alpha)

  if(alN == 1)
  {
    EGCVlab <- "EGCV"
  } else
  {
    allab <- names(alpha)

    if(is.null(allab))
    {
      names(alpha) <- allab <- 1:alN
    }

    EGCVlab <- paste0("EGCV_", allab)
  }

  LOG <- matrix(0, delN*lamN, 6+alN) %>% set_colnames(c(
    "delta", "lambda", "rss", "df", EGCVlab, "PRSS", "iter"
  ))

  THETA.h <- map(1:alN, ~{LSE})
  egcv.m <- rep(Inf, alN); optD <- optL <- opt <- numeric(alN)
  it <- 0

  for(deli in 1:delN)
  {
    del <- Delta[deli]

    coefW <- c(rep(del, q-1), 1)
    penW <- penW0 %>% map(~{.x*coefW})

    Lambda <- lamCand(A.UV, penW, q, lamN, lamR)

    ############################################################################
    ###   lambda
    ############################################################################

    THETA <- LSE

    for(lami in 1:lamN)
    {
      it <- it + 1
      lam <- Lambda[lami]
      THETA0 <- THETA

      fit <- GMANOVA.HOGL.fit(
        U0, A, V, A.A, V.V, A.UV, n, k, p, q, lam, penW, L, THETA0, alpha,
        tol, cw.update, norm2
      )

      THETA <- fit$THETA; df <- fit$df; egcv <- fit$EGCV

      LOG[it,] <- c(
        del, lam, fit$rss, df, egcv, fit$PRSS, fit$iter
      )

      for(ali in 1:alN)
      {
        if(egcv[ali] < egcv.m[ali])
        {
          THETA.h[[ali]] <- THETA; egcv.m[ali] <- egcv[ali]
          optD[ali] <- deli; optL[ali] <- lami; opt[ali] <- it
        }
      }

      if(df == p){break}

      # print(c(deli, lami))
    } #end for lam
  } #end for del

  t2 <- proc.time()[3] %>% set_names(NULL)

  ##############################################################################
  ###   output
  ##############################################################################

  LOG <- LOG[1:it,] %>% as.data.frame

  if(ONB)
  {
    QQ <- solve(t(Q))
    THETA.h <- THETA.h %>% map(~{.x%*%QQ})
  }

  if(alN == 1){THETA.h <- THETA.h[[1]]}

  if(alN == 1)
  {
    XI <- scale(THETA.h/scaleA, center=F, scale=scaleX) %>%
      `attr<-`("scaled:scale", NULL)
    AXIX. <- mprod3(.A, XI, .X.)
    Eta <- MU[1,] - colMeans(AXIX.)

    res <- list(
      intercepts = Eta,
      coefficients = XI %>% set_rownames(labA) %>%
        set_colnames(paste0("d", (q-1):0)),
      fitted.values = matrix(Eta, n, p, byrow=T) + AXIX.
    )
  } else
  {
    XI <- THETA.h %>% map(~{
      scale(.x/scaleA, center=F, scale=scaleX) %>%
        `attr<-`("scaled:scale", NULL) %>% set_rownames(labA) %>%
        set_colnames(paste0("d", (q-1):0))
    }) %>% set_names(allab)
    Eta <- XI %>% map(~{
      MU[1,] - colMeans(mprod3(.A, .x, .X.))
    })

    res <- list(
      intercepts = Eta,
      coefficients = XI,
      fitted.values = map(1:alN, ~{
        matrix(Eta[[.x]], n, p, byrow=T) + mprod3(.A, XI[[.x]], .X.)
      }) %>% set_names(allab)
    )
  }

  tp <- seq(min(time.points), max(time.points), length=1000)

  if(alN == 1)
  {
    res$polynomial.degrees <- apply(THETA.h, 1, function(x){
      zeroN <- sum(x == 0)
      return(
        ifelse(zeroN == q, 0, q-1-zeroN)
      )
    }) %>% set_names(labA)

    res$variable.selection <- list(
      active = apply(THETA.h, 1, function(x){any(x != 0)}) %>% which,
      inactive = apply(THETA.h, 1, function(x){all(x == 0)}) %>% which
    )

    res$summary <- data.frame(
      delta.idx = optD, lambda.idx = optL,
      delta = LOG[opt, "delta"], lambda = LOG[opt, "lambda"],
      rss = LOG[opt, "rss"], df = LOG[opt, "df"],
      EGCV = LOG[opt, "EGCV"], PRSS = LOG[opt, "PRSS"], time = t2 - t1
    )

    act <- res$variable.selection$active

    if(is.null(act))
    {
      res$trend <- NULL
    } else
    {
      res$trend <- map(act, function(l){
        tp %>% map_dbl(~sum(
          XI[l,]*(.x^((q-1):0))
        ))
      }) %>% exec(cbind, !!!.) %>% set_colnames(paste0("V", act)) %>%
        cbind(time.points=tp) %>% as.data.frame %>% pivot_longer(
          all_of(1:length(act)), names_to = "variable"
        ) %>% ggplot() +
        geom_line(aes(x=time.points, y=.data$value, col=.data$variable), linewidth=1.5) +
        xlab("time points") +
        scale_color_discrete(
          labels = map(act, ~eval(parse(text=paste0("expression(beta[", .x, "])"))))
        ) + theme(
          axis.title.y = element_blank(), legend.title = element_blank(),
          text = element_text(size=20)
        )
    }

  } else
  {
    res$polynomial.degrees <- THETA.h %>% map(~{
      apply(.x, 1, function(x){
        zeroN <- sum(x == 0)
        return(
          ifelse(zeroN == q, 0, q-1-zeroN)
        )
      }) %>% set_names(labA)
    }) %>% exec(rbind, !!!.) %>% set_rownames(allab)

    res$variable.selection <- THETA.h %>% map(~{
      list(
        active = apply(.x, 1, function(x){any(x != 0)}) %>% which,
        inactive = apply(.x, 1, function(x){all(x == 0)}) %>% which
      )
    }) %>% set_names(allab)

    res$summary <- data.frame(
      optD, optL, LOG[opt, -ncol(LOG)], t2 - t1
    ) %>% set_names(c(
      "delta.idx", "lambda.idx", colnames(LOG)[-ncol(LOG)], "time"
    )) %>% set_rownames(allab)

    res$trend <- map(1:alN, function(ali){
      act <- res$variable.selection[[ali]]$active

      if(is.null(act))
      {
        return("none")
      } else
      {
        return(
          map(act, function(l){
            tp %>% map_dbl(~sum(
              XI[[ali]][l,]*(.x^((q-1):0))
            ))
          }) %>% exec(cbind, !!!.) %>% set_colnames(paste0("V", act)) %>%
            cbind(time.points=tp) %>% as.data.frame %>% pivot_longer(
              all_of(1:length(act)), names_to = "variable"
            ) %>% ggplot() +
            geom_line(aes(x=time.points, y=.data$value, col=.data$variable), linewidth=1.5) +
            xlab("time points") +
            scale_color_discrete(
              labels = map(act, ~eval(parse(text=paste0("expression(beta[", .x, "])"))))
            ) + theme(
              axis.title.y = element_blank(), legend.title = element_blank(),
              text = element_text(size=20)
            )
        )
      }
    }) %>% set_names(allab)
  }

  if(output.all)
  {
    res$log = LOG
  }

  return(res)

}
