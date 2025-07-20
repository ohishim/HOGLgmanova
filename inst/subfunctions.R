# by GMANOVA_HOGL2/GMANOVA.HOGL v1.1.0

#---   lambda candidates   -----------------------------------------------------

lamCand <- function(A.UV, penW, q, lamN, lamR=NULL){

  W <- penW %>% exec(rbind, !!!.)
  lam.max <- max(abs(A.UV) / W)

  if(is.null(lamR))
  {
    lam1 <- sqrt(.Machine$double.eps)
    lamR <- lam1 / lam.max
  }

  return(
    lam.max * exp(
      seq(0, log(lamR), by = log(lamR)/(lamN-1))
    ) %>% rev
  )
}

#---   MM   --------------------------------------------------------------------

GMANOVA.HOGL.fit.MM <- function(
    U0, A, V, A.A, V.V, A.UV, n, k, p, q, lam, penW, L,
    THETA0, alpha, tol
){

  lamW <- penW %>% map(~{lam*.x/L})

  run <- TRUE
  it <- 0
  THETA.aft <- THETA0

  while(run)
  {
    it <- it + 1
    THETA.bef <- THETA.aft

    R <- A.A%*%THETA.bef%*%V.V - A.UV

    THETA.aft <- map(1:k, ~cw.update(
      THETA.bef[.x,] - (R[.x,]/L),
      lamW[[.x]]
    )) %>% exec(rbind, !!!.)

    if(all(THETA.bef == 0))
    {
      dif <- ifelse(all(THETA.aft == 0), 0, 1)
    } else
    {
      dif <- max(((THETA.aft - THETA.bef)^2) / THETA.bef^2, na.rm=T)
    }

    run <- dif > tol
  } #end while run

  rss <- sum((U0 - A%*%tcrossprod(THETA.aft, V))^2)/2
  df <- sum(THETA.aft != 0) + p

  return(list(
    THETA = THETA.aft,
    rss = rss,
    df = df,
    EGCV = rss / ((1 - df/(n*p))^alpha),
    PRSS = rss + sum(map_dbl(1:k, function(l){
      lamw <- lamW[[l]]
      Theta <- THETA.aft[l,]
      out <- map_dbl(1:q, function(j){
        norm2(Theta[1:j])
      }) %>% multiply_by(lamw) %>% sum
      return(out)
    })),
    iter = it
  ))
}

#---   BCDM  -------------------------------------------------------------------

GMANOVA.HOGL.fit.BCDM <- function(
    U0, A, V, A.A, V.V, A.UV, n, k, p, q, lam, penW, L=NULL,
    THETA0, alpha, tol
){

  lamW <- penW %>% map(~{lam*.x})

  run <- TRUE
  it <- 0
  THETA.aft <- THETA0
  A.UtV <- A.UV - A.A%*%THETA0%*%V.V

  while(run)
  {
    it <- it + 1
    THETA.bef <- THETA.aft

    for(l in 1:k)
    {
      .Theta <- THETA.aft[l,]

      THETA.aft[l,] <- Theta <- cw.update(
        A.UtV[l,] + A.A[l,l]*crossprod(.Theta, V.V) %>% drop,
        lamW[[l]]
      )

      A.UtV <- A.UtV - A.A[,l]%*%crossprod(Theta - .Theta, V.V)
    } #end for l

    if(all(THETA.bef == 0))
    {
      dif <- ifelse(all(THETA.aft == 0), 0, 1)
    } else
    {
      dif <- max(((THETA.aft - THETA.bef)^2) / THETA.bef^2, na.rm=T)
    }

    run <- dif > tol
  } #end while run

  rss <- sum((U0 - A%*%tcrossprod(THETA.aft, V))^2)/2
  df <- sum(THETA.aft != 0) + p

  return(list(
    THETA = THETA.aft,
    rss = rss,
    df = df,
    EGCV = rss / ((1 - df/(n*p))^alpha),
    PRSS = rss + sum(map_dbl(1:k, function(l){
      lamw <- lamW[[l]]
      Theta <- THETA.aft[l,]
      out <- map_dbl(1:q, function(j){
        norm2(Theta[1:j])
      }) %>% multiply_by(lamw) %>% sum
      return(out)
    })),
    iter = it
  ))
}

#---   block wise update   -----------------------------------------------------

cw.update <- function(b, Lam){

  q <- length(b)

  get.d <- function(alpha, j, d0=NULL){

    if(j < alpha)
    {
      d <- 0
    } else
    {
      .d <- ifelse(
        j == alpha,
        abs(b[alpha]) - Lam[alpha],
        sqrt((d0^2) + (b[j]^2)) - Lam[j]
      )

      d <- ifelse(.d > 0, .d, 0)
    }

    return(d)
  }

  ##############################################################################
  ###   alpha search
  ##############################################################################

  alpha <- 1
  d.alp <- numeric(q)

  for(j in 1:q)
  {
    l <- j - alpha + 1

    d.alp[l] <- d <- get.d(alpha, j, d.alp[l-1])

    if(d == 0)
    {
      d.alp <- d.alp[-(1:l)]
      alpha <- j + 1
    }
  } #end for j

  ##############################################################################
  ###   get solution
  ##############################################################################

  if(alpha == q + 1)
  {
    Gamma <- numeric(q)
  } else
  {
    idx <- alpha:q

    s <- rev(d.alp / (d.alp + Lam[idx])) %>% cumprod %>% rev %>%
      multiply_by(b[idx])

    if(alpha == 1)
    {
      Gamma <- s
    } else
    {
      Gamma <- c(numeric(alpha-1), s)
    }
  }

  return(Gamma)
}

