# define variables
#   by GMANOVA_HOGL2/functions v1.2.0

defVari <- function(
    Y, A, time.points, n, p, k, q, alpha, ONB
){

  ##############################################################################
  ###   data
  ##############################################################################

  #---   A   -------------------------------------------------------------------

  .A <- A
  A1 <- scale(.A, scale=F)
  scaleA1 <- apply(A1, 2, norm2)
  A <- scale(A1, center=F, scale=scaleA1) %>%
    `attr<-`("scaled:center", NULL) %>% `attr<-`("scaled:scale", NULL)

  labA <- colnames(A)
  if(is.null(labA)){labA <- paste0("A", 1:k)}

  #---   X   -------------------------------------------------------------------

  qidx <- (q-1):0
  .X <- map(time.points, ~{.x^qidx}) %>% exec(rbind, !!!.)
  scaleX <- apply(.X, 2, norm2)
  X <- scale(.X, center=F, scale=scaleX) %>%
    `attr<-`("scaled:center", NULL) %>% `attr<-`("scaled:scale", NULL)

  #---   S   -------------------------------------------------------------------

  A.A <- crossprod(A)
  A.Ainv <- solve(A.A)

  S <- t(Y)%*%(
    diag(n) - matrix(1/n, n, n) - A%*%tcrossprod(A.Ainv, A)
  )%*%Y / (n - k - 1)

  #---   transformation   ------------------------------------------------------

  Sinv <- solve(S)
  Sinvsq <- expm::sqrtm(Sinv)

  U <- Y%*%Sinvsq; V <- Sinvsq%*%X
  U0 <- scale(U, scale=F) %>% `attr<-`("scaled:center", NULL)

  ##############################################################################
  ###   orthonormal basis
  ##############################################################################

  if(ONB)
  {
    QR <- pracma::gramSchmidt(V[,q:1])
    H <- QR$Q[,q:1]
    Q <- QR$R[q:1, q:1]

    V.V <- diag(q)

    #---   LSE   ---------------------------------------------------------------

    A.UV <- mprod3(t(A), U, H)
    LSE <- A.Ainv%*%A.UV
  } else
  {
    #---   LSE   ---------------------------------------------------------------

    V.V <- crossprod(V)
    A.UV <- mprod3(t(A), U, V)
    LSE <- A.Ainv%*%A.UV%*%solve(V.V)

    #---   L   -----------------------------------------------------------------

    L <- max(eigen(A.A)$values) * max(eigen(V.V)$values)
  }

  ##############################################################################
  ###   penalty weight
  ##############################################################################

  penW0 <- map(1:k, function(l){
    Theta <- LSE[l,]
    out <- map_dbl(1:q, function(j){
      1/norm2(Theta[1:j])
    })
    return(out)
  })

  ##############################################################################
  ###   other
  ##############################################################################

  MU <- colMeans(Y) %>% matrix(n, p, byrow=T)

  #---   alpha   ---------------------------------------------------------------

  if(is.null(alpha))
  {
    alpha <- log(n*p)
  } else if(!is.numeric(alpha))
  {
    alpL <- data.frame(
      type = c(2, "log-n", "log-np", "sqrt-n", "sqrt-np"),
      value = c(2, log(n), log(n*p), sqrt(n), sqrt(n*p))
    )

    Alp <- data.frame(type = alpha) %>% left_join(
      alpL, by="type"
    )

    alpha <- Alp$value %>% set_names(Alp$type)
  }

  ##############################################################################
  ###   output
  ##############################################################################

  out <- list(
    A = A, X. = t(X), U0 = U0,
    A.UV = A.UV, A.A = A.A, V.V = V.V,
    penW0 = penW0, MU = MU, LSE = LSE, labA = labA, alpha = alpha,
    scaleA = scaleA1, scaleX = scaleX, .A = .A, .X. = t(.X)
  )

  if(ONB)
  {
    out$V <- H
    out$Q <- Q
    out$L <- NULL
  } else
  {
    out$V <- V
    out$L <- L
  }

  return(out)
}
