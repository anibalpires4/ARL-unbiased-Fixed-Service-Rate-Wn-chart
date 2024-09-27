##################################
### M/Ek/1, fixed service rate ###
##################################

library(spc)

# CDF and PDF of V-U, and QF of V-U for the M/Ek/1 system
# k <- 2
cdf.VU1_m <- Vectorize(function(x, lambda, mu, k) {
  if ( x<=0 ) result <- exp(lambda*x)*(k*mu/(k*mu+lambda))^k
  if ( x>0 ) result <- pgamma(x, k, k*mu) + exp(lambda*x)*(k*mu/(k*mu+lambda))^k * ( 1 - pgamma(x, k, k*mu+lambda) )
  result
})  
pdf.VU1_m <- function(x, lambda, mu, k) as.numeric(x<=0)*lambda*cdf.VU1_m(x,lambda,mu,k) + as.numeric(x>0)*( dgamma(x, k, k*mu) + lambda*exp(lambda*x)*(k*mu/(k*mu+lambda))^k * ( 1 - pgamma(x, k, k*mu+lambda) ) - exp(lambda*x)*(k*mu/(k*mu+lambda))^k*dgamma(x, k, k*mu+lambda) )
qf.VU1_m <- function(y, lambda, mu, k) {
  ff <- function(x) y - cdf.VU1_m(x, lambda, mu, k)
  x <- uniroot(ff, c(-50, 50))$root
  x
}


# classical Brook/Evans with one half size and (r-1) full size intervals
BE.arl_m <- function(h, gamma, lambda, mu, k, r=250) {
  ii <- 1:(r-1)
  w <- 2*h/(2*r-1)  
  qij <- function(i,j) cdf.VU1_m( (j-i)*w+w/2, lambda, mu, k) - cdf.VU1_m( (j-i)*w-w/2, lambda, mu, k)
  Qi <- function(i)    cdf.VU1_m(   -i*w+w/2, lambda, mu, k)
  Q <- rbind( cbind( Qi(0), t(qij(0,ii)) ), 
              cbind( Qi(ii), outer(ii,ii,qij) ) )
  Q[,1] <- (1-gamma) * Q[,1]
  one <- array(1, r)
  I <- diag(1, r)
  ARL <- solve(I-Q, one)  
  arl <- 1 + (1-gamma)*Qi(0)*ARL[1] + sum( qij(0,ii) * ARL[-1] )
  arl
}


# Brook/Evans with atom and r full size intervals
BE2.arl_m <- function(h, gamma, lambda, mu, k, r=250) {
  ii <- 1:r
  w <- h/r  
  qij <- function(i,j) cdf.VU1_m( (j-i)*w+w/2, lambda, mu, k) - cdf.VU1_m( (j-i)*w-w/2, lambda, mu, k)
  Qi <- function(i)    cdf.VU1_m(-(i-1)*w-w/2, lambda, mu, k)
  Qj <- function(j)    cdf.VU1_m( j*w, lambda, mu) - cdf.VU1_m( (j-1)*w, lambda, mu, k)
  Q <- rbind( cbind( cdf.VU1_m(0, lambda, mu, k), t(Qj(ii)) ), 
              cbind( Qi(ii), outer(ii,ii,qij) ) )      
  Q[,1] <- (1-gamma) * Q[,1]
  one <- array(1, r+1)
  I <- diag(1, r+1)
  ARL <- solve(I-Q, one)  
  #arl <- 1 + (1-gamma)*cdf.VU1(0, lambda, mu)*ARL[1] + sum( Qj(ii) * ARL[-1] )
  arl <- ARL[1]
  arl
}


# Chebyshev polynomials
Tn <- function(z, n) {
  if ( n==0 ) result <- 1
  if ( n==1 ) result <- z
  if ( n==2 ) result <- 2*z^2 - 1
  if ( n==3 ) result <- 4*z^3 - 3*z
  if ( n==4 ) result <- 8*z^4 - 8*z^2 + 1
  if ( n==5 ) result <- 16*z^5 - 20*z^3 + 5*z
  if ( n>5 )  result <- cos( n*acos(z) ) 
  result
}


# derivatives of Chebyshev polynomials
dTn <- function(z, n) {
  if ( n==0 ) result <- 0
  if ( n==1 ) result <- 1
  if ( n==2 ) result <- 4*z
  if ( n==3 ) result <- 12*z^2 - 3
  if ( n==4 ) result <- 32*z^3 - 16*z
  if ( n==5 ) result <- 80*z^4 - 60*z^2 + 5
  if ( n>5  ) result <- n * ( Tn(z,n-1) - z*Tn(z,n) ) / (1-z^2)
  result
}


# solve ARL integral equation with collocation
coll.arl_m <- function(h, gamma, lambda, mu, k, z0=0, r=40, qm=30) {
  GQ <- quadrature.nodes.weights(qm, type="GL", x1=-1, x2=1)
  z <- GQ$nodes
  w <- GQ$weights
  zch <- h * ( 1 + cos( pi*(2*(r:1)-1)/2/r ) )/2
  A <- matrix(NA, nrow=r, ncol=r)
  for ( i in 1:r ) {
    zi <- zch[i]
    xm1 <- zi/2
    h_ <- qf.VU1_m(1-1e-9, lambda, mu, k) + zi
    if ( h_ > h ) h_ <- h
    xm2 <- (zi+h_)/2
    xw1 <- zi/2
    xw2 <- (h_-zi)/2
    w1 <- xw1 * w
    w2 <- xw2 * w
    z1 <- xm1 + xw1 * z
    z2 <- xm2 + xw2 * z
    for ( j in 1:r ) {
      integral <- sum( w1 * pdf.VU1_m( z1-zi, lambda, mu, k) * Tn(2*z1/h-1, j-1)  ) + sum( w2 * pdf.VU1_m( z2-zi, lambda, mu, k) * Tn(2*z2/h-1, j-1)  )
      A[i,j] <- Tn(2*zi/h-1, j-1) - (1-gamma)*cdf.VU1_m(-zi, lambda, mu, k)*Tn(-1, j-1) - integral
    }
  }
  one <- rep(1, r) 
  g <- solve(A, one)
  result <- 0
  for ( i in 1:r ) result <- result + g[i] * Tn(2*z0/h - 1, i-1)
  result
}


# solve ARL integral equation with collocation
coll2.arl_m <- function(h, gamma, lambda, mu, k, z0=0, r=40, qm=30) {
  GQ <- quadrature.nodes.weights(qm, type="GL", x1=-1, x2=1)
  z <- GQ$nodes
  w <- GQ$weights
  zch <- h * ( 1 + cos( pi*(2*(r:1)-1)/2/r ) )/2
  A <- matrix(NA, nrow=r, ncol=r)
  for ( i in 1:r ) {
    zi <- zch[i]
    xm1 <- zi/2
    xm2 <- (zi+h)/2
    xw1 <- zi/2
    xw2 <- (h-zi)/2
    w1 <- xw1 * w
    w2 <- xw2 * w
    z1 <- xm1 + xw1 * z
    z2 <- xm2 + xw2 * z
    for ( j in 1:r ) {
      integral <- cdf.VU1_m(h-zi, lambda, mu, k) - cdf.VU1_m(-zi, lambda, mu, k)*(-1)^(j-1) - sum( w1 * cdf.VU1_m( z1-zi, lambda, mu, k) * dTn(2*z1/h-1, j-1)*2/h  ) - sum( w2 * cdf.VU1_m( z2-zi, lambda, mu, k) * dTn(2*z2/h-1, j-1)*2/h  )
      A[i,j] <- Tn(2*zi/h-1, j-1) - (1-gamma)*cdf.VU1_m(-zi, lambda, mu, k)*(-1)^(j-1) - integral
    }
  }
  one <- rep(1, r) 
  g <- solve(A, one)
  result <- 0
  for ( i in 1:r ) result <- result + g[i] * Tn(2*z0/h - 1, i-1)
  result
}


# calibrate (determine threshold h) without randomization (aka gamma=0), collocation based
get.h_m <- function(L0, lambda, mu, k, z0=0, r=40, qm=30, OUTPUT=FALSE) {
  h1 <- 1
  L1 <- coll.arl_m(h1, 0, lambda, mu, k, z0=z0, r=r, qm=qm)
  if ( L1 > L0 ) {
    while ( L1 > L0 ) {
      h2 <- h1
      L2 <- L1
      h1 <- h1 * 0.8
      L1 <- coll.arl_m(h1, 0, lambda, mu, k, z0=z0, r=r, qm=qm)
      if ( OUTPUT ) cat(paste("1:\t", h1, "\t", L1, "\n"))
    }
  } else {
    while ( L1 <= L0 ) {
      h2 <- h1
      L2 <- L1
      h1 <- h1 * 1.2
      L1 <- coll.arl_m(h1, 0, lambda, mu, k, z0=z0, r=r, qm=qm)
      if ( OUTPUT ) cat(paste("1:\t", h1, "\t", L1, "\n"))
    }
  }
  h.error <- 1
  L.error <- 1
  while ( abs(h.error) > 1e-9 & abs(L.error) > 1e-9 ) {
    h3 <- h1 + (L0 - L1) / (L2 - L1) * (h2 - h1)
    L3 <- coll.arl_m(h3, 0, lambda, mu, k, z0=z0, r=r, qm=qm)
    if ( OUTPUT ) cat(paste("2:\t", h3, "\t", L3, "\n"))
    h1 <- h2
    h2 <- h3
    L1 <- L2
    L2 <- L3
    h.error <- h2 - h1
    L.error <- L2 - L1
  }
  h3
}


# calibrate (determine threshold h) without randomization (aka gamma=0), MC based
be.get.h_m <- function(L0, lambda, mu, k, r=100, OUTPUT=FALSE) {
  h1 <- 1
  L1 <- BE.arl_m(h1, 0, lambda, mu, k, r=r)
  if ( L1 > L0 ) {
    while ( L1 > L0 ) {
      h2 <- h1
      L2 <- L1
      h1 <- h1 * 0.8
      L1 <- BE.arl_m(h1, 0, lambda, mu, k, r=r)
      if ( OUTPUT ) cat(paste("1:\t", h1, "\t", L1, "\n"))
    }
  } else {
    while ( L1 <= L0 ) {
      h2 <- h1
      L2 <- L1
      h1 <- h1 * 1.2
      L1 <- BE.arl_m(h1, 0, lambda, mu, k, r=r)
      if ( OUTPUT ) cat(paste("1:\t", h1, "\t", L1, "\n"))
    }
  }
  h.error <- 1
  L.error <- 1
  while ( abs(h.error) > 1e-9 & abs(L.error) > 1e-9 ) {
    h3 <- h1 + (L0 - L1) / (L2 - L1) * (h2 - h1)
    L3 <- BE.arl_m(h3, 0, lambda, mu, k, r=r)
    if ( OUTPUT ) cat(paste("2:\t", h3, "\t", L3, "\n"))
    h1 <- h2
    h2 <- h3
    L1 <- L2
    L2 <- L3
    h.error <- h2 - h1
    L.error <- L2 - L1
  }
  h3
}


# search for gamma so that for fixed threshold the ARL becomes pre-defined L0, collocation based
get.gamma_m <- function(h, L0, lambda, mu, k, z0=0, r=40, qm=30, OUTPUT=FALSE) {
  L0.full_m <- coll.arl_m(h, 0, lambda, mu, k, z0=z0, r=r, qm=qm)
  if ( OUTPUT ) cat(paste("0:\t", 0, "\t", L0.full_m, "\n"))
  if ( L0.full_m <= L0 ) stop("h too small or L0 too large")
  g1 <- .05
  L1 <- coll.arl_m(h, g1, lambda, mu, k, z0=z0, r=r, qm=qm)
  if ( L1 > L0 ) {
    while ( L1 > L0 ) {
      g2 <- g1
      L2 <- L1
      g1 <- g1 * 1.2
      L1 <- coll.arl_m(h, g1, lambda, mu, k, z0=z0, r=r, qm=qm)
      if ( OUTPUT ) cat(paste("1:\t", g1, "\t", L1, "\n"))
    }
  } else {
    while ( L1 <= L0 ) {
      g2 <- g1
      L2 <- L1
      g1 <- g1 * 0.8
      L1 <- coll.arl_m(h, g1, lambda, mu, k, z0=z0, r=r, qm=qm)
      if ( OUTPUT ) cat(paste("1:\t", g1, "\t", L1, "\n"))
    }
  }
  g.error <- 1
  L.error <- 1
  while ( abs(g.error) > 1e-9 & abs(L.error) > 1e-9 ) {
    g3 <- g1 + (L0 - L1) / (L2 - L1) * (g2 - g1)
    L3 <- coll.arl_m(h, g3, lambda, mu, k, z0=z0, r=r, qm=qm)
    if ( OUTPUT ) cat(paste("2:\t", g3, "\t", L3, "\n"))
    g1 <- g2
    g2 <- g3
    L1 <- L2
    L2 <- L3
    g.error <- g2 - g1
    L.error <- L2 - L1
  }
  g3
}


# search for gamma so that for fixed threshold the ARL becomes pre-defined L0, MC based
be.get.gamma_m <- function(h, L0, lambda, mu, k, r=100, OUTPUT=FALSE) {
  L0.full <- BE2.arl_m(h, 0, lambda, mu, k, r=r)
  if ( OUTPUT ) cat(paste("0:\t", 0, "\t", L0.full_m, "\n"))
  if ( L0.full_m <= L0 ) stop("h too small or L0 too large")
  g1 <- .05
  L1 <- BE2.arl_m(h, g1, lambda, mu, k, r=r)
  if ( L1 > L0 ) {
    while ( L1 > L0 ) {
      g2 <- g1
      L2 <- L1
      g1 <- g1 * 1.2
      L1 <- BE2.arl_m(h, g1, lambda, mu, k, r=r)
      if ( OUTPUT ) cat(paste("1:\t", g1, "\t", L1, "\n"))
    }
  } else {
    while ( L1 <= L0 ) {
      g2 <- g1
      L2 <- L1
      g1 <- g1 * 0.8
      L1 <- BE2.arl_m(h, g1, lambda, mu, k, r=r)
      if ( OUTPUT ) cat(paste("1:\t", g1, "\t", L1, "\n"))
    }
  }
  g.error <- 1
  L.error <- 1
  while ( abs(g.error) > 1e-9 & abs(L.error) > 1e-9 ) {
    g3 <- g1 + (L0 - L1) / (L2 - L1) * (g2 - g1)
    L3 <- BE2.arl_m(h, g3, lambda, mu, k, r=r)
    if ( OUTPUT ) cat(paste("2:\t", g3, "\t", L3, "\n"))
    g1 <- g2
    g2 <- g3
    L1 <- L2
    L2 <- L3
    g.error <- g2 - g1
    L.error <- L2 - L1
  }
  g3
}


# calculate (h, gamma) so that (i) ARL = L0 and (ii) ARL as function of lambda takes max at given lambda value, collocation based
get.gammaAh.lambda_m <- function(L0, lambda, mu, k, z0=0, r=40, qm=30, eps=1e-6, OUTPUT=FALSE) {
  lm <- lambda - eps
  lp <- lambda + eps
  
  h1 <- get.h_m(2*L0, lambda, mu, k, z0=z0, r=r, qm=qm)
  g1 <- get.gamma_m(h1, L0, lambda, mu, k, z0=z0, r=r, qm=qm)
  ARLm <- coll.arl_m(h1, g1, lm, mu, k, z0=z0, r=r, qm=qm)
  ARLp <- coll.arl_m(h1, g1, lp, mu, k, z0=z0, r=r, qm=qm)
  dratio1 <- ( ARLm - ARLp ) / ( 2*eps )
  if ( OUTPUT ) cat(paste("0\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
  
  if ( dratio1 > 0 ) {
    while ( dratio1 > 0 ) {
      h2 <- h1
      dratio2 <- dratio1
      h1 <- h1 * 1.2
      g1 <- get.gamma_m(h1, L0, lambda, mu, k, z0=z0, r=r, qm=qm)
      ARLm <- coll.arl_m(h1, g1, lm, mu, k, z0=z0, r=r, qm=qm)
      ARLp <- coll.arl_m(h1, g1, lp, mu, k, z0=z0, r=r, qm=qm)
      dratio1 <- ( ARLm - ARLp ) / ( 2*eps )
      if ( OUTPUT ) cat(paste("1\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
    }
  } else {
    while ( dratio1 <= 0 ) {
      h2 <- h1
      dratio2 <- dratio1
      h1 <- h1 * 0.8
      g1 <- get.gamma_m(h1, L0, lambda, mu, k, z0=z0, r=r, qm=qm)
      ARLm <- coll.arl_m(h1, g1, lm, mu, k, z0=z0, r=r, qm=qm)
      ARLp <- coll.arl_m(h1, g1, lp, mu, k, z0=z0, r=r, qm=qm)
      dratio1 <- ( ARLm - ARLp ) / ( 2*eps )
      if ( OUTPUT ) cat(paste("1\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
    }
  }
  
  h.error <- 1
  dr.error <- 1
  while ( abs(h.error) > 1e-10 & abs(dr.error) > 1e-8 ) {
    h3 <- h1 - dratio1 / (dratio2 - dratio1) * (h2 - h1)
    g3 <- get.gamma_m(h3, L0, lambda, mu, k, z0=z0, r=r, qm=qm)
    ARLm <- coll.arl_m(h3, g3, lm, mu, k, z0=z0, r=r, qm=qm)
    ARLp <- coll.arl_m(h3, g3, lp, mu, k, z0=z0, r=r, qm=qm)
    dratio3 <- ( ARLm - ARLp ) / ( 2*eps )
    if ( OUTPUT ) cat(paste("2\th3 =", h3, "\tg3 =", g3, "\tdratio1 =", dratio3, "\n"))
    h1 <- h2
    h2 <- h3
    dratio1 <- dratio2
    dratio2 <- dratio3
    h.error <- h2 - h1
    dr.error <- dratio2 - dratio1
  }
  data.frame(h=h3, gamma=g3)
}


# calculate (h, gamma) so that (i) ARL = L0 and (ii) ARL as function of lambda takes max at given lambda value, MC based
be.get.gammaAh.lambda_m <- function(L0, lambda, mu, k, r=100, eps=1e-6, OUTPUT=FALSE) {
  lm <- lambda - eps
  lp <- lambda + eps
  
  h1 <- be.get.h_m(2*L0, lambda, mu, k, r=r)
  g1 <- be.get.gamma_m(h1, L0, lambda, mu, k, r=r)
  ARLm <- BE2.arl_m(h1, g1, lm, mu, k, r=r)
  ARLp <- BE2.arl_m(h1, g1, lp, mu, k, r=r)
  dratio1 <- ( ARLm - ARLp ) / ( 2*eps )
  if ( OUTPUT ) cat(paste("0\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
  
  if ( dratio1 > 0 ) {
    while ( dratio1 > 0 ) {
      h2 <- h1
      dratio2 <- dratio1
      h1 <- h1 * 1.2
      g1 <- be.get.gamma_m(h1, L0, lambda, mu, k, r=r)
      ARLm <- BE2.arl_m(h1, g1, lm, mu, k, r=r)
      ARLp <- BE2.arl_m(h1, g1, lp, mu, k, r=r)
      dratio1 <- ( ARLm - ARLp ) / ( 2*eps )
      if ( OUTPUT ) cat(paste("1\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
    }
  } else {
    while ( dratio1 <= 0 ) {
      h2 <- h1
      dratio2 <- dratio1
      h1 <- h1 * 0.8
      g1 <- be.get.gamma_m(h1, L0, lambda, mu, k, r=r)
      ARLm <- BE2.arl_m(h1, g1, lm, mu, k, r=r)
      ARLp <- BE2.arl_m(h1, g1, lp, mu, k, r=r)
      dratio1 <- ( ARLm - ARLp ) / ( 2*eps )
      if ( OUTPUT ) cat(paste("1\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
    }
  }
  
  h.error <- 1
  dr.error <- 1
  while ( abs(h.error) > 1e-10 & abs(dr.error) > 1e-8 ) {
    h3 <- h1 - dratio1 / (dratio2 - dratio1) * (h2 - h1)
    g3 <- be.get.gamma_m(h3, L0, lambda, mu, k, r=r)
    ARLm <- BE2.arl_m(h3, g3, lm, mu, k, r=r)
    ARLp <- BE2.arl_m(h3, g3, lp, mu, k, r=r)
    dratio3 <- ( ARLm - ARLp ) / ( 2*eps )
    if ( OUTPUT ) cat(paste("2\th3 =", h3, "\tg3 =", g3, "\tdratio1 =", dratio3, "\n"))
    h1 <- h2
    h2 <- h3
    dratio1 <- dratio2
    dratio2 <- dratio3
    h.error <- h2 - h1
    dr.error <- dratio2 - dratio1
  }
  data.frame(h=h3, gamma=g3)
}


# calculate (h, gamma) so that (i) ARL = L0 and (ii) ARL as function of mu takes max at given mu value, collocation based
get.gammaAh.mu_m <- function(L0, lambda, mu, k, z0=0, r=40, qm=30, eps=1e-6, OUTPUT=FALSE) {
  mm <- mu - eps
  mp <- mu + eps
  
  h1 <- get.h_m(2*L0, lambda, mu, k, z0=z0, r=r, qm=qm)
  g1 <- get.gamma_m(h1, L0, lambda, mu, k, z0=z0, r=r, qm=qm)
  ARLm <- coll.arl_m(h1, g1, lambda, mm, k, z0=z0, r=r, qm=qm)
  ARLp <- coll.arl_m(h1, g1, lambda, mp, k, z0=z0, r=r, qm=qm)
  dratio1 <- ( ARLm - ARLp ) / ( 2*eps )
  if ( OUTPUT ) cat(paste("0\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
  
  if ( dratio1 > 0 ) {
    while ( dratio1 > 0 ) {
      h2 <- h1
      dratio2 <- dratio1
      h1 <- h1 * 0.8
      g1 <- get.gamma_m(h1, L0, lambda, mu, k, z0=z0, r=r, qm=qm)
      ARLm <- coll.arl_m(h1, g1, lambda, mm, k, z0=z0, r=r, qm=qm)
      ARLp <- coll.arl_m(h1, g1, lambda, mp, k, z0=z0, r=r, qm=qm)
      dratio1 <- ( ARLm - ARLp ) / ( 2*eps )
      if ( OUTPUT ) cat(paste("1\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
    }
  } else {
    while ( dratio1 <= 0 ) {
      h2 <- h1
      dratio2 <- dratio1
      h1 <- h1 * 1.2
      g1 <- get.gamma_m(h1, L0, lambda, mu, k, z0=z0, r=r, qm=qm)
      ARLm <- coll.arl_m(h1, g1, lambda, mm, k, z0=z0, r=r, qm=qm)
      ARLp <- coll.arl_m(h1, g1, lambda, mp, k, z0=z0, r=r, qm=qm)
      dratio1 <- ( ARLm - ARLp ) / ( 2*eps )
      if ( OUTPUT ) cat(paste("1\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
    }
  }
  
  h.error <- 1
  dr.error <- 1
  while ( abs(h.error) > 1e-10 & abs(dr.error) > 1e-8 ) {
    h3 <- h1 - dratio1 / (dratio2 - dratio1) * (h2 - h1)
    g3 <- get.gamma_m(h3, L0, lambda, mu, k, z0=z0, r=r, qm=qm)
    ARLm <- coll.arl_m(h3, g3, lambda, mm, k, z0=z0, r=r, qm=qm)
    ARLp <- coll.arl_m(h3, g3, lambda, mp, k, z0=z0, r=r, qm=qm)
    dratio3 <- ( ARLm - ARLp ) / ( 2*eps )
    if ( OUTPUT ) cat(paste("2\th3 =", h3, "\tg3 =", g3, "\tdratio1 =", dratio3, "\n"))
    h1 <- h2
    h2 <- h3
    dratio1 <- dratio2
    dratio2 <- dratio3
    h.error <- h2 - h1
    dr.error <- dratio2 - dratio1
  }
  data.frame(h=h3, gamma=g3)
}


# calculate (h, gamma) so that (i) ARL = L0 and (ii) ARL as function of mu takes max at given mu value, MC based
be.get.gammaAh.mu_m <- function(L0, lambda, mu, k, z0=0, r=100, eps=1e-6, OUTPUT=FALSE) {
  mm <- mu - eps
  mp <- mu + eps
  
  h1 <- be.get.h_m(2*L0, lambda, mu, r=r)
  g1 <- be.get.gamma_m(h1, L0, lambda, mu, k, r=r)
  ARLm <- BE2.arl_m(h1, g1, lambda, mm, k, r=r)
  ARLp <- BE2.arl_m(h1, g1, lambda, mp, k, r=r)
  dratio1 <- ( ARLm - ARLp ) / ( 2*eps )
  if ( OUTPUT ) cat(paste("0\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
  
  if ( dratio1 > 0 ) {
    while ( dratio1 > 0 ) {
      h2 <- h1
      dratio2 <- dratio1
      h1 <- h1 * 0.8
      g1 <- be.get.gamma_m(h1, L0, lambda, mu, k, r=r)
      ARLm <- BE2.arl_m(h1, g1, lambda, mm, k, r=r)
      ARLp <- BE2.arl_m(h1, g1, lambda, mp, k, r=r)
      dratio1 <- ( ARLm - ARLp ) / ( 2*eps )
      if ( OUTPUT ) cat(paste("1\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
    }
  } else {
    while ( dratio1 <= 0 ) {
      h2 <- h1
      dratio2 <- dratio1
      h1 <- h1 * 1.2
      g1 <- be.get.gamma_m(h1, L0, lambda, mu, k, r=r)
      ARLm <- BE2.arl_m(h1, g1, lambda, mm, k, r=r)
      ARLp <- BE2.arl_m(h1, g1, lambda, mp, k, r=r)
      dratio1 <- ( ARLm - ARLp ) / ( 2*eps )
      if ( OUTPUT ) cat(paste("1\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
    }
  }
  
  h.error <- 1
  dr.error <- 1
  while ( abs(h.error) > 1e-10 & abs(dr.error) > 1e-8 ) {
    h3 <- h1 - dratio1 / (dratio2 - dratio1) * (h2 - h1)
    g3 <- get.gamma_m(h3, L0, lambda, mu, k, z0=z0, r=r)
    ARLm <- BE2.arl_m(h3, g3, lambda, mm, k, r=r)
    ARLp <- BE2.arl_m(h3, g3, lambda, mp, k, r=r)
    dratio3 <- ( ARLm - ARLp ) / ( 2*eps )
    if ( OUTPUT ) cat(paste("2\th3 =", h3, "\tg3 =", g3, "\tdratio1 =", dratio3, "\n"))
    h1 <- h2
    h2 <- h3
    dratio1 <- dratio2
    dratio2 <- dratio3
    h.error <- h2 - h1
    dr.error <- dratio2 - dratio1
  }
  data.frame(h=h3, gamma=g3)
}


##################################
### Ek/M/1, fixed service rate ###
##################################


# CDF and PDF of V-U, and QF of V-U for the Ek/M/1 system
#k <- 2
cdf.VU1 <- Vectorize(function(x, lambda, mu, k) {
  if ( x<=0 ) result <- 1-pgamma(-x, k, k*lambda) - exp(-mu*x)*(k*lambda/(k*lambda+mu))^k*(1-pgamma(-x, k, k*lambda+mu)) 
  if ( x>0 ) result <- 1 - exp(-mu*x)*(k*lambda/(k*lambda+mu))^k
  result
})  
pdf.VU1 <- function(x, lambda, mu, k) as.numeric(x<=0)*( dgamma(-x, k, k*lambda) + mu*exp(-mu*x)*(k*lambda/(k*lambda+mu))^k*(1-pgamma(-x, k, k*lambda+mu)) - exp(-mu*x)*(k*lambda/(k*lambda+mu))^k*dgamma(-x, k, k*lambda+mu) ) + as.numeric(x>0)*mu*exp(-mu*x)*(k*lambda/(k*lambda+mu))^k  

qf.VU1 <- function(y, lambda, mu, k) {
  ff <- function(x) y - cdf.VU1(x, lambda, mu, k)
  x <- uniroot(ff, c(-50, 50))$root
  x
}


# classical Brook/Evans with one half size and (r-1) full size intervals
BE.arl <- function(h, gamma, lambda, mu, k, r=250) {
  ii <- 1:(r-1)
  w <- 2*h/(2*r-1)  
  qij <- function(i,j) cdf.VU1( (j-i)*w+w/2, lambda, mu, k) - cdf.VU1( (j-i)*w-w/2, lambda, mu, k)
  Qi <- function(i)    cdf.VU1(   -i*w+w/2, lambda, mu, k)
  Q <- rbind( cbind( Qi(0), t(qij(0,ii)) ), 
              cbind( Qi(ii), outer(ii,ii,qij) ) )
  Q[,1] <- (1-gamma) * Q[,1]
  one <- array(1, r)
  I <- diag(1, r)
  ARL <- solve(I-Q, one)  
  arl <- 1 + (1-gamma)*Qi(0)*ARL[1] + sum( qij(0,ii) * ARL[-1] )
  arl
}


# Brook/Evans with atom and r full size intervals
BE2.arl <- function(h, gamma, lambda, mu, k, r=250) {
  ii <- 1:r
  w <- h/r  
  qij <- function(i,j) cdf.VU1( (j-i)*w+w/2, lambda, mu, k) - cdf.VU1( (j-i)*w-w/2, lambda, mu, k)
  Qi <- function(i)    cdf.VU1(-(i-1)*w-w/2, lambda, mu, k)
  Qj <- function(j)    cdf.VU1( j*w, lambda, mu, k) - cdf.VU1( (j-1)*w, lambda, mu, k)
  Q <- rbind( cbind( cdf.VU1(0, lambda, mu, k), t(Qj(ii)) ), 
              cbind( Qi(ii), outer(ii,ii,qij) ) )      
  Q[,1] <- (1-gamma) * Q[,1]
  one <- array(1, r+1)
  I <- diag(1, r+1)
  ARL <- solve(I-Q, one)  
  #arl <- 1 + (1-gamma)*cdf.VU1(0, lambda, mu)*ARL[1] + sum( Qj(ii) * ARL[-1] )
  arl <- ARL[1]
  arl
}


# Chebyshev polynomials
Tn <- function(z, n) {
  if ( n==0 ) result <- 1
  if ( n==1 ) result <- z
  if ( n==2 ) result <- 2*z^2 - 1
  if ( n==3 ) result <- 4*z^3 - 3*z
  if ( n==4 ) result <- 8*z^4 - 8*z^2 + 1
  if ( n==5 ) result <- 16*z^5 - 20*z^3 + 5*z
  if ( n>5 )  result <- cos( n*acos(z) ) 
  result
}


# derivatives of Chebyshev polynomials
dTn <- function(z, n) {
  if ( n==0 ) result <- 0
  if ( n==1 ) result <- 1
  if ( n==2 ) result <- 4*z
  if ( n==3 ) result <- 12*z^2 - 3
  if ( n==4 ) result <- 32*z^3 - 16*z
  if ( n==5 ) result <- 80*z^4 - 60*z^2 + 5
  if ( n>5  ) result <- n * ( Tn(z,n-1) - z*Tn(z,n) ) / (1-z^2)
  result
}


# solve ARL integral equation with collocation
coll.arl <- function(h, gamma, lambda, mu, k, z0=0, r=40, qm=30) {
  GQ <- quadrature.nodes.weights(qm, type="GL", x1=-1, x2=1)
  z <- GQ$nodes
  w <- GQ$weights
  zch <- h * ( 1 + cos( pi*(2*(r:1)-1)/2/r ) )/2
  A <- matrix(NA, nrow=r, ncol=r)
  for ( i in 1:r ) {
    zi <- zch[i]
    xm1 <- zi/2
    h_ <- qf.VU1(1-1e-9, lambda, mu, k) + zi
    if ( h_ > h ) h_ <- h
    xm2 <- (zi+h_)/2
    xw1 <- zi/2
    xw2 <- (h_-zi)/2
    w1 <- xw1 * w
    w2 <- xw2 * w
    z1 <- xm1 + xw1 * z
    z2 <- xm2 + xw2 * z
    for ( j in 1:r ) {
      integral <- sum( w1 * pdf.VU1( z1-zi, lambda, mu, k) * Tn(2*z1/h-1, j-1)  ) + sum( w2 * pdf.VU1( z2-zi, lambda, mu, k) * Tn(2*z2/h-1, j-1)  )
      A[i,j] <- Tn(2*zi/h-1, j-1) - (1-gamma)*cdf.VU1(-zi, lambda, mu, k)*Tn(-1, j-1) - integral
    }
  }
  one <- rep(1, r) 
  g <- solve(A, one)
  result <- 0
  for ( i in 1:r ) result <- result + g[i] * Tn(2*z0/h - 1, i-1)
  result
}


# solve ARL integral equation with collocation
coll2.arl <- function(h, gamma, lambda, mu, k, z0=0, r=40, qm=30) {
  GQ <- quadrature.nodes.weights(qm, type="GL", x1=-1, x2=1)
  z <- GQ$nodes
  w <- GQ$weights
  zch <- h * ( 1 + cos( pi*(2*(r:1)-1)/2/r ) )/2
  A <- matrix(NA, nrow=r, ncol=r)
  for ( i in 1:r ) {
    zi <- zch[i]
    xm1 <- zi/2
    xm2 <- (zi+h)/2
    xw1 <- zi/2
    xw2 <- (h-zi)/2
    w1 <- xw1 * w
    w2 <- xw2 * w
    z1 <- xm1 + xw1 * z
    z2 <- xm2 + xw2 * z
    for ( j in 1:r ) {
      integral <- cdf.VU1(h-zi, lambda, mu, k) - cdf.VU1(-zi, lambda, mu, k)*(-1)^(j-1) - sum( w1 * cdf.VU1( z1-zi, lambda, mu, k) * dTn(2*z1/h-1, j-1)*2/h  ) - sum( w2 * cdf.VU1( z2-zi, lambda, mu, k) * dTn(2*z2/h-1, j-1)*2/h  )
      A[i,j] <- Tn(2*zi/h-1, j-1) - (1-gamma)*cdf.VU1(-zi, lambda, mu, k)*(-1)^(j-1) - integral
    }
  }
  one <- rep(1, r) 
  g <- solve(A, one)
  result <- 0
  for ( i in 1:r ) result <- result + g[i] * Tn(2*z0/h - 1, i-1)
  result
}


# calibrate (determine threshold h) without randomization (aka gamma=0), collocation based
get.h <- function(L0, lambda, mu, k, z0=0, r=40, qm=30, OUTPUT=FALSE) {
  h1 <- 1
  L1 <- coll.arl(h1, 0, lambda, mu, k, z0=z0, r=r, qm=qm)
  if ( L1 > L0 ) {
    while ( L1 > L0 ) {
      h2 <- h1
      L2 <- L1
      h1 <- h1 * 0.8
      L1 <- coll.arl(h1, 0, lambda, mu, k, z0=z0, r=r, qm=qm)
      if ( OUTPUT ) cat(paste("1:\t", h1, "\t", L1, "\n"))
    }
  } else {
    while ( L1 <= L0 ) {
      h2 <- h1
      L2 <- L1
      h1 <- h1 * 1.2
      L1 <- coll.arl(h1, 0, lambda, mu, k, z0=z0, r=r, qm=qm)
      if ( OUTPUT ) cat(paste("1:\t", h1, "\t", L1, "\n"))
    }
  }
  h.error <- 1
  L.error <- 1
  while ( abs(h.error) > 1e-9 & abs(L.error) > 1e-9 ) {
    h3 <- h1 + (L0 - L1) / (L2 - L1) * (h2 - h1)
    L3 <- coll.arl(h3, 0, lambda, mu, k, z0=z0, r=r, qm=qm)
    if ( OUTPUT ) cat(paste("2:\t", h3, "\t", L3, "\n"))
    h1 <- h2
    h2 <- h3
    L1 <- L2
    L2 <- L3
    h.error <- h2 - h1
    L.error <- L2 - L1
  }
  h3
}


# calibrate (determine threshold h) without randomization (aka gamma=0), MC based
be.get.h <- function(L0, lambda, mu, k, r=100, OUTPUT=FALSE) {
  h1 <- 1
  L1 <- BE.arl(h1, 0, lambda, mu, k, r=r)
  if ( L1 > L0 ) {
    while ( L1 > L0 ) {
      h2 <- h1
      L2 <- L1
      h1 <- h1 * 0.8
      L1 <- BE.arl(h1, 0, lambda, mu, k, r=r)
      if ( OUTPUT ) cat(paste("1:\t", h1, "\t", L1, "\n"))
    }
  } else {
    while ( L1 <= L0 ) {
      h2 <- h1
      L2 <- L1
      h1 <- h1 * 1.2
      L1 <- BE.arl(h1, 0, lambda, mu, k, r=r)
      if ( OUTPUT ) cat(paste("1:\t", h1, "\t", L1, "\n"))
    }
  }
  h.error <- 1
  L.error <- 1
  while ( abs(h.error) > 1e-9 & abs(L.error) > 1e-9 ) {
    h3 <- h1 + (L0 - L1) / (L2 - L1) * (h2 - h1)
    L3 <- BE.arl(h3, 0, lambda, mu, k, r=r)
    if ( OUTPUT ) cat(paste("2:\t", h3, "\t", L3, "\n"))
    h1 <- h2
    h2 <- h3
    L1 <- L2
    L2 <- L3
    h.error <- h2 - h1
    L.error <- L2 - L1
  }
  h3
}


# search for gamma so that for fixed threshold the ARL becomes pre-defined L0, collocation based
get.gamma <- function(h, L0, lambda, mu, k, z0=0, r=40, qm=30, OUTPUT=FALSE) {
  L0.full <- coll.arl(h, 0, lambda, mu, k, z0=z0, r=r, qm=qm)
  if ( OUTPUT ) cat(paste("0:\t", 0, "\t", L0.full, "\n"))
  if ( L0.full <= L0 ) stop("h too small or L0 too large")
  g1 <- .05
  L1 <- coll.arl(h, g1, lambda, mu, k, z0=z0, r=r, qm=qm)
  if ( L1 > L0 ) {
    while ( L1 > L0 ) {
      g2 <- g1
      L2 <- L1
      g1 <- g1 * 1.2
      L1 <- coll.arl(h, g1, lambda, mu, k, z0=z0, r=r, qm=qm)
      if ( OUTPUT ) cat(paste("1:\t", g1, "\t", L1, "\n"))
    }
  } else {
    while ( L1 <= L0 ) {
      g2 <- g1
      L2 <- L1
      g1 <- g1 * 0.8
      L1 <- coll.arl(h, g1, lambda, mu, k, z0=z0, r=r, qm=qm)
      if ( OUTPUT ) cat(paste("1:\t", g1, "\t", L1, "\n"))
    }
  }
  g.error <- 1
  L.error <- 1
  while ( abs(g.error) > 1e-9 & abs(L.error) > 1e-9 ) {
    g3 <- g1 + (L0 - L1) / (L2 - L1) * (g2 - g1)
    L3 <- coll.arl(h, g3, lambda, mu, k, z0=z0, r=r, qm=qm)
    if ( OUTPUT ) cat(paste("2:\t", g3, "\t", L3, "\n"))
    g1 <- g2
    g2 <- g3
    L1 <- L2
    L2 <- L3
    g.error <- g2 - g1
    L.error <- L2 - L1
  }
  g3
}


# search for gamma so that for fixed threshold the ARL becomes pre-defined L0, MC based
be.get.gamma <- function(h, L0, lambda, mu, k, r=100, OUTPUT=FALSE) {
  L0.full <- BE2.arl(h, 0, lambda, mu, k, r=r)
  if ( OUTPUT ) cat(paste("0:\t", 0, "\t", L0.full, "\n"))
  if ( L0.full <= L0 ) stop("h too small or L0 too large")
  g1 <- .05
  L1 <- BE2.arl(h, g1, lambda, mu, k, r=r)
  if ( L1 > L0 ) {
    while ( L1 > L0 ) {
      g2 <- g1
      L2 <- L1
      g1 <- g1 * 1.2
      L1 <- BE2.arl(h, g1, lambda, mu, k, r=r)
      if ( OUTPUT ) cat(paste("1:\t", g1, "\t", L1, "\n"))
    }
  } else {
    while ( L1 <= L0 ) {
      g2 <- g1
      L2 <- L1
      g1 <- g1 * 0.8
      L1 <- BE2.arl(h, g1, lambda, mu, k, r=r)
      if ( OUTPUT ) cat(paste("1:\t", g1, "\t", L1, "\n"))
    }
  }
  g.error <- 1
  L.error <- 1
  while ( abs(g.error) > 1e-9 & abs(L.error) > 1e-9 ) {
    g3 <- g1 + (L0 - L1) / (L2 - L1) * (g2 - g1)
    L3 <- BE2.arl(h, g3, lambda, mu, k, r=r)
    if ( OUTPUT ) cat(paste("2:\t", g3, "\t", L3, "\n"))
    g1 <- g2
    g2 <- g3
    L1 <- L2
    L2 <- L3
    g.error <- g2 - g1
    L.error <- L2 - L1
  }
  g3
}


# calculate (h, gamma) so that (i) ARL = L0 and (ii) ARL as function of lambda takes max at given lambda value, collocation based
get.gammaAh.lambda <- function(L0, lambda, mu, k, z0=0, r=40, qm=30, eps=1e-6, OUTPUT=FALSE) {
  lm <- lambda - eps
  lp <- lambda + eps
  
  h1 <- get.h(2*L0, lambda, mu, k, z0=z0, r=r, qm=qm)
  g1 <- get.gamma(h1, L0, lambda, mu, k, z0=z0, r=r, qm=qm)
  ARLm <- coll.arl(h1, g1, lm, mu, k, z0=z0, r=r, qm=qm)
  ARLp <- coll.arl(h1, g1, lp, mu, k, z0=z0, r=r, qm=qm)
  dratio1 <- ( ARLm - ARLp ) / ( 2*eps )
  if ( OUTPUT ) cat(paste("0\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
  
  if ( dratio1 > 0 ) {
    while ( dratio1 > 0 ) {
      h2 <- h1
      dratio2 <- dratio1
      h1 <- h1 * 1.2
      g1 <- get.gamma(h1, L0, lambda, mu, k, z0=z0, r=r, qm=qm)
      ARLm <- coll.arl(h1, g1, lm, mu, k, z0=z0, r=r, qm=qm)
      ARLp <- coll.arl(h1, g1, lp, mu, k, z0=z0, r=r, qm=qm)
      dratio1 <- ( ARLm - ARLp ) / ( 2*eps )
      if ( OUTPUT ) cat(paste("1\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
    }
  } else {
    while ( dratio1 <= 0 ) {
      h2 <- h1
      dratio2 <- dratio1
      h1 <- h1 * 0.8
      g1 <- get.gamma(h1, L0, lambda, mu, k, z0=z0, r=r, qm=qm)
      ARLm <- coll.arl(h1, g1, lm, mu, k, z0=z0, r=r, qm=qm)
      ARLp <- coll.arl(h1, g1, lp, mu, k, z0=z0, r=r, qm=qm)
      dratio1 <- ( ARLm - ARLp ) / ( 2*eps )
      if ( OUTPUT ) cat(paste("1\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
    }
  }
  
  h.error <- 1
  dr.error <- 1
  while ( abs(h.error) > 1e-10 & abs(dr.error) > 1e-8 ) {
    h3 <- h1 - dratio1 / (dratio2 - dratio1) * (h2 - h1)
    g3 <- get.gamma(h3, L0, lambda, mu, k, z0=z0, r=r, qm=qm)
    ARLm <- coll.arl(h3, g3, lm, mu, k, z0=z0, r=r, qm=qm)
    ARLp <- coll.arl(h3, g3, lp, mu, k, z0=z0, r=r, qm=qm)
    dratio3 <- ( ARLm - ARLp ) / ( 2*eps )
    if ( OUTPUT ) cat(paste("2\th3 =", h3, "\tg3 =", g3, "\tdratio1 =", dratio3, "\n"))
    h1 <- h2
    h2 <- h3
    dratio1 <- dratio2
    dratio2 <- dratio3
    h.error <- h2 - h1
    dr.error <- dratio2 - dratio1
  }
  data.frame(h=h3, gamma=g3)
}


# calculate (h, gamma) so that (i) ARL = L0 and (ii) ARL as function of lambda takes max at given lambda value, MC based
be.get.gammaAh.lambda <- function(L0, lambda, mu, k, r=100, eps=1e-6, OUTPUT=FALSE) {
  lm <- lambda - eps
  lp <- lambda + eps
  
  h1 <- be.get.h(2*L0, lambda, mu, k, r=r)
  g1 <- be.get.gamma(h1, L0, lambda, mu, k, r=r)
  ARLm <- BE2.arl(h1, g1, lm, mu, k, r=r)
  ARLp <- BE2.arl(h1, g1, lp, mu, k, r=r)
  dratio1 <- ( ARLm - ARLp ) / ( 2*eps )
  if ( OUTPUT ) cat(paste("0\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
  
  if ( dratio1 > 0 ) {
    while ( dratio1 > 0 ) {
      h2 <- h1
      dratio2 <- dratio1
      h1 <- h1 * 1.2
      g1 <- be.get.gamma(h1, L0, lambda, mu, k, r=r)
      ARLm <- BE2.arl(h1, g1, lm, mu, k, r=r)
      ARLp <- BE2.arl(h1, g1, lp, mu, k, r=r)
      dratio1 <- ( ARLm - ARLp ) / ( 2*eps )
      if ( OUTPUT ) cat(paste("1\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
    }
  } else {
    while ( dratio1 <= 0 ) {
      h2 <- h1
      dratio2 <- dratio1
      h1 <- h1 * 0.8
      g1 <- be.get.gamma(h1, L0, lambda, mu, k, r=r)
      ARLm <- BE2.arl(h1, g1, lm, mu, k, r=r)
      ARLp <- BE2.arl(h1, g1, lp, mu, k, r=r)
      dratio1 <- ( ARLm - ARLp ) / ( 2*eps )
      if ( OUTPUT ) cat(paste("1\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
    }
  }
  
  h.error <- 1
  dr.error <- 1
  while ( abs(h.error) > 1e-10 & abs(dr.error) > 1e-8 ) {
    h3 <- h1 - dratio1 / (dratio2 - dratio1) * (h2 - h1)
    g3 <- be.get.gamma(h3, L0, lambda, mu, k, r=r)
    ARLm <- BE2.arl(h3, g3, lm, mu, k, r=r)
    ARLp <- BE2.arl(h3, g3, lp, mu, k, r=r)
    dratio3 <- ( ARLm - ARLp ) / ( 2*eps )
    if ( OUTPUT ) cat(paste("2\th3 =", h3, "\tg3 =", g3, "\tdratio1 =", dratio3, "\n"))
    h1 <- h2
    h2 <- h3
    dratio1 <- dratio2
    dratio2 <- dratio3
    h.error <- h2 - h1
    dr.error <- dratio2 - dratio1
  }
  data.frame(h=h3, gamma=g3)
}


# calculate (h, gamma) so that (i) ARL = L0 and (ii) ARL as function of mu takes max at given mu value, collocation based
get.gammaAh.mu <- function(L0, lambda, mu, k, z0=0, r=40, qm=30, eps=1e-6, OUTPUT=FALSE) {
  mm <- mu - eps
  mp <- mu + eps
  
  h1 <- get.h(2*L0, lambda, mu, z0=z0, r=r, qm=qm)
  g1 <- get.gamma(h1, L0, lambda, mu, k, z0=z0, r=r, qm=qm)
  ARLm <- coll.arl(h1, g1, lambda, mm, k, z0=z0, r=r, qm=qm)
  ARLp <- coll.arl(h1, g1, lambda, mp, k, z0=z0, r=r, qm=qm)
  dratio1 <- ( ARLm - ARLp ) / ( 2*eps )
  if ( OUTPUT ) cat(paste("0\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
  
  if ( dratio1 > 0 ) {
    while ( dratio1 > 0 ) {
      h2 <- h1
      dratio2 <- dratio1
      h1 <- h1 * 0.8
      g1 <- get.gamma(h1, L0, lambda, mu, k, z0=z0, r=r, qm=qm)
      ARLm <- coll.arl(h1, g1, lambda, mm, k, z0=z0, r=r, qm=qm)
      ARLp <- coll.arl(h1, g1, lambda, mp, k, z0=z0, r=r, qm=qm)
      dratio1 <- ( ARLm - ARLp ) / ( 2*eps )
      if ( OUTPUT ) cat(paste("1\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
    }
  } else {
    while ( dratio1 <= 0 ) {
      h2 <- h1
      dratio2 <- dratio1
      h1 <- h1 * 1.2
      g1 <- get.gamma(h1, L0, lambda, mu, k, z0=z0, r=r, qm=qm)
      ARLm <- coll.arl(h1, g1, lambda, mm, k, z0=z0, r=r, qm=qm)
      ARLp <- coll.arl(h1, g1, lambda, mp, k, z0=z0, r=r, qm=qm)
      dratio1 <- ( ARLm - ARLp ) / ( 2*eps )
      if ( OUTPUT ) cat(paste("1\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
    }
  }
  
  h.error <- 1
  dr.error <- 1
  while ( abs(h.error) > 1e-10 & abs(dr.error) > 1e-8 ) {
    h3 <- h1 - dratio1 / (dratio2 - dratio1) * (h2 - h1)
    g3 <- get.gamma(h3, L0, lambda, mu, k, z0=z0, r=r, qm=qm)
    ARLm <- coll.arl(h3, g3, lambda, mm, k, z0=z0, r=r, qm=qm)
    ARLp <- coll.arl(h3, g3, lambda, mp, k, z0=z0, r=r, qm=qm)
    dratio3 <- ( ARLm - ARLp ) / ( 2*eps )
    if ( OUTPUT ) cat(paste("2\th3 =", h3, "\tg3 =", g3, "\tdratio1 =", dratio3, "\n"))
    h1 <- h2
    h2 <- h3
    dratio1 <- dratio2
    dratio2 <- dratio3
    h.error <- h2 - h1
    dr.error <- dratio2 - dratio1
  }
  data.frame(h=h3, gamma=g3)
}


# calculate (h, gamma) so that (i) ARL = L0 and (ii) ARL as function of mu takes max at given mu value, MC based
be.get.gammaAh.mu <- function(L0, lambda, mu, k, z0=0, r=100, eps=1e-6, OUTPUT=FALSE) {
  mm <- mu - eps
  mp <- mu + eps
  
  h1 <- be.get.h(2*L0, lambda, mu, k, r=r)
  g1 <- be.get.gamma(h1, L0, lambda, mu, k, r=r)
  ARLm <- BE2.arl(h1, g1, lambda, mm, k, r=r)
  ARLp <- BE2.arl(h1, g1, lambda, mp, k, r=r)
  dratio1 <- ( ARLm - ARLp ) / ( 2*eps )
  if ( OUTPUT ) cat(paste("0\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
  
  if ( dratio1 > 0 ) {
    while ( dratio1 > 0 ) {
      h2 <- h1
      dratio2 <- dratio1
      h1 <- h1 * 0.8
      g1 <- be.get.gamma(h1, L0, lambda, mu, k, r=r)
      ARLm <- BE2.arl(h1, g1, lambda, mm, k, r=r)
      ARLp <- BE2.arl(h1, g1, lambda, mp, k, r=r)
      dratio1 <- ( ARLm - ARLp ) / ( 2*eps )
      if ( OUTPUT ) cat(paste("1\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
    }
  } else {
    while ( dratio1 <= 0 ) {
      h2 <- h1
      dratio2 <- dratio1
      h1 <- h1 * 1.2
      g1 <- be.get.gamma(h1, L0, lambda, mu, k, r=r)
      ARLm <- BE2.arl(h1, g1, lambda, mm, k, r=r)
      ARLp <- BE2.arl(h1, g1, lambda, mp, k, r=r)
      dratio1 <- ( ARLm - ARLp ) / ( 2*eps )
      if ( OUTPUT ) cat(paste("1\th1 =", h1, "\tg1 =", g1, "\tdratio1 =", dratio1, "\n"))
    }
  }
  
  h.error <- 1
  dr.error <- 1
  while ( abs(h.error) > 1e-10 & abs(dr.error) > 1e-8 ) {
    h3 <- h1 - dratio1 / (dratio2 - dratio1) * (h2 - h1)
    g3 <- get.gamma(h3, L0, lambda, mu, k, z0=z0, r=r)
    ARLm <- BE2.arl(h3, g3, lambda, mm, k, r=r)
    ARLp <- BE2.arl(h3, g3, lambda, mp, k, r=r)
    dratio3 <- ( ARLm - ARLp ) / ( 2*eps )
    if ( OUTPUT ) cat(paste("2\th3 =", h3, "\tg3 =", g3, "\tdratio1 =", dratio3, "\n"))
    h1 <- h2
    h2 <- h3
    dratio1 <- dratio2
    dratio2 <- dratio3
    h.error <- h2 - h1
    dr.error <- dratio2 - dratio1
  }
  data.frame(h=h3, gamma=g3)
}

#########################################
#########################################
#########################################

# Load necessary libraries
library(shiny)
library(ggplot2)
library(shinycssloaders)


# Define UI
ui <- fluidPage(
  titlePanel("ARL-unbiased Wn chart - Fixed Service Rate"),
  
  sidebarLayout(
    sidebarPanel(
      radioButtons("distribution", "Select queue:",
                   choices = c("Ek/M/1", "M/Ek/1")),
      conditionalPanel(
        condition = "input.distribution == 'M/Ek/1'",
        sliderInput("lambda", HTML("Target Traffic Intensity (&rho;0):"), min = 0.01, max = 0.99, value = 0.1, step = 0.01),
        sliderInput("L0", "Target ARL:", min = 100, max = 1000, value = 200, step = 1),
        sliderInput("k", "Erlang Parameter (k):", min = 1, max = 100, value = 1, step = 1)
      ),
      conditionalPanel(
        condition = "input.distribution == 'Ek/M/1'",
        sliderInput("lambda", HTML("Target Traffic Intensity (&rho;0):"), min = 0.01, max = 0.99, value = 0.1, step = 0.01),
        sliderInput("L0", "Target ARL:", min = 100, max = 1000, value = 200, step = 1),
        sliderInput("k", "Erlang Parameter (k):", min = 1, max = 100, value = 1, step = 1)
      ),
      downloadButton("downloadPlot", "Download Plot")
    ),
    
    mainPanel(
      withSpinner(plotOutput("ARLPlot"), type = 8, color = "#0D6EFD")
    )
  )
)


server <- function(input, output) {
  
  # Store the calculated plot data for re-use in downloadHandler
  rv <- reactiveValues()
  
  output$ARLPlot <- renderPlot({
    # Create a progress object
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    
    distribution <- input$distribution
    
    if (distribution == "M/Ek/1") {
      progress$set(message = "Calculating", value = 0)
      progress$inc(0.2, message = "Initializing variables")
      
      lambda0 <- input$lambda
      mu <- 1
      L0 <- input$L0
      k <- input$k
      
      progress$inc(0.3, message = "Calculating limits")
      calibration_results_m <- get.gammaAh.lambda_m(L0, lambda0, mu, k)
      h <- calibration_results_m$h
      gamma <- calibration_results_m$gamma
      
      progress$inc(0.3, message = "Generating plot data")
      lambda_seq <- seq(0.01, 0.99, by = 0.01)
      ARL_values <- sapply(lambda_seq, function(lambda) coll.arl_m(h, gamma, lambda, mu, k))
      
      # Cache results in reactiveValues for downloading
      rv$lambda_seq <- lambda_seq
      rv$ARL_values <- ARL_values
      rv$h <- h
      rv$gamma <- gamma
      
      progress$inc(0.2, message = "Rendering plot")
      legend_pos <- if (lambda0 > 0.5) "topleft" else "topright"
      
      plot(lambda_seq, ARL_values, type = 'l', col = 'black', xlab = "ρ", ylab = "ARL")
      abline(v=lambda0, h=L0, lty=4, col="grey")
      legend(legend_pos, legend = c(paste("UCL =", round(h, 6)), paste("γL =", round(gamma, 6))), cex = 0.8, inset = c(0.025, 0.1))
      
    } else if (distribution == "Ek/M/1") {
      progress$set(message = "Calculating", value = 0)
      progress$inc(0.2, message = "Initializing variables")
      
      lambda0 <- input$lambda
      mu <- 1
      L0 <- input$L0
      k <- input$k
      
      progress$inc(0.3, message = "Calculating limits")
      calibration_results <- get.gammaAh.lambda(L0, lambda0, mu, k)
      h <- calibration_results$h
      gamma <- calibration_results$gamma
      
      progress$inc(0.3, message = "Generating plot data")
      lambda_seq <- seq(0.01, 0.99, by = 0.01)
      ARL_values <- sapply(lambda_seq, function(lambda) coll.arl(h, gamma, lambda, mu, k))
      
      # Cache results in reactiveValues for downloading
      rv$lambda_seq <- lambda_seq
      rv$ARL_values <- ARL_values
      rv$h <- h
      rv$gamma <- gamma
      
      progress$inc(0.2, message = "Rendering plot")
      legend_pos <- if (lambda0 > 0.5) "topleft" else "topright"
      
      plot(lambda_seq, ARL_values, type = 'l', col = 'black', xlab = "ρ", ylab = "ARL")
      abline(v=lambda0, h=L0, lty=4, col="grey")
      legend(legend_pos, legend = c(paste("UCL =", round(h, 6)), paste("γL =", round(gamma, 6))), cex = 0.8, inset = c(0.025, 0.1))
    }
  })
  
  # Add downloadHandler for downloading the plot
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("ARLPlot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      png(file)
      # Reuse the cached plot data from reactiveValues
      lambda_seq <- rv$lambda_seq
      ARL_values <- rv$ARL_values
      h <- rv$h
      gamma <- rv$gamma
      
      legend_pos <- if (input$lambda > 0.5) "topleft" else "topright"
      
      plot(lambda_seq, ARL_values, type = 'l', col = 'black', xlab = "ρ", ylab = "ARL")
      abline(v=input$lambda, h=input$L0, lty=4, col="grey")
      legend(legend_pos, legend = c(paste("UCL =", round(h, 6)), paste("γL =", round(gamma, 6))), cex = 0.8, inset = c(0.025, 0.1))
      
      dev.off()
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
