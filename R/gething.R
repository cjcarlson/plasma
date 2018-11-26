#' Title
#' @description
#' Text
#'
#' @param value text
#'
#' @keywords text
#'
#' @import raster

gething <- function(temps,malaria,burnin=1,maxtime=1,tau=30,phi=NA,Tmin=NA) {

  # What kind of malaria are we working with here?

  # temps <- stack(lapply(1:(180*12), function(i) T))

  blank <- temps[[1]]*0

  if(malaria=='Pf') {
    phi <- 111
    Tmin <- 16
  }

  if(malaria=='Pv') {
    phi <-105
    Tmin <- 14.5
  }

  if(malaria=='custom') {
    phi <- phi
    Tmin <- Tmin
  }

  # initiate the blank rasters

  pT <- blank
  dT <- blank
  M <- blank
  Z <- blank

  # initiate blank raster stacks

  Y <- stack(lapply(1:(tau*12), function(i) blank))
  S <- Y
  I <- Y
  Zdummy <- Y

  # This is the burn-in time loop
  for (i in 1:(12*60*burnin)) {

    T <- temps[[i]]
    pT <- exp(-1/(-4.4+1.31*T-0.03*T^2))
    dT <- (T-Tmin)/12
    M <- 1 + M*pT

    for (k in nlayers(S):2) {
      S[[k]] <- S[[k-1]] + dT
    }

    for (k in 1:(nlayers(I)-1)) {
      I[[k]] <- 1-sign(sign(S[[k]]-phi)+1)
    }

    Y[[1]] <- M

    for (k in (nlayers(Y)-1):1) {
      Y[[k+1]] <- Y[[k]]*pT*I[[k]]
    }

    for (k in 1:nlayers(Y)) {
      Zdummy[[k]] <- Y[[k]]*pT*(1-I[[k]])
    }

    Z <- Z*pT + sum(Zdummy)
    print(i)
    plot(Z)

  }

  print("burn-in complete")

  # This is the actual time loop

  Zstack <- stack() # stores values

  for (i in 1:(12*60*burnin)) {

    T <- temps[[i]]
    pT <- exp(-1/(-4.4+1.31*T-0.03*T^2))
    dT <- (T-Tmin)/12
    M <- 1 + M*pT

    for (k in nlayers(S):2) {
      S[[k]] <- S[[k-1]] + dT
    }

    for (k in 1:(nlayers(I)-1)) {
      I[[k]] <- 1-sign(sign(S[[k]]-phi)+1)
    }

    Y[[1]] <- M

    for (k in (length(Y)-1):1) {
      Y[[k+1]] <- Y[[k]]*pT*I[[k]]
    }

    for (k in 1:nlayers(Y)) {
      Zdummy[[k]] <- Y[[k]]*pT*(1-I[[k]])
    }

    Z <- Z*pT + sum(Zdummy)
    Zstack <- stack(Zstack, Z)

  }

  return(Zstack)

}
