library(rstan)
library(Matrix)
library(spdep)

# Load data
load("ck_data.RData")
data <- ck_data
period <- with(data, congress - 45)
N <- nrow(data)
C <- length(unique(data$congress))

make_rw_matrix <- function(dim, which = c("penalty", "adjacency", "degree"), 
                           type = c("breslow", "rw1", "rw2")) {
  which <- match.arg(which)
  type <- match.arg(type)
  
  require(Matrix)
  require(spdep)
  omega <- 0.99
  if (type == "breslow") {
    diag <- c(1,5, rep(6, dim-4), 5, 1)
    band1 <- c(-2, rep(-4, dim-3), -2)
    band2 <- rep(1, dim-2)
    Diags <- list(diag, band1, band2)
    k <- c(0,1,2)
    R <- bandSparse(n = dim, k = k, diag = Diags, symm=TRUE)
    return(as.matrix(R))
  }
  if (type == "rw1") {
    nb <- cell2nb(1, dim)
  }
  if (type == "rw2") {
    nb <- cell2nb(1, dim)
    nb[[1]] <- c(nb[[1]], 3)
    nb[[2]] <- c(nb[[2]], 4)
    for(i in 3:(dim-2)) {
      nb[[i]] <- c(i - 2, nb[[i]], i + 2)
    }
    nb[[dim-1]] <- c(nb[[dim-1]], dim-3)
    nb[[dim]] <- c(nb[[dim]], dim-2) 
  }
  A <- mat.or.vec(dim, dim)
  for(i in 1:dim) {
    for(j in 1:dim) {
      if (i %in% nb[[j]]) A[i,j] <- 1
    }
  }
  if (which == "adjacency") return(A)
  if (which == "degree") return(diag(rowSums(A)))
  diag(rowSums(A)) - omega * A
}

# Penalty matrix 
P <- make_rw_matrix(dim = C, which = "penalty", type = "rw1")

# Data for Stan
stan_data <- list(N = N,
                  C = C,
                  Pinverse = solve(P),
                  cong = period,
                  nVotes = data$nvotes, 
                  nWins = data$majps, 
                  lvRatio = data$lnmajvavg)


fit_compile <- stan(file = "code_for_paper.stan", 
                    data = stan_data, chains = 1, iter = 10)

stanfit <- stan(fit = fit_compile, 
                data = stan_data, 
                chains = 8, 
                iter = 2000)

