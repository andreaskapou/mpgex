create_diff_bpr_data <- function(N = 300, pi.c = c(0.45, 0.35, 0.2), max_L = 30,
                            xmin = -100, xmax=100, fmin = -1, fmax = 1){
  # Create a list to store data for each methylation region for both control and
  # treated samples
  X       <- list(control   = list(),
                  treatment = list())
  # A list for storing corresponding gene expression data for both control and
  # treated samples
  Y       <- list(control   = vector(mode = "numeric", length = N),
                  treatment = vector(mode = "numeric", length = N))

  # For each of the N objects
  for (i in 1:N){
    # L is the number of CpGs found in the ith region
    L <- rbinom(n = 1, size = max_L, prob = .8)
    X$control[[i]] <- matrix(0, nrow = L, ncol = 3)
    X$treatment[[i]] <- X$control[[i]]
    # Randomly sample locations for the CpGs
    obs <- sort(sample(xmin:xmax, L))
    # Scale them, so the data lie in the (fmin, fmax) range
    X$control[[i]][ ,1] <- minmax_scaling(data = obs,
                                          xmin = xmin,
                                          xmax = xmax,
                                          fmin = fmin,
                                          fmax = fmax)
    X$treatment[[i]][ ,1] <- X$control[[i]][ ,1]


    if (i < N * pi.c[1]){   # First methylation profile
      lb <- round(L / 4)

      # Control
      X$control[[i]][1:lb,2] <- rbinom(lb, 20, .9)
      repeat{
        X$control[[i]][1:lb,3] <- rbinom(lb, 14, .9)
        if(all(X$control[[i]][1:lb,2] > X$control[[i]][1:lb,3]))
          break
      }

      X$control[[i]][(lb + 1):L,2] <- rbinom(L - lb, 20, .9)
      repeat{
        X$control[[i]][(lb + 1):L,3] <- rbinom(L - lb, 2, .9)
        if (all(X$control[[i]][(lb + 1):L,2] > X$control[[i]][(lb + 1):L,3]))
          break
      }
      Y$control[i] <- rpois(1, lambda=200)

      # Treatment
      X$treatment[[i]][1:lb,2] <- rbinom(lb, 20, .9)
      repeat{
        X$treatment[[i]][1:lb,3] <- rbinom(lb, 2, .9)
        if(all(X$treatment[[i]][1:lb,2] > X$treatment[[i]][1:lb,3]))
          break
      }

      X$treatment[[i]][(lb + 1):L,2] <- rbinom(L - lb, 20, .9)
      repeat{
        X$treatment[[i]][(lb + 1):L,3] <- rbinom(L - lb, 10, .9)
        if (all(X$treatment[[i]][(lb + 1):L,2] > X$treatment[[i]][(lb + 1):L,3]))
          break
      }
      Y$treatment[i] <- rpois(1, lambda=100)


    }else if (i < (N * pi.c[2] + N * pi.c[1])){ # Second methylation profile
      lb <- round(L / 1.5)

      # Control
      X$control[[i]][1:lb,2] <- rbinom(lb, 20, .9)
      repeat{
        X$control[[i]][1:lb,3] <- rbinom(lb, 2, .8)
        if(all(X$control[[i]][1:lb,2] > X$control[[i]][1:lb,3]))
          break
      }

      X$control[[i]][(lb + 1):L,2] <- rbinom(L - lb, 20, .9)
      repeat{
        X$control[[i]][(lb + 1):L,3] <- rbinom(L-lb, 14, .9)
        if (all(X$control[[i]][(lb + 1):L,2] > X$control[[i]][(lb + 1):L,3]))
          break
      }
      Y$control[i] <- rpois(1, lambda=100)

      # Treatment
      X$treatment[[i]][1:lb,2] <- rbinom(lb, 20, .9)
      repeat{
        X$treatment[[i]][1:lb,3] <- rbinom(lb, 12, .8)
        if(all(X$treatment[[i]][1:lb,2] > X$treatment[[i]][1:lb,3]))
          break
      }

      X$treatment[[i]][(lb + 1):L,2] <- rbinom(L - lb, 20, .9)
      repeat{
        X$treatment[[i]][(lb + 1):L,3] <- rbinom(L-lb, 4, .9)
        if (all(X$treatment[[i]][(lb + 1):L,2] > X$treatment[[i]][(lb + 1):L,3]))
          break
      }
      Y$treatment[i] <- rpois(1, lambda=200)


    }else{                  # Third methylation profile
      lb <- round(L / 2.5)
      mb <- round(L / 3.5)

      # Control
      X$control[[i]][1:lb,2] <- rbinom(lb, 20, .9)
      repeat{
        X$control[[i]][1:lb,3] <- rbinom(lb, 2, .9)
        if(all(X$control[[i]][1:lb,2] > X$control[[i]][1:lb,3]))
          break
      }

      X$control[[i]][(lb + 1):(lb + mb),2] <- rbinom(mb, 20, .9)
      repeat{
        X$control[[i]][(lb + 1):(lb + mb),3] <- rbinom(mb, 14, .9)
        if (all(X$control[[i]][(lb + 1):(lb + mb),2] > X$control[[i]][(lb + 1):(lb + mb),3]))
          break
      }

      X$control[[i]][(lb + 1 + mb):L,2] <- rbinom(L - mb - lb, 20, .9)
      repeat{
        X$control[[i]][(lb + 1 + mb):L,3] <- rbinom(L - mb - lb, 2, .9)
        if (all(X$control[[i]][(lb + 1 + mb):L,2] > X$control[[i]][(lb + 1 + mb):L],3))
          break
      }
      Y$control[i] <- rpois(1, lambda=60)

      # Treatment
      X$treatment[[i]][1:lb,2] <- rbinom(lb, 20, .9)
      repeat{
        X$treatment[[i]][1:lb,3] <- rbinom(lb, 2, .9)
        if(all(X$treatment[[i]][1:lb,2] > X$treatment[[i]][1:lb,3]))
          break
      }

      X$treatment[[i]][(lb + 1):(lb + mb),2] <- rbinom(mb, 20, .9)
      repeat{
        X$treatment[[i]][(lb + 1):(lb + mb),3] <- rbinom(mb, 14, .9)
        if (all(X$treatment[[i]][(lb + 1):(lb + mb),2] > X$treatment[[i]][(lb + 1):(lb + mb),3]))
          break
      }

      X$treatment[[i]][(lb + 1 + mb):L,2] <- rbinom(L - mb - lb, 20, .9)
      repeat{
        X$treatment[[i]][(lb + 1 + mb):L,3] <- rbinom(L - mb - lb, 2, .9)
        if (all(X$treatment[[i]][(lb + 1 + mb):L,2] > X$treatment[[i]][(lb + 1 + mb):L],3))
          break
      }
      Y$treatment[i] <- rpois(1, lambda=60)
    }
  }
  return(list(X = X, Y = Y))
}

set.seed(1)
bpr <- create_diff_bpr_data(N=600)
bpr_control_data <- bpr$X$control
bpr_treatment_data <- bpr$X$treatment
gex_control_data <- bpr$Y$control
gex_treatment_data <- bpr$Y$treatment
devtools::use_data(bpr_control_data,
                   bpr_treatment_data,
                   gex_control_data,
                   gex_treatment_data, overwrite = TRUE)
