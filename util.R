library(lubridate)

read_precip_data <- function(file) {
  df <- read.csv(file, skip = 3, na.strings = "-99")
  colnames(df) <- c("Date", "Value")
  return(df)
}

convert_to_date <- function(numeric_date) {
  year <- as.integer(numeric_date / 100)
  month <- sprintf("%02d", numeric_date %% 100)
  return(paste0(year, "-", month, "-01"))
}

station_data_get = function (ss_ls) {
  p = length(ss_ls)
  for (i in 1:p) {
    state_name = names(ss_ls)[i]
    for (j in 1:ss_ls[[i]]) {
      dir_name = paste0("./precipitation/", state_name, j, ".csv")
      temp = read_precip_data(dir_name)
      if (i == 1 && j == 1) {
        n = nrow(temp)
        Dates = convert_to_date(temp[,1])
        res = matrix(0, nrow = n, ncol = ss_ls[[1]])
      } 
      res[,j] = temp[,2]
    }
  }
  
  return (list("dta" = res, "Dates" = Dates))
}

yearly_obs = function (x) {
  n = length(x)
  n_full_yr = n %/% 12
  x = x[1:(n_full_yr * 12)]
  res = matrix(x, ncol = 12, byrow = T)
  
  return(res)
}

#' exponential map on the sphere
#' 
#' @param x an array of tangent vectors
#' @param mu base point 
Exp_sphere = function (x, mu) {
  
  if (is.matrix(x)) {
    
    if (sum(abs(x)) == 0) {
      return (matrix(mu, nrow = nrow(x), ncol = ncol(x), byrow = T))
    }
    
    x_norm = sqrt(rowSums(x^2))
    x_norm[which(x_norm == 0)] = 1
    std_x = x / x_norm
    res = outer(cos(x_norm), mu) + sin(x_norm) * std_x
    
    res = res / sqrt(rowSums(res^2)) # normalize again to avoid numerical instability
  } else {
    
    if (sum(abs(x)) == 0) {
      return (mu)
    }
    
    x_norm = sqrt(sum(x^2))
    res = cos(x_norm) * mu + sin(x_norm) * x / x_norm
    
    res = res / sqrt(sum(res^2))
  }
  
  return (res)
}

#' logarithmic map for the sphere
#' 
#' @param x an array of points on the sphere
#' @param mu a base point
Log_sphere = function (x, mu) {
  
  if (is.matrix(x)) {
    w = x - matrix(mu, nrow = nrow(x), ncol = ncol(x), byrow = T)
    Proj = w - outer(c(w %*% mu), mu)
    Proj = Proj / sqrt(rowSums(Proj^2))
    
    res = acos(c(x %*% mu)) * Proj
    
    if (any(rowSums(abs(w)) < 1e-7)) {
      res[which(rowSums(abs(w)) < 1e-7),] = 0
    }
    
    
  } else {
    w = x - mu
    Proj = w - mu * c(t(mu) %*% w)
    
    if (sum(abs(w)) < 1e-7) {
      res = rep(0, length(mu))
    } else {
      res = acos(c(x %*% mu)) * Proj / sqrt(sum(Proj^2))
    }
  }
  
  return (res)
}

#' computes the Frechet mean on the sphere
#' 
#' @param x (n by q) array of data
#' 
mean_on_sphere = function (x, tau = 0.1, tol = 1e-8, max.iter = 1000, verbose = FALSE) {
  
  n = nrow(x)
  mu = x[sample(n, 1),]
  for (i in 1:max.iter) {
    grad = colMeans(Log_sphere(x, mu))
    mu_new = Exp_sphere(tau * grad, mu)
    
    temp = c(x %*% mu_new / sqrt(rowSums(x^2)))
    if (any(temp > 1.01) || any(temp < -1.01)) {
      stop("mean_on_sphere: something must be wrong")
    } else {
      temp = pmin(pmax(temp, -1), 1)
    }
    loss = mean(acos(temp))
    if (i > 1 && (loss_old - loss < tol)) {
      mu = mu_new
      break
    }
    if (verbose) {
      cat("mean_on_sphere: iter", i, "; loss", round(loss, 4), "\n")
    }
    mu = mu_new
    loss_old = loss
  }
  
  return (mu)
}

LYB_fm = function (x, r, h, demean = TRUE) {
  # x: (n by d) observation matrix
  # r: number of factors
  # h: number of lags to use
  # Estimate the factor model of Lam, Yao, and Bathia
  
  if (demean) {
    mean = colMeans(x)
    x = t(t(x) - colMeans(x))
  }
  
  # Compute the auxiliary positive definite matrix
  pd = 0
  H = h + 1
  n = nrow(x)
  for (i in 1:h) {
    temp = t(x[H:n,]) %*% x[(H - i):(n - i),] / n
    pd = pd + temp %*% t(temp)
  }
  
  # Eigenanalysis
  model = eigen(pd)
  Evec = model$vectors
  V = Evec[,1:r] # (d by r)
  evals = model$values
  
  # Extract factors and residuals
  f_hat = t(t(V) %*% t(x)) # (n by r)
  e_hat = x - f_hat %*% t(V)
  
  # Estimate the number of factors
  ratios = evals[2:r] / evals[1:(r - 1)]
  r_hat = which.min(ratios)
  
  return (list("V" = V, "f_hat" = f_hat, "e_hat" = e_hat, 
               "fitted.val" = f_hat %*% t(V), "r_hat" = r_hat,
               "mean" = mean))
}

predict_fm = function (V, mu, new_x) {
  
  if (is.matrix(new_x)) {
    z_temp = t(t(new_x) - mu)
    Factor = z_temp %*% V
    z_hat = Factor %*% t(V)
    x_hat = t(t(z_hat) + mu)
  } else if (is.vector(new_x)) {
    z_temp = new_x - mu
    Factor = c(t(z_temp) %*% V)
    z_hat = V %*% Factor
    x_hat = mu + z_hat
  } else {
    stop("predict_fm: new_x must be matrix or vector")
  }
  
  return (x_hat)
}

tan_basis_sphere = function (mu) {
  n = length(mu)
  U = svd(mu, nu = n, nv = n)$u
  
  return(U[,-1])
}

#' Helper that slices the correct indices
q_index = function (qs, j, type = "ambient") {
  if (type == "ambient") {
    NULL
  } else if (type == "intrinsic") {
    qs = qs - 1
  } else {
    stop("q_index: unsupported slicing type")
  }
  
  if (j == 1) {
    res = 1:qs[j]
  } else {
    res = sum(qs[1:(j - 1)]) + c(1:qs[j])
  }
  
  return (res)
}

#' estimate the linear factor model for the (product of spheres)
#' 
lfm_sphere = function (x, r, h = 6) {
  d = length(x)
  n = dim(x[[1]])[1]
  qs = rep(NA, d)
  for (j in 1:d) {
    qs[j] = dim(x[[j]])[2]
  }
  
  # flatten the list
  x_train = NULL
  for (j in 1:d) {
    x_train = cbind(x_train, x[[j]])
  }
  
  # estimate the factor model
  model = LYB_fm(x_train, r = r, h = h)
  
  return (list("A" = model$V, "f_hat" = model$f_hat, "factor_model" = model,
               "r_hat" = model$r_hat))
}

#' predict from the lfm_sphere output
#' the predictions live on the sphere
#' 
predict_lfm = function (x_test, model) {
  d = length(x_test)
  n = dim(x_test[[1]])[1]
  qs = rep(0, d)
  for (i in 1:d) {
    qs[i] = dim(x_test[[i]])[2]
  }
  
  V = model$A
  z_mean = model$factor_model$mean
  if (is.vector(V)) {
    r = 1
  } else {
    r = dim(V)[2]
  }
  
  # flatten data
  x_test_flat = NULL
  for (i in 1:d) {
    x_test_flat = cbind(x_test_flat, x_test[[i]])
  }
  
  # make predictions
  z_hat = predict_fm(V, z_mean, x_test_flat)
  x_hat = vector("list", d)
  for (i in 1:d) {
    x_hat[[i]] = z_hat[,q_index(qs, i, "ambient")]
  }
  
  return (x_hat)
}

#' estimate the RFM for the (product of spheres)
#' 
#' @param x a list of (n by q_j) arrays, j = 1,2,...,d
#' 
#' @export
rfm_sphere = function (x, r, h = 6, tau = 0.5, max.iter = 100) {
  d = length(x)
  n = dim(x[[1]])[1]
  qs = rep(NA, d)
  for (j in 1:d) {
    qs[j] = dim(x[[j]])[2]
  }
  
  # estimate the Frechet mean
  mus = vector("list", length = d)
  for (j in 1:d) {
    mus[[j]] = mean_on_sphere(x[[j]], tau = tau, max.iter = max.iter,
                              verbose = FALSE) 
  }
  
  # construct log-mapped data
  log_x_vec = NULL
  bases = vector("list", d)
  for (j in 1:d) {
    bases[[j]] = tan_basis_sphere(mus[[j]])
    
    log_x = Log_sphere(x[[j]], mus[[j]])
    temp = log_x %*% bases[[j]]
    
    log_x_vec = cbind(log_x_vec, temp)
  }
  if (ncol(log_x_vec) != (sum(qs) - d)) {
    stop("rfm_sphere: incorrect log_x_vec dimension") # for prototyping
  }
  
  # estimate factor model
  model = LYB_fm(log_x_vec, r = r, h = h)
  
  return (list("A" = model$V, "f_hat" = model$f_hat, "E" = bases,
               "mu_hat" = mus, "factor_model" = model,
               "r_hat" = model$r_hat))  
  
}

#' predict from an rfm_sphere output
#' the predictions live on the sphere
#' 
predict_rfm = function (x_test, model) {
  d = length(x_test)
  n = dim(x_test[[1]])[1]
  qs = rep(0, d)
  
  for (i in 1:d) {
    qs[i] = dim(x_test[[i]])[2]
  }
  
  V = model$A
  z_mean = model$factor_model$mean
  if (is.vector(V)) {
    r = 1
  } else {
    r = dim(V)[2]
  }
  
  mu_hat = model$mu_hat
  E = model$E
  factor_model = model$factor_model
  
  # Construct log mapped data
  log_x_vec = NULL
  for (i in 1:d) {
    log_x = Log_sphere(x_test[[i]], mu_hat[[i]])
    temp = log_x %*% E[[i]]
    
    log_x_vec = cbind(log_x_vec, temp)
  }
  
  # make predictions
  x_hat = vector("list", d)
  z_hat = predict_fm(V, z_mean, log_x_vec)
  
  for (j in 1:d) {
    idx = q_index(qs, j, "intrinsic")
    
    temp = z_hat[,idx]
    x_hat[[j]] = Exp_sphere(temp %*% t(E[[j]]), mu_hat[[j]])
  }
  
  return (x_hat)
}

#' computes Log Ratio transformation for compositional data
#' 
#' @param x (n by q) matrix
#' @param g index for the reference group (if "alr" is used)
#' @param type type of transformations ("alr", "clr", "ilr")
#' @param inv whether to compute the inverse transformation. Default is FALSE
#' @param lr_obj log ratio object for meta data, only needed for inverse transforms
#' 
LR_trans = function (x, type = "alr", inv = FALSE, g = NULL, lr_obj = NULL) {
  n = nrow(x)
  q = ncol(x)
  
  if (inv == FALSE) {
    if (type == "alr") {
      lr_obj = alr(acomp(x), ivar = g)
    } else if (type == "clr") {
      lr_obj = clr(acomp(x))
    } else if (type == "ilr") {
      lr_obj = ilr(acomp(x))
    } else {
      stop("LR_trans: transformation type not supported")
    }
    
    lr_clean = matrix(as.numeric(lr_obj), nrow = n, ncol = q - 1)
    
    return (list("x" = lr_clean, "x.lr.obj" = lr_obj))
  } else {
    if (type == "alr") {
      inv_obj = alrInv(x, orig = lr_obj)
    } else if (type == "clr") {
      inv_obj = clrInv(x, orig = lr_obj)
    } else if (type == "ilr") {
      inv_obj = ilrInv(x, orig = lr_obj)
    } else {
      stop("LR_trans: transformation type not supported")
    }
    
    lr_clean = matrix(as.numeric(inv_obj), nrow = n, ncol = q + 1)
    
    return (list("x" = lr_clean, "x.lr.obj" = inv_obj))
  }
}

#' estimate the factor model after log ratio transformations
#' 
#' @param x a list of (n by q_j) arrays, j = 1,2,...,d
#' @param r number of factors
#' @param h lags used in factor estimation
#' @param type type of log ratio transformation ("alr", "clr", "ilr")
#' @param reference_group a vector of indices (for each j) which are used as the 
#'                        reference group in additive log ratio transformation
#'                        
logR_fm = function (x, r, h = 6, type = "alr", reference_group = NULL) {
  d = length(x)
  n = dim(x[[1]])[1]
  qs = rep(NA, d)
  for (i in 1:d) {
    qs[i] = dim(x[[i]])[2]
  }
  if (type == "alr") {
    if (is.null(reference_group)) {
      ref_g = qs
    } else {
      ref_g = reference_group
    } 
  }
  
  x_lr = NULL
  x_lr_obj = vector("list", length = d)
  for (i in 1:d) {
    temp = LR_trans(x[[i]], type = type, g = ref_g[i])
    x_lr = cbind(x_lr, temp$x)
    x_lr_obj = temp$x.lr.obj
  }
  
  model = LYB_fm(x_lr, r = r, h = h)
  
  return (list("A" = model$V, "f_hat" = model$f_hat, "x" = x_lr,
               "x.lr.obj" = x_lr_obj, "factor_model" = model,
               "r_hat" = model$r_hat, "ref_g" = ref_g))
  
}

#' predict from an logR_fm output
#' 
predict_logR = function (x_test, model) {
  
}

#' predict from an alr_fm output
#' 
predict_alr = function (x_test, model) {
  d = length(x_test)
  n = dim(x_test[[1]])[1]
  qs = rep(NA, d)
  for (i in 1:d) {
    qs[i] = dim(x_test[[i]])[2]
  }
  ref_g = model$ref_g
  V = model$A
  z_mean = model$factor_model$mean
  
  x_alr = NULL
  for (i in 1:d) {
    x_alr = cbind(x_alr, alr(x_test[[i]], g = ref_g[[i]]))
  }
  
  z_hat = predict_fm(V, z_mean, x_alr)
  
  res = vector("list", d)
  for (i in 1:d) {
    temp = z_hat[,q_index(qs, i, "intrinsic")]
    res[[i]] = inv_alr(temp, ref_g[i])
  }
  
  return (res)
}

#' compute the KL divergence
#' @param p reference distribution
#' @param q candidate distribution
#' 
KL_div = function (p, q) {
  if (abs(sum(p) - 1) > 1e-10 || abs(sum(q) - 1) > 1e-10) {
    stop("Both p and q must sum to 1")
  }
  if (any(q == 0 & p > 0)) {
    return(Inf)  
  }
  sum(ifelse(p == 0, 0, p * log(p / q)))
}

comp_barplot <- function(X, years = NULL, col = NULL, legend = TRUE, ...) {
  X <- as.matrix(X)
  if (is.null(col)) {
    col <- rainbow(ncol(X))
  }
  
  bar_args <- list(
    height = t(X),
    col = col,
    beside = FALSE,
    legend.text = legend,
    args.legend = list(x = "topright", bty = "L")
  )
  
  if (!is.null(years)) {
    bar_args$names.arg <- years
  }
  
  bp <- do.call(barplot, c(bar_args, list(...)))
  invisible(bp)
}




### Below codes are defunct
#' computes ALR transformation for a compositional time series
#' 
#' @param x (n by q) matrix
#' @param g index for the reference group
#' 
# alr = function (x, g) {
#   n = nrow(x)
#   q = ncol(x)
#   
#   res = matrix(0, nrow = n, ncol = q - 1)
#   counter = 1
#   for (i in 1:q) {
#     if (i == g) {
#       next
#     } else {
#       res[,counter] = log(x[,i] / x[,g])
#       counter = counter + 1
#     }
#   }
#   
#   return (res)
# }

#' computes inverse alr transformation
#' 
#' @param z (n by (q - 1)) matrix
#' @param g index for the reference group
#' 
# inv_alr = function (z, g) {
#   z = exp(z)
#   q_minus = ncol(z)
#   if (g == 1) {
#     res = cbind(1, z)
#   } else if (g == q_minus + 1) {
#     res = cbind(z, 1)
#   } else {
#     res = cbind(z[,1:(g - 1)], 1, z[g:q_minus])
#   }
#   
#   res = res / rowSums(res)
#   
#   return (res)
# }

#' estimate the factor model after log ratio transformation
#' 
#' @param x a list of (n by q_j) arrays, j = 1,2,...,d
#' @param r number of factors
#' @param h lags used in factor estimation
#' @param reference_group a vector of indices (for each j) which are used as the 
#'                        reference group in log transformation
# alr_fm = function (x, r, h = 6, type = "alr", reference_group = NULL) {
#   d = length(x)
#   n = dim(x[[1]])[1]
#   qs = rep(NA, d)
#   for (i in 1:d) {
#     qs[i] = dim(x[[1]])[2]
#   }
#   if (!is.null(reference_group)) {
#     ref_g = reference_group
#   } else {
#     ref_g = qs # last group serves as reference in each composition
#   }
#   
#   x_alr = NULL
#   for (i in 1:d) {
#     x_alr = cbind(x_alr, alr(x[[i]], g = ref_g[i]))
#   }
#   
#   model = LYB_fm(x_alr, r = r, h = h)
#   
#   return (list("A" = model$V, "f_hat" = model$f_hat, "x" = x_alr,
#                "factor_model" = model, "r_hat" = model$r_hat,
#                "ref_g" = ref_g))
# } 






