source("util.R")

temp_CA = station_data_get(list("CA" = 7))
x_CA = yearly_obs(rowSums(temp_CA$dta))
x_CA = x_CA / rowSums(x_CA)

# exploratory analysis
par(mfrow = c(1, 1))
comp_barplot(x_CA, years = 1895:2024)

x_CA_train = x_CA[1:94,]
x_CA_test = x_CA[95:130,]

# train factor models
model_RFM = rfm_sphere(list(sqrt(x_CA_train)), r = 10, h = 6)

model_ALR = logR_fm(list(x_CA_train), r = 10, h = 6, type = "alr")
model_ILR = logR_fm(list(x_CA_train), r = 10, h = 6, type = "ilr")

model_sph_sq = lfm_sphere(list(sqrt(x_CA_train)), r = 10, h = 6)
model_sph_raw = lfm_sphere(list(x_CA_train), r = 10, h = 6)

# par(mfrow = c(1, 2))
# plot(model_RFM$f_hat[,1], type = "l", main = "RFM 1st factor", 
#      xlab = "", ylab = "")
# plot(model_RFM$f_hat[,2], type = "l", main = "RFM 2nd factor", 
#      xlab = "", ylab = "")


# In-sample fit
n_test = 94

RFM_preds = array(0, dim = c(n_test, 12, 10))
RFM_loss = rep(0, 10)
RFM_KL = rep(0, 10)
for (r in 1:10) {
  model_RFM = rfm_sphere(list(sqrt(x_CA_train)), r = r, h = 6)
  RFM_preds[,,r] = (predict_rfm(list(sqrt(x_CA_train)), model_RFM)[[1]])^2
  RFM_loss[r] = sqrt((norm(x_CA_train - RFM_preds[,,r], "F")^2) / n_test)
  
  for (i in 1:n_test) {
    RFM_KL[r] = RFM_KL[r] + KL_div(x_CA_train[i,], RFM_preds[i,,r])
  }
}

ALR_preds = array(0, dim = c(n_test, 12, 10))
ALR_loss = rep(0, 10)
ALR_KL = rep(0, 10)
for (r in 1:10) {
  model_ALR = logR_fm(list(x_CA_train), r = r, h = 6, type = "alr")
  ALR_preds[,,r] = predict_logR(list(x_CA_train), model_ALR)$pred[[1]]
  ALR_loss[r] = sqrt((norm(x_CA_train - ALR_preds[,,r], "F")^2) / n_test)
  
  for (i in 1:n_test) {
    ALR_KL[r] = ALR_KL[r] + KL_div(x_CA_train[i,], ALR_preds[i,,r])
  }
}

ILR_preds = array(0, dim = c(n_test, 12, 10))
ILR_loss = rep(0, 10)
ILR_KL = rep(0, 10)
for (r in 1:10) {
  model_ILR = logR_fm(list(x_CA_train), r = r, h = 6, type = "ilr")
  ILR_preds[,,r] = predict_logR(list(x_CA_train), model_ILR)$pred[[1]]
  ILR_loss[r] = sqrt((norm(x_CA_train - ILR_preds[,,r], "F")^2) / n_test)
  
  for (i in 1:n_test) {
    ILR_KL[r] = ILR_KL[r] + KL_div(x_CA_train[i,], ILR_preds[i,,r])
  }
}

Lin_sph_preds = array(0, dim = c(n_test, 12, 10))
Lin_sph_loss = rep(0, 10)
Lin_sph_KL = rep(0, 10)
for (r in 1:10) {
  model_linsph = lfm_sphere(list(sqrt(x_CA_train)), r = r, h = 6)
  Lin_sph_preds[,,r] = (predict_lfm(list(sqrt(x_CA_train)), model_linsph)[[1]])^2
  Lin_sph_loss[r] = sqrt((norm(x_CA_train - Lin_sph_preds[,,r], "F")^2) / n_test)
  
  for (i in 1:n_test) {
    Lin_sph_KL[r] = Lin_sph_KL[r] + KL_div(x_CA_train[i,], 
                                           Lin_sph_preds[i,,r] / sum(Lin_sph_preds[i,,r]))
  }
}

Lin_raw_preds = array(0, dim = c(n_test, 12, 10))
Lin_raw_loss = rep(0, 10)
Lin_raw_KL = rep(0, 10)
for (r in 1:10) {
  model_linraw = lfm_sphere(list(x_CA_train), r = r, h = 6)
  Lin_raw_preds[,,r] = (predict_lfm(list(x_CA_train), model_linraw)[[1]])
  Lin_raw_loss[r] = sqrt((norm(x_CA_train - Lin_raw_preds[,,r], "F")^2) / n_test)
  
  cat("r = ", r, "\n")
  for (i in 1:n_test) {
    cat(" i = ", i,  round(KL_div(x_CA_train[i,], Lin_raw_preds[i,,r] / sum(Lin_raw_preds[i,,r])), 3), 
        "\n")
    Lin_raw_KL[r] = Lin_raw_KL[r] + KL_div(x_CA_train[i,], 
                                           Lin_raw_preds[i,,r] / sum(Lin_raw_preds[i,,r]))
  }
}


par(mfrow = c(1, 2))

ylim = range(c(RFM_loss, ALR_loss, ILR_loss, Lin_sph_loss, Lin_raw_loss))
ylim[1] = 0
plot(RFM_loss, col = "steelblue", type = "l", xlab = "number of factors",
     ylab = "RMSE", main = "In-sample prediction errors", ylim = ylim, lwd = 2)
points(RFM_loss, col = "steelblue", pch = 19)

lines(ALR_loss, col = "red", lty = 2)
points(ALR_loss, col = "red", pch = 18)
lines(ILR_loss, col = "red4", lty = 2)
points(ILR_loss, col = "red4", pch = 17)

lines(Lin_sph_loss, col = "green", lty = 3)
points(Lin_sph_loss, col = "green", pch = 15)
lines(Lin_raw_loss, col = "green4", lty = 3)
points(Lin_raw_loss, col = "green4", pch = 14)

ylim = range(c(RFM_KL, ALR_KL, ILR_KL, Lin_sph_KL))
ylim[1] = 0
plot(RFM_KL, col = "steelblue", type = "l", xlab = "number of factors",
     ylab = "KL divergence", main = "In-sample prediction errors", ylim = ylim, lwd = 2)
points(RFM_KL, col = "steelblue", pch = 19)

lines(ALR_KL, col = "red", lty = 2)
points(ALR_KL, col = "red", pch = 18)
lines(ILR_KL, col = "red4", lty = 2)
points(ILR_KL, col = "red4", pch = 17)

lines(Lin_sph_KL, col = "green", lty = 3)
points(Lin_sph_KL, col = "green", pch = 15)
# lines(Lin_raw_loss, col = "green4", lty = 3)
# points(Lin_raw_loss, col = "green4", pch = 14)
legend("topright", 
       lty = c(1, 2, 2, 3, 3), pch = c(19, 18, 17, 15, 15),
       legend = c("RFM", "ALR", "ILR", "Lin-Sph", "Lin-Raw"),
       col = c("steelblue", "red", "red4", "green", "green4"))


# Pseudo-prediction
n_test = 36

RFM_preds = array(0, dim = c(n_test, 12, 10))
RFM_loss = rep(0, 10)
RFM_KL = rep(0, 10)
for (r in 1:10) {
  model_RFM = rfm_sphere(list(sqrt(x_CA_train)), r = r, h = 6)
  
  RFM_preds[,,r] = (predict_rfm(list(sqrt(x_CA_test)), model_RFM)[[1]])^2
  RFM_loss[r] = sqrt((norm(x_CA_test - RFM_preds[,,r], "F")^2) / n_test)
  for (i in 1:n_test) {
    RFM_KL[r] = RFM_KL[r] + KL_div(x_CA_test[i,], RFM_preds[i,,r])
  }
}

ALR_preds = array(0, dim = c(n_test, 12, 10))
ALR_loss = rep(0, 10)
ALR_KL = rep(0, 10)
for (r in 1:10) {
  model_ALR = logR_fm(list(x_CA_train), r = r, h = 6, type = "alr")
  
  ALR_preds[,,r] = predict_logR(list(x_CA_test), model_ALR)$pred[[1]]
  ALR_loss[r] = sqrt((norm(x_CA_test - ALR_preds[,,r], "F")^2) / n_test)
  for (i in 1:n_test) {
    ALR_KL[r] = ALR_KL[r] + KL_div(x_CA_test[i,], ALR_preds[i,,r])
  }
}

ILR_preds = array(0, dim = c(n_test, 12, 10))
ILR_loss = rep(0, 10)
ILR_KL = rep(0, 10)
for (r in 1:10) {
  model_ILR = logR_fm(list(x_CA_train), r = r, h = 6, type = "ilr")
  
  ILR_preds[,,r] = predict_logR(list(x_CA_test), model_ILR)$pred[[1]]
  ILR_loss[r] = sqrt((norm(x_CA_test - ILR_preds[,,r], "F")^2) / n_test)
  for (i in 1:n_test) {
    ILR_KL[r] = ILR_KL[r] + KL_div(x_CA_test[i,], ILR_preds[i,,r])
  }
}

Lin_sph_preds = array(0, dim = c(n_test, 12, 10))
Lin_sph_loss = rep(0, 10)
Lin_sph_KL = rep(0, 10)
for (r in 1:10) {
  model_linsph = lfm_sphere(list(sqrt(x_CA_train)), r = r, h = 6)
  
  Lin_sph_preds[,,r] = (predict_lfm(list(sqrt(x_CA_test)), model_linsph)[[1]])^2
  Lin_sph_loss[r] = sqrt((norm(x_CA_test - Lin_sph_preds[,,r], "F")^2) / n_test)
  for (i in 1:n_test) {
    Lin_sph_KL[r] = Lin_sph_KL[r] + KL_div(x_CA_test[i,], 
                                           Lin_sph_preds[i,,r] / sum(Lin_sph_preds[i,,r]))
  }
}

Lin_raw_preds = array(0, dim = c(n_test, 12, 10))
Lin_raw_loss = rep(0, 10)
Lin_raw_KL = rep(0, 10)
for (r in 1:10) {
  model_linraw = lfm_sphere(list(x_CA_train), r = r, h = 6)
  
  Lin_raw_preds[,,r] = (predict_lfm(list(x_CA_test), model_linraw)[[1]])
  Lin_raw_loss[r] = sqrt((norm(x_CA_test - Lin_raw_preds[,,r], "F")^2) / n_test)
  for (i in 1:n_test) {
    Lin_raw_KL[r] = Lin_raw_KL[r] + KL_div(x_CA_test[i,], 
                                           Lin_raw_preds[i,,r] / sum(Lin_raw_preds[i,,r]))
  }
}

par(mfrow = c(1, 2))

ylim = range(c(RFM_loss, ALR_loss, ILR_loss, Lin_sph_loss, Lin_raw_loss))
ylim[1] = 0
plot(RFM_loss, col = "steelblue", type = "l", xlab = "number of factors",
     ylab = "RMSE", main = "Pseudo-prediction errors", ylim = ylim, lwd = 2)
points(RFM_loss, col = "steelblue", pch = 19)

lines(ALR_loss, col = "red", lty = 2)
points(ALR_loss, col = "red", pch = 18)
lines(ILR_loss, col = "red4", lty = 2)
points(ILR_loss, col = "red4", pch = 17)

lines(Lin_sph_loss, col = "green", lty = 3)
points(Lin_sph_loss, col = "green", pch = 15)
lines(Lin_raw_loss, col = "green4", lty = 3)
points(Lin_raw_loss, col = "green4", pch = 14)

ylim = range(c(RFM_KL, ALR_KL, ILR_KL, Lin_sph_KL))
ylim[1] = 0
plot(RFM_KL, col = "steelblue", type = "l", xlab = "number of factors",
     ylab = "KL divergence", main = "Pseudo-prediction errors", ylim = ylim, lwd = 2)
points(RFM_KL, col = "steelblue", pch = 19)

lines(ALR_KL, col = "red", lty = 2)
points(ALR_KL, col = "red", pch = 18)
lines(ILR_KL, col = "red4", lty = 2)
points(ILR_KL, col = "red4", pch = 17)

lines(Lin_sph_KL, col = "green", lty = 3)
points(Lin_sph_KL, col = "green", pch = 15)
# lines(Lin_raw_loss, col = "green4", lty = 3)
# points(Lin_raw_loss, col = "green4", pch = 14)
legend("topright", 
       lty = c(1, 2, 2, 3, 3), pch = c(19, 18, 17, 15, 14),
       legend = c("RFM", "ALR", "ILR", "Lin-Sph", "Lin-Raw"),
       col = c("steelblue", "red", "red4", "green", "green4"))






# Out-of-sample prediction
# Pseudo-prediction
n_test = 36

RFM_preds = array(0, dim = c(n_test, 12, 10))
RFM_loss = rep(0, 10)
RFM_KL = rep(0, 10)
for (r in 1:10) {
  for (i in 1:n_test) {
    if (i == 1) {
      train_data = x_CA_train
    } else {
      train_data = rbind(x_CA_train, x_CA_test[(i - 1),])
    }
  }
  model_RFM = rfm_sphere(list(sqrt(train_data)), r = r, h = 6)
  fac = model_RFM$f_hat
  alpha_hat = simple_AR(fac)
  
  RFM_preds[,,r] = (predict_rfm(list(sqrt(x_CA_test)), model_RFM)[[1]])^2
  RFM_loss[r] = sqrt((norm(x_CA_test - RFM_preds[,,r], "F")^2) / n_test)
  for (i in 1:n_test) {
    RFM_KL[r] = RFM_KL[r] + KL_div(x_CA_test[i,], RFM_preds[i,,r])
  }
}

# Factor interpretation





















