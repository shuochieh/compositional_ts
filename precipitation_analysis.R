temp_CA = station_data_get(list("CA" = 7))
x_CA = yearly_obs(rowSums(temp_CA$dta))
x_CA = x_CA / rowSums(x_CA)

# exploratory analysis
par(mfrow = c(1, 1))
comp_barplot(x_CA, years = 1895:2024)

x_CA_train = x_CA[1:94,]
x_CA_test = x_CA[95:130,]

model_RFM = rfm_sphere(list(sqrt(x_CA_train)), r = 10, h = 6)

par(mfrow = c(1, 2))
plot(model_RFM$f_hat[,1], type = "l", main = "RFM 1st factor", 
     xlab = "", ylab = "")
plot(model_RFM$f_hat[,2], type = "l", main = "RFM 2nd factor", 
     xlab = "", ylab = "")

model_ALR = alr_fm(list(x_CA_train), r = 10, h = 6)
plot(model_ALR$f_hat[,1], type = "l", main = "ALR 1st factor", 
     xlab = "", ylab = "")
plot(model_ALR$f_hat[,2], type = "l", main = "ALR 2nd factor", 
     xlab = "", ylab = "")

# In-sample fit
mean_Euc = colMeans(x_CA_train)

RFM_preds = array(0, dim = c(94, 12, 11))
RFM_loss = rep(0, 11)
for (r in 1:11) {
  model_RFM = rfm_sphere(list(sqrt(x_CA_train)), r = r, h = 6)
  RFM_preds[,,r] = (predict_rfm(list(sqrt(x_CA_train)), model_RFM)[[1]])^2
  
  RFM_loss[r] = norm(x_CA_train - RFM_preds[,,r], "F")^2
}

ALR_preds = array(0, dim = c(94, 12, 11))
ALR_loss = rep(0, 11)
for (r in 1:11) {
  model_ALR = alr_fm(list(x_CA_train), r = r, h = 6)
  ALR_preds[,,r] = predict_alr(list(x_CA_train), model_ALR)[[1]]
  
  ALR_loss[r] = norm(x_CA_train - ALR_preds[,,r], "F")^2
}

RFM_loss_KL = rep(0, 11)
ALR_loss_KL = rep(0, 11)
mean_loss_KL = 0
for (r in 1:11) {
  for (i in 1:94) {
    RFM_loss_KL[r] = RFM_loss_KL[r] + KL_div(RFM_preds[i,,r], x_CA_train[i,])
    ALR_loss_KL[r] = ALR_loss_KL[r] + KL_div(ALR_preds[i,,r], x_CA_train[i,])
    if (r == 1) {
      mean_loss_KL = mean_loss_KL + KL_div(mean_Euc, x_CA_train[i,])
    }
  }
}

par(mfrow = c(1, 2))
plot(sqrt(RFM_loss / 94), type = "l", col = "steelblue", xlab = "number of factors",
     ylab = "RMSE", main = "In-sample prediction errors",
     pch = 19, ylim = c(0, 0.2))
points(sqrt(RFM_loss / 94), col = "steelblue", pch = 19)
lines(sqrt(ALR_loss / 94), col = "firebrick", lty = 2)
points(sqrt(ALR_loss / 94), col = "firebrick", pch = 19)
legend("topright", col = c("steelblue", "firebrick"),
       lty = c(1, 2), legend = c("Sph", "ALR"))

plot(RFM_loss_KL / 94, type = "l", col = "steelblue", xlab = "number of factors",
     ylab = "KL divergence", main = "In-sample prediction errors",
     pch = 19, ylim = c(0, 0.2))
points(RFM_loss_KL / 94, col = "steelblue", pch = 19)
lines(ALR_loss_KL / 94, col = "firebrick", lty = 2)
points(ALR_loss_KL / 94, col = "firebrick", pch = 19)


# Prediction

# Factor interpretation
