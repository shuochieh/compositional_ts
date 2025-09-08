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

temp_CA = station_data_get(list("CA" = 7))
x_CA = yearly_obs(rowSums(temp_CA$dta))