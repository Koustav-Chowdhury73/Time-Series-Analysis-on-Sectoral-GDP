library(aTSA)
library(strucchange)

data = read.csv("data.csv")
l = 10
n = 50

for (k in 2:25) {
  y1 = (data[, k] / data$GDP.at.factor.cost) * 100
  y = ts(y1[1:n])
  
  m1 = lm(y ~ seq(n))
  y_hat = m1$fitted.values
  x = y - y_hat
  
  m2 = estimate(
    x,
    p = 3,
    d = 1,
    q = 3,
    PDQ = c(0, 1, 0),
    S = n %/% 2
  )
  
  f = forecast(m2, lead = l)
  
  trend_forecast = cbind(rep(1, (n + l)), seq(n + l)) %*% matrix(coefficients(m1))
  err_pred = x - m2[["residuals"]]
  
  t1 = min(y1[1:(n + l)])
  t2 = max(y1[1:(n + l)])
  temp = c(1.25 * t1 - 0.25 * t2, 1.25 * t2 - 0.25 * t1)
  plot(
    y1[1:(n + l)] ~ seq(1950 + 1, 1950 + n + l),
    type = "l",
    col = "blue",
    ylim = temp,
    ylab = colnames(data)[k],
    xlab = "year"
  )
  points((trend_forecast + c(err_pred, f[, 2])) ~ seq(1950 + 1, 1950 + n + l),
         type = "l",
         col = "red")
  points((trend_forecast + c(err_pred, f[, 4])) ~ seq(1950 + 1, 1950 + n + l), type = "l")
  points((trend_forecast + c(err_pred, f[, 5])) ~ seq(1950 + 1, 1950 + n + l), type = "l")
  
  d = 0
  temp = 10000
  for (i in 1:(n - 1)) {
    l1 = lm(y[1:i] ~ seq(i))
    l2 = lm(y[(i + 1):n] ~ seq(n - i))
    sse = sum((y[1:i] - l1$fitted.values)^2 + (y[(i + 1):n] - l2$fitted.values)^2)
    if (sse < temp) {
      temp = sse
      d = i
    }
  }
  
  w = min(d - 1, n - d)
  chow = sctest(y[(d - w):(d + w)] ~ seq(2 * w + 1), type = "Chow", point = w)
  cat(colnames(data)[k], w, 1950 + d, chow$p.value, "\n")
  
  l1 = lm(y[1:d] ~ seq(d))
  l2 = lm(y[(d + 1):n] ~ seq(n - d))
  y_hat = c(l1$fitted.values, l2$fitted.values)
  x = y - y_hat
  
  m = estimate(
    x,
    p = 3,
    d = 1,
    q = 3,
    PDQ = c(0, 1, 0),
    S = n %/% 2
  )
  
  f = forecast(m, lead = l)
  
  trend_forecast = c(l1$fitted.values,
                     cbind(rep(1, (n - d + l)), seq(n - d + l)) %*% matrix(coefficients(l2)))
  err_pred = x - m[["residuals"]]
  
  t1 = min(y1[1:(n + l)])
  t2 = max(y1[1:(n + l)])
  temp = c(1.25 * t1 - 0.25 * t2, 1.25 * t2 - 0.25 * t1)
  plot(
    y1[1:(n + l)] ~ seq(1950 + 1, 1950 + n + l),
    type = "l",
    col = "blue",
    ylim = temp,
    ylab = colnames(data)[k],
    xlab = "year"
  )
  points((trend_forecast + c(err_pred, f[, 2])) ~ seq(1950 + 1, 1950 + n + l),
         type = "l",
         col = "red")
  points((trend_forecast + c(err_pred, f[, 4])) ~ seq(1950 + 1, 1950 + n + l), type = "l")
  points((trend_forecast + c(err_pred, f[, 5])) ~ seq(1950 + 1, 1950 + n + l), type = "l")
  
  d1 = 0
  d2 = 0
  temp = 10000
  for (i in 1:(n - 2)) {
    for (j in (i + 1):(n - 1)) {
      l1 = lm(y[1:i] ~ seq(i))
      l2 = lm(y[(i + 1):j] ~ seq(j - i))
      l3 = lm(y[(j + 1):n] ~ seq(n - j))
      sse = sum((y[1:i] - l1$fitted.values)^2 + (y[(i + 1):j] - l2$fitted.values)^2 + (y[(j + 1):n] - l3$fitted.values)^2)
      if (sse < temp) {
        temp = sse
        d1 = i
        d2 = j
      }
    }
  }
  
  w = min(d1 - 1, (d2 - d1 - 1) %/% 2, n - d2)
  if (w > 1) {
    chow1 = sctest(y[(d1 - w):(d1 + w)] ~ seq(2 * w + 1), type = "Chow", point = w)
    chow2 = sctest(y[(d2 - w):(d2 + w)] ~ seq(2 * w + 1), type = "Chow", point = w)
    cat(colnames(data)[k],
        w,
        1950 + d1,
        1950 + d2,
        chow1$p.value,
        chow2$p.value,
        "\n")
  } else {
    cat(colnames(data)[k], w, 1950 + d1, 1950 + d2, "\n")
  }
  
  l1 = lm(y[1:d1] ~ seq(d1))
  l2 = lm(y[(d1 + 1):d2] ~ seq(d2 - d1))
  l3 = lm(y[(d2 + 1):n] ~ seq(n - d2))
  
  y_hat = c(l1$fitted.values, l2$fitted.values, l3$fitted.values)
  x = y - y_hat
  
  m = estimate(
    x,
    p = 3,
    d = 1,
    q = 3,
    PDQ = c(0, 1, 0),
    S = n %/% 2
  )
  
  f = forecast(m, lead = l)
  
  trend_forecast = c(l1$fitted.values,
                     l2$fitted.values,
                     cbind(rep(1, (n - d2 + l)), seq(n - d2 + l)) %*% matrix(coefficients(l3)))
  err_pred = x - m[["residuals"]]
  
  t1 = min(y1[1:(n + l)])
  t2 = max(y1[1:(n + l)])
  temp = c(1.25 * t1 - 0.25 * t2, 1.25 * t2 - 0.25 * t1)
  plot(
    y1[1:(n + l)] ~ seq(1950 + 1, 1950 + n + l),
    type = "l",
    col = "blue",
    ylim = temp,
    ylab = colnames(data)[k],
    xlab = "year"
  )
  points((trend_forecast + c(err_pred, f[, 2])) ~ seq(1950 + 1, 1950 + n + l),
         type = "l",
         col = "red")
  points((trend_forecast + c(err_pred, f[, 4])) ~ seq(1950 + 1, 1950 + n + l), type = "l")
  points((trend_forecast + c(err_pred, f[, 5])) ~ seq(1950 + 1, 1950 + n + l), type = "l")
}