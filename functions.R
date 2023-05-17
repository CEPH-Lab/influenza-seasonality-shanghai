
####################################################################
## Code for all the functions used in seasonal-influenza-shanghai ##
####################################################################

## Function for estimating descriptive statistics in Table 1 ##

descriptive_table = function(x, y){
  table = tableby(x ~ y, data = master_df, 
                  numeric.stats = c("mean", "iqr", "q1q3"),
                  na.tableby(FALSE))
  
  return(summary(table))
}




## Negative binomial regression function ##

my.nbreg <- function(yname, d, out = NA){
  idx <- setdiff(1:nrow(d), out)
  d$y <- d[, yname]
  model <- glm.nb(y ~ normage + seasonal_trend + daily_temp_var + normhh_size + 
                    newocctype + normliveyear + adj_weekday + 
                    factor(irregular) + 
                    gender,
                  data = d[idx, ])
  return(model)
}


## Function for multinomial sampling

mymultinom <- function(n, q, x){
  u = runif(n, 0, 1)
  idx.v <- as.numeric(cut(u, q))
  res = unlist(lapply(c(1:n), function(i)
    runif(1, x[idx.v[i]], x[idx.v[i] + 1])))
  return(res)
}



## Function for putting regression coefficients into a matrix ##

my.data.matrix <- function(list.length, col.length, nsample, df){
  dm = list()
  for(k in 1:length(list.length)){
    dm[[k]] = matrix(NA, ncol = ncol(col.length), 
                     nrow = length(nsample))
    colnames(dm[[k]]) <- colnames(col.length)
  }
  for(k in 1:length(list.length)){
    for (j in 1:length(nsample)){
      dm[[k]][j, ] <- as.matrix(df[k, ])
    }
  }
  return(dm)
}




## Function for replacing the coefficients of interest ##

my.replace <- function(nsample1, var1, var2 = NA, nsample2 = NA, model){
  rep.list = list()
    for (j in 1:length(nsample1)){
      model$coef.tr[[j]] = list()
      TR = which(names(model$coefficients) == var1)
      if(!is.na(var2)){
        RE = which(names(model$coefficients) == var2)
      }
      model$coef.tr[[j]] = sapply(nsample1[j], function(x)
        replace(model$coefficients, TR, x))
      if(!is.na(nsample2[j])){
        model$coef.tr[[j]] = sapply(nsample2[j], function(y)
          replace(model$coef.tr[[j]], RE, y))
      }
      rownames(model$coef.tr[[j]]) <- names(model$coefficients)
    }
    rep.list <- model$coef.tr

  return(rep.list)
}




## Function for using the new coefficients to predict the outcome of interest ##

my.predict <- function(loop, model, columns, nsample, rep.df, newdata){

    y.fit.link = matrix(NA, ncol = length(columns),
                             nrow = length(nsample))
    colnames(y.fit.link) <- format(columns)
    ilink.out = matrix(NA, ncol = length(columns),
                            nrow = length(nsample))
    colnames(ilink.out) <- format(columns)
    for(k in 1:length(columns)){
      for (j in 1:length(nsample)){
        y.fit.link[j, k] = apply(rep.df[[j]], 2, function(y)
          apply(newdata[[k]][j, , drop = F], 1, function(x) as.numeric(y) %*% c(1, as.numeric(x))))
        ilink <- family(model)$linkinv
        ilink.out[j, k] <- ilink(y.fit.link[j, k])
      }
    }

  return(ilink.out)
}



## Function for estimating Means and Quantiles ##
my.mean_quantiles <- function(loop, prediction){

    mean.pred = NA
    q2.5 = NA
    q97.5 = NA
    for(k in 1:length(loop)){
      mean.pred[k] <- mean(prediction[, k], na.rm = T)
      q2.5[k] <- quantile(prediction[, k], 0.025, na.rm = T)
      q97.5[k] <- quantile(prediction[, k], 0.975, na.rm = T)
    }
    df <- cbind(mean.pred, q2.5, q97.5)
    colnames(df) <- c("mean", "2.5%", "97.5%")

  return(df)
}


## Synthesized Prediction Calculations using previous functions ##

PredictionProcess =  function(dataframe, dummydata, dates, sample1, sample2,
                              trend_var, daily_var, mod){
  dm = my.data.matrix(list.length = dates,
                      col.length = dummydata,
                      nsample = sample1,
                      df = dataframe)
  
  rep = my.replace(nsample1 = sample1, 
                   var1 = trend_var,
                   var2 = daily_var,  nsample2 = sample2, 
                   model = mod)
  
  pred = my.predict(model = mod, 
                    columns = as.Date(dates),
                    nsample = sample1, rep.df = rep, 
                    newdata = dm)
  
  mq <- my.mean_quantiles(loop = dates,
                          prediction = pred)
  
  return(mq)
}

PredictionProcessNomeans =  function(dataframe, dummydata, dates, sample1, sample2,
                                     ind_var, trend_var, daily_var, mod){
  dm = my.data.matrix(list.length = dates,
                      col.length = dummydata,
                      nsample = sample1,
                      df = dataframe)
  
  rep = my.replace(nsample1 = sample1, 
                   var1 = trend_var,
                   var2 = daily_var,  nsample2 = sample2, 
                   model = mod)
  
  pred = my.predict(model = mod, 
                    columns = as.Date(dates),
                    nsample = sample1, rep.df = rep, 
                    newdata = dm)
  
  
  return(pred)
}



## Conversion from daily to weekly 
daily_to_weekly_conversion = function(daily_estimates, n_weeks,
                                      days = 7, fn = sum){
  z = NULL
  for(i in 1:n_weeks){
    idx1 = max(1, (i*days - days))
    idx2 = min(i*days, length(daily_estimates))
    z[i] = fn(daily_estimates[idx1:idx2], na.rm = T)
  }
  
  return(z)
}




## Function for creating a boxplot 

myboxplot = function(data, Color = "lightblue", iteration, middle = NA){
  if(is.na(middle)){
    middle = data[3]
  }
  x = iteration
  width = 0.25
  
  rect(xleft = x - width, ybottom = 0, 
       xright = x + width, ytop = middle, border = T, col = Color)
  lines(c(x - width/3, x + width/3), c(data[1], data[1]))
  lines(c(x - width/3, x + width/3), c(data[5], data[5]))
  lines(c(x, x), c(data[1], data[5]))
}






