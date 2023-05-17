rm(list = ls())
gc()

#########################################################################################
## Code for calculating the descriptive statistics and running the regression analysis ##
#########################################################################################

## Libraries needed ##
INSTALL = F

if(INSTALL == T){
  install.packages("arsenal")
  install.packages("table1")
  install.packages("ggplot2")
  install.packages("MASS")
  install.packages("prettyGraphs")
}

library(arsenal)
library(table1)
library(ggplot2)
library(MASS)
library(prettyGraphs)


## Source file for functions ##
source("functions.R")


################################################
## Read in files for the statistical analysis ##
################################################

cnt_data = read.csv("statistical_analysis/files/clean_data.csv") # Processed data from contact diaries

daily_temp = read.csv("statistical_analysis/files/Meteorological_data_2017_2018.csv") # Daily meteorological data in Shanghai
names(daily_temp)[names(daily_temp) == "Date"] <- "diary_date"



#############################################################
## Estimate seasonal trend and daily temperature variation ##
#############################################################

## Seasonal Trend calculations
spline_max_temp <- smooth.spline(daily_temp$maxtemp, nknots = 4)

summary(spline_max_temp)

daily_temp$seasonal_trend <- spline_max_temp$y

## Daily temperature variation calculations
daily_temp$daily_temp_var <- NA

for(i in 1:length(daily_temp$diary_date)){
  daily_temp$daily_temp_var[i] <- daily_temp$maxtemp[which(daily_temp$diary_date ==
                                                               daily_temp$diary_date[i])] - daily_temp$seasonal_trend[i]
}

## Merge with contact data
master_df <- merge(cnt_data, daily_temp[, c("diary_date","seasonal_trend", "daily_temp_var")], 
                   by = "diary_date")


############################
## Descriptive Statistics ##
############################

## Estimates for Table 1
descriptive_table(master_df$gender, master_df$cnt_n)

descriptive_table(master_df$agegroup, master_df$cnt_n)

descriptive_table(master_df$newocctype, master_df$cnt_n)
summary(master_df$cnt_n[which(is.na(master_df$newocctype))])

descriptive_table(master_df$normliveyear, master_df$cnt_n)
summary(master_df$cnt_n[which(is.na(master_df$normliveyear))])

descriptive_table(master_df$adj_weekday, master_df$cnt_n)

descriptive_table(master_df$irregular, master_df$cnt_n)
summary(master_df$cnt_n[which(is.na(master_df$irregular))])

total_contacts_stats = table1(~ cnt_n, data = master_df)

covariates_n = table1(~ gender + agegroup + 
                        newocctype + adj_weekday +
                        factor(irregular) + normliveyear, 
                      data = master_df)

## Figure 1A p
mycol = c("#c86d3f", "#eaaf47", "black")

ggplot(daily_temp, aes(x = diary_date, y = maxtemp)) +
  geom_line(aes(x = diary_date, y = seasonal_trend, group = 1), 
            color = mycol[1],
            show.legend = T) +
  geom_segment(aes(xend = diary_date, yend = seasonal_trend, 
                   colour = "Daily variation"), 
               color = mycol[2],
               show.legend = T) + 
  geom_point(aes(x = diary_date, y = maxtemp, fill = "Maximum temperature"),
             fill = mycol[3],
             color = mycol[3], show.legend = T) +
  ylab("Temperature (°C)") +
  xlab("Date") +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        axis.line = element_line(),
        legend.position = "bottom") +
  scale_y_continuous(expand = c(0,0), limits = c(-1,40)) +
  scale_colour_manual("", values = c("Maximum temperature" =  mycol[3],
                                     "Seasonal trend" = mycol[1], 
                                     "Daily variaton" = mycol[2])) +
  guides(color = guide_legend(override.aes = 
                                list(shape = c(16, NA, NA),
                                     linetype = c(0,1,1),
                                     color = c(mycol[3], mycol[1], mycol[2]))))




#########################
## Regression Analysis ##
#########################
## Put variables used in the regression into new data frame
y <- c("cnt_n")
x <- c("normage", "seasonal_trend",  "daily_temp_var", "normhh_size", 
       "newocctype", "normliveyear", "adj_weekday", "irregular", "gender") 

regression_df <- master_df[, c(y, x)]

## Run the full model
full.model <- my.nbreg(yname = names(regression_df)[1], d = regression_df)
summary(full.model)
  
## Check for outliers
cd <- cooks.distance(full.model)
id <- which(as.numeric(cd) > mean(cd)*20)
outlier_idx <- id

## Run model without outliers
final.model <- my.nbreg(yname = names(regression_df)[1], 
                        d = regression_df,  
                        out = outlier_idx)
summary(final.model)


## Save into CSV
#Save Regression into Excel file
pvalue = coef(summary(final.model))[,4]

asterisk = NULL
asterisk[pvalue< 0.1] <- "."
asterisk[pvalue < 0.05] <- "*"  
asterisk[pvalue < 0.01] <- "**"
asterisk[pvalue < 0.001] <- "***"
asterisk[is.na(pvalue)] <- NA

variability =  confint(final.model)

regression_output = data.frame(final.model$coefficients, pvalue, asterisk, 
                               variability[, 1], variability[, 2])
colnames(regression_output) = c("Coefficients", "p-value", "significance",
                                "2.5%", "97.5%")




#####################################
## Predict  the number of contacts ##
#####################################

## Read in regression output 
regression_estimates = read.csv("statistical_analysis/files/regression_output.csv")


## Multinomial sampling for seasonal trend and daily temperature variation
b1 = confint(final.model, level = 0.95)[3,] 
mystep = 0.01
my.levels <- c(0.00001, seq(mystep, .99999, by = mystep), 0.9999)

ci.list <- lapply(my.levels, function(x)
  confint(final.model, level = x)[c(3,4),])
head(ci.list)

my.names <- as.numeric(gsub( " %", "", c(colnames(ci.list[[1]]),
                                         colnames(ci.list[[length(ci.list)]]))))

Names.v <- c(seq(my.names[3], my.names[1], by = 100*(mystep/2)),
             seq(my.names[2], my.names[4], by = 100*(mystep/2)))


## Seasonal Trend
tmp_trend1 <-  rev(unlist(lapply(ci.list, function(x)
  x[1,1])))
tmp_trend2 <-  unlist(lapply(ci.list, function(x)
  x[1,2]))

tmptrend <- c(tmp_trend1, tmp_trend2)
names(tmptrend) <- c(Names.v)

## Daily temperature variation
tmp_var1 <-  rev(unlist(lapply(ci.list, function(x)
  x[2,1])))
tmp_var1 <-  unlist(lapply(ci.list, function(x)
  x[2,2]))

tmpvar <- c(tmp_var1, tmp_var1)
names(tmpvar) <- c(Names.v)

## MULTINOMIAL SAMPLE
nexperiment <- 1000 #Number of samples

trend_multi <- mymultinom(nexperiment, q = my.levels, x = tmptrend)
var_multi <- mymultinom(nexperiment, q = my.levels, x = tmpvar)

## These will vary based on the sample multinomials

mulnom.trend <- read.csv("statistical_analysis/files/multinomial_trends.csv")[, -1] 
mulnom.var <- read.csv("statistical_analysis/files/multinomial_variations.csv")[, -1]


## Data frame for prediction coefficients
## Rename variables to match regression output
names(master_df)[which(names(master_df) == "newocctype_Not.Employed")] = "newocctype_Not Employed"
names(master_df)[which(names(master_df) == "normliveyear_..6.years")] = "normliveyear_< 6 years"
names(master_df)[which(names(master_df) == "normliveyear_6.10.years")] = "normliveyear_6-10 years"


covariates <- c("normage", "seasonal_trend", "daily_temp_var", "normhh_size", 
              "newocctype_Not Employed", "newocctype_Students", 
              "normliveyear_< 6 years", "normliveyear_6-10 years", 
              "adj_weekday_SAT", "adj_weekday_SUN", "irregular_1", 
              "irregular_2", "gender_Female") 

prediction_coefs <- master_df[, covariates] 
head(prediction_coefs)

values <- c(mean(master_df$normage, na.rm = T), NA, NA, 
            mean(master_df$normhh_size, na.rm = T), 
            rep(0, times = 9)) 


prediction_reg_values <- prediction_coefs[1, , drop = F]

for (i in 1:length(daily_temp$diary_date)) {
  values[2] <- daily_temp$seasonal_trend[which(as.Date(daily_temp$diary_date,
                                           na.rm = T) == daily_temp$diary_date[i])][1]
  values[3] <- daily_temp$daily_temp_var[which(as.Date(daily_temp$diary_date,
                                           na.rm = T) == daily_temp$diary_date[i])][1]
  prediction_reg_values[i, ] <- c(values)
}

prediction_reg_values$normage <- as.numeric(prediction_reg_values$normage)
prediction_reg_values$seasonal_trend <- as.numeric(prediction_reg_values$seasonal_trend)
prediction_reg_values$daily_temp_var <- as.numeric(prediction_reg_values$daily_temp_var)
prediction_reg_values$normhh_size <- as.numeric(prediction_reg_values$normhh_size)


## Adjust for weekend days
daily_temp$day <- weekdays(as.Date(daily_temp$diary_date))

wd <- c("Sunday", "Monday", "Tuesday", "Wednesday",
        "Thursday", "Friday", "Saturday")

id.w = list()
for(w in 1:length(wd)){
  id.w[[w]] <- which(daily_temp$day == wd[w])
}
for(w in 1:length(wd)) {
  for(l in 1:length(id.w [[1]])){
    for(k in 1:length(id.w [[7]])){
      Su <- id.w [[1]][l]
      Sa <- id.w [[7]][k]
      prediction_reg_values[Su,][, "adj_weekday_SUN"] <- 1
      prediction_reg_values[Sa,][, "adj_weekday_SAT"] <- 1
    }
  }
}

## Include the extended holiday
extendedHoliday = which(daily_temp$diary_date >= "2018-01-22" &
                          daily_temp$diary_date <= "2018-02-22")

for(i in 1:length(extendedHoliday)){
  prediction_reg_values[extendedHoliday[i], "irregular_2"] <- 1
}

## Predict the number of contacts
predicted_cnt = PredictionProcessNomeans(dataframe = prediction_reg_values, dummydata = prediction_coefs, 
                                         dates = daily_temp$diary_date, sample1 = mulnom.trend, 
                                         sample2 = mulnom.var, trend_var = "seasonal_trend", 
                                         daily_var = "daily_temp_var", mod = final.model)


## Convert to Weekly number of contacts
weekly_cnt = matrix(NA, nrow = nrow(predicted_cnt), ncol = 53)

for(i in 1:nrow(predicted_cnt)){
  for(j in 1:ncol(weekly_cnt)){
    if(j == 1){
      weekly_cnt[i, j] = predicted_cnt[i, j]
    } 
    
    weekly_cnt[i, 2:ncol(weekly_cnt)] = daily_to_weekly_conversion(predicted_cnt[i, -1], 
                                                                   n_weeks = round(365/7), 
                                                                   fn = mean)
  }
}

mean_weekly_cnt = apply(weekly_cnt, MARGIN = 2, function(x)
  mean(x, na.rm = T))
ci_weekly_cnt = apply(weekly_cnt, MARGIN = 2, function(x)
  quantile(x, probs = c(0.025, 0.975), na.rm = T))


## Figure 1B
wk1 <- c(1, seq(2, length(daily_temp$diary_date), by = 7))
wk2 <- c(seq(1, length(daily_temp$diary_date), by = 7))

x_at = 1:length(daily_temp$diary_date[wk1][-1])

plot(mean_weekly_cnt[-1], type = "n", axes = F,
     xlab = "", ylab = "Average Daily Number of Contacts", ylim = c(5, 25), cex = 0.5)
box(bty = "l")
axis(1, at = c(x_at), 
     label = daily_temp$diary_date[wk1][-1],  cex.axis = 0.8, las = 2)
axis(2)

polygon(c(rev(x_at), x_at), c(rev(ci_weekly_cnt[2, ][-1]), ci_weekly_cnt[1, ][-1]), 
        col = add.alpha("turquoise3", alpha = 0.5), border = NA)
lines(mean_weekly_cnt[-1], col = "turquoise3")





