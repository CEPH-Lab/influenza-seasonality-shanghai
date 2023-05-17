rm(list = ls())
gc()

####################################
## Code for the modeling analysis ##
####################################

## Libraries needed ##
INSTALL = F

if(INSTALL == T){
  install.packages("aweek")
  install.packages("prettyGraphs")
}

library(aweek)
library(prettyGraphs)


## Source file for functions ##
source("functions.R")

##############################
## Running MCMC Instructions##
##############################

# To run MCMC, change the mcmc parameter (params$MCMC) equal to 1 (line 141)
# and set simulations equal to 0 (line 155).

######################################
## Plotting Posterior distributions ##
######################################

## These file take a long time to read into the program 

parameter_distributions = read.csv("modeling_analysis/files/calibrated_param_distributions.csv", header = T) ## Posterior distributions from MCMC

n_infections = read.csv("modeling_analysis/files/estimated_ili_pos.csv", header = T)

net_reproduction_number = read.csv("modeling_analysis/files/net_reproduction_number.csv", header = T)

ili_pos_points <- read.csv("modeling_analysis/files/postive_ili.csv")


## Remove first half of the simulations
idx_remove = c(1:(nrow(parameter_distributions) / 2))

population = 24860000 ## Population of Shanghai
mycol = c("#2A91A9", "#2D747D")


## Incidence of ILI
n_infections = n_infections[, -c(1,2, ncol(n_infections))]

estimated_ili = (n_infections[-idx_remove, ]/population)*10000

mean_ili = apply(estimated_ili, MARGIN = 2, function(x)
  mean(x, na.rm = T))

var_ili = apply(estimated_ili, MARGIN = 2, function(x) 
  quantile(x, probs = c(0.025, 0.975), na.rm = T))

n_ili = ili_pos_points$est_Inc

ili_inc = (n_ili / population)*10000


## Plot only the weeks of the influenza season
idx_estimates = c(1:18)
idx_data = c(49:(49+17))

dates = seq(as.Date("2017-12-04"), as.Date("2018-09-30"), by = "1 week")

## Figure 2A
plot(dates[idx_estimates], var_ili[2, idx_estimates], 
     type = "n", axes = F, xlab = "", ylab = "Incidence of ILI+")
axis(1, at = dates[idx_estimates], labels = dates[idx_estimates], las = 2)
axis(2)
polygon(c(dates[idx_estimates], rev(dates[idx_estimates])),
        c(var_ili[1, idx_estimates], rev(var_ili[2, idx_estimates])),
        border = F, col = add.alpha(mycol[1], alpha = 0.5))
lines(dates[idx_estimates], mean_ili[idx_estimates], lwd = 3, col = mycol[1])
points(dates[idx_estimates], ili_inc[idx_data], col = mycol[2], pch = 19)


## Net reproduction number
net_Rt = net_reproduction_number[-idx_remove, -c(1, 2,  ncol(net_reproduction_number))]

net_Rt_weekly = matrix(NA, ncol = round(ncol(net_Rt)/7), nrow = nrow(net_Rt))

## This takes a long time to run
for(i in 1:nrow(net_Rt)){
  net_Rt_weekly[i, ] = daily_to_weekly_conversion(as.numeric(net_Rt[i, ]), 
                                                  n_weeks = round(ncol(net_Rt)/7),
                                                  fn = mean)
}

mean_netRt = apply(net_Rt_weekly, MARGIN = 2, function(x)
  mean(x, na.rm = T))

quantiles_netRt = apply(net_Rt_weekly, MARGIN = 2, function(x)
  quantile(x, probs = c(0.025, 0.975), na.rm = T))

# Figure 2B
plot(dates[idx_estimates], quantiles_netRt[2, idx_estimates], 
     type = "n", axes = F, ylab = "Net Reproduction Number",
     xlab = "", ylim = c(0.6, 1.4))
axis(1, at = dates[idx_estimates], labels = dates[idx_estimates], las = 2)
axis(2)
polygon(c(dates[idx_estimates], rev(dates[idx_estimates])),
        c(quantiles_netRt[1, idx_estimates], rev(quantiles_netRt[2, idx_estimates])),
        col = add.alpha("salmon3", alpha = 0.5), 
        border = F)
lines(dates[idx_estimates], mean_netRt[idx_estimates], col = "salmon3")



## Final infection attack rate
quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975)

AR_quantiles = quantile(parameter_distributions$AR[-idx_remove], probs = quantiles)*100

## Figure 2C
boxplot(AR_quantiles, col = mycol[1])




##############################################
## Running the Influenza Transmission Model ##
##############################################

model.dir = "modeling_analysis/seasonal_influenza_model/"

RUN_SIM = F ## To run the transmission model, set RUN_SIM = T

params = read.csv(paste0(model.dir, "inputs/parameters.csv"), header = T)

## MCMC parameter - 0 = no mcmc runs; 1 = run mcmc
params$MCMC = 0

location.names = c("Shanghai")
locations = c("Shanghai_estimated")
populations = c(24860000)

beta_adjustment = c(1)
scales = c(1) 

t0 = format(c(seq(as.Date("2017-11-01"), as.Date("2018-03-31"), by = "1 month")), "%Y-%m-%d")
t0names = c("Nov.01", "Dec.01", "Jan.01", "Feb.01", "Mar.01")


## Number of simulations that will run 
simulations = 1000


## Use the posterior distributions for the parameters
sample_mcmc = parameter_distributions[-idx_remove, ]

sample_n = sample(which(sample_mcmc$Accept == 1), simulations)

used_dist = sample_mcmc[sample_n, ]


## Vector of the number of contacts predicted in the statistical analysis
### Adjusted for different starting days

orig = "2017-10-01"
dates = seq(as.Date(orig), as.Date("2019-09-30"), by = "1 day")
wks = as.integer(format(as.Date(as.aweek(dates, week_start = 1)), "%W"))
date_df = data.frame(dates, wks)

cntMatrix = list()
contacts = list()
# Loop over locations
for(i in 1:length(locations)){ 
  contacts[[i]] = list()
  cntMatrix[[i]] = read.csv(paste0("modeling_analysis/files/", location.names[i], "_meancontacts.csv"))[, -1]
  
  # Loop over starting date
  for(j in 1:length(t0)){
    contacts[[i]][[j]] = list()
    idx = NULL
    idx = which(dates == t0[j])
    contacts[[i]][[j]] = unlist(c(cntMatrix[[i]][idx:length(cntMatrix[[i]])], 
                                  cntMatrix[[i]][1:(idx - 1)]))
  }
}

# Adjust the parameters as needed
parameters = list()

# Loop over locations and beta adjustments
for(k in 1:length(beta_adjustment)){
  parameters[[k]] = list()
  for(i in 1:length(locations)){
    parameters[[k]][[i]] = list()
    for(q in 1:simulations){
      parameters[[k]][[i]][[q]] = params
      parameters[[k]][[i]][[q]][which(names(params) == "Population")] = populations[i]
      parameters[[k]][[i]][[q]][which(names(params) == "beta")] = (used_dist$Beta_0[q]*scales[i]*beta_adjustment[i])
      parameters[[k]][[i]][[q]][which(names(params) == "I0")] = (used_dist$I0_0[q])
      parameters[[k]][[i]][[q]][which(names(params) == "Dispersion")] = (used_dist$Size_0[q])
      parameters[[k]][[i]][[q]][which(names(params) == "Report_Rate")] = (used_dist$RR_0[q])
      parameters[[k]][[i]][[q]][which(names(params) == "holiday_Report_rate")] = (used_dist$RRholiday_0[q])
      parameters[[k]][[i]][[q]][which(names(params) == "start_holiday")] = -1
      parameters[[k]][[i]][[q]][which(names(params) == "end_holiday")] = -1
      parameters[[k]][[i]][[q]][which(names(params) == "first_week")] = 49
    }
  }
}

if(RUN_SIM == T){
  setwd(model.dir)
  
  main = paste(model.dir, "Differential_main", sep = "")
  
  for(i in 1:length(locations)){
    for(j in 1:length(t0)){
      for(q in 1:simulations){
        
        ## files for running system
        tmpcnt_file = paste0(model.dir, "inputs/tmp/tmpcontacts.csv")
        tmpparam_file = paste0(model.dir, "inputs/tmp/TmpParams.csv")
        tmpoutput = paste0(model.dir, "outputs/", locations[i], "/", t0[j], "/Sim_", q, "/") 
        
        if(dir.exists(paste0(model.dir, "inputs/tmp/")) == F){
          dir.create(file.path(paste0(model.dir, "inputs/tmp/")), recursive = T)
        }
        
        if(dir.exists(tmpoutput) == F){
          dir.create(file.path(tmpoutput), recursive = T)
        }

        write.csv(contacts[[i]][[j]], tmpcnt_file)
        write.csv(parameters[[k]][[i]][[q]], tmpparam_file)
        
        
        system(paste0('"', main, '" "', tmpparam,'" "', tmpcnt, '" "', tmpoutput, '"'))
        
        print(paste0("simulation ", q, "is done"))
      }
      
      print(paste0("t = ", j, "is done"))
    }
    
    print(paste0(locations[i], "is done"))
  }
}



##################################################
## Estimating the Potential Reproduction Number ##
##################################################

params

parameter_distributions

cntMatrix

beta = parameter_distributions[(nrow(parameter_distributions)/2):nrow(parameter_distributions), ]

# summary(beta$RRholiday_0)

y = 1/params$Tg

cnt_gamma = (cntMatrix[[1]]/y)

beta_adjustment = c(0.9, 1, 1.1)

dates = seq(as.Date("2017/10/01", origin = "2017/01/01"), as.Date("2018/09/30", origin = "2018/01/01"),
            by = "1 week")[-1]

reproduction = list()
weekly = list()
for(k in 1:length(beta_adjustment)){
  reproduction[[k]] = matrix(NA, nrow = length(cnt_gamma), ncol = length(beta$Beta_0))
  weekly[[k]] = matrix(NA, nrow = round(length(cnt_gamma)/7), ncol = length(beta$Beta_0))
  for(i in 1:length(beta$Beta_0)){
    reproduction[[k]][, i] = cnt_gamma*(beta$Beta_0[i])*beta_adjustment[k]
    weekly[[k]][, i] = daily_to_weekly_conversion(reproduction[[k]][, i][-1], 
                                                  n_weeks = round(length(cnt_gamma)/7),
                                                  fn = mean)
  }
}

dec_mean_Rt = apply(weekly[[1]], MARGIN = 1, FUN = mean)
dec_lower_Rt = apply(weekly[[1]], MARGIN = 1, function(x)
  quantile(x, 0.025))
dec_upper_Rt = apply(weekly[[1]], MARGIN = 1, function(x)
  quantile(x, 0.975))

inc_mean_Rt = apply(weekly[[3]], MARGIN = 1, FUN = mean)
inc_lower_Rt = apply(weekly[[3]], MARGIN = 1, function(x)
  quantile(x, 0.025))
inc_upper_Rt = apply(weekly[[3]], MARGIN = 1, function(x)
  quantile(x, 0.975))

mean_Rt = apply(weekly[[2]], MARGIN = 1, FUN = mean)
lower_Rt = apply(weekly[[2]], MARGIN = 1, function(x)
  quantile(x, 0.025))
upper_Rt = apply(weekly[[2]], MARGIN = 1, function(x)
  quantile(x, 0.975))


mycol = c("#89C97F", "#2A91A9", "#F2C55A")

## Figure 3A
plot(dates, inc_upper_Rt, type = "n", ylim = c(0.5, 1.6))
polygon(c(dates, rev(dates)),
        c(lower_Rt, rev(upper_Rt)),
        col = add.alpha(mycol[2], alpha = 0.5),
        border = F)

polygon(c(dates, rev(dates)),
        c(inc_lower_Rt, rev(inc_upper_Rt)),
        col = add.alpha(mycol[3], alpha = 0.5),
        border = F)

polygon(c(dates, rev(dates)),
        c(dec_lower_Rt, rev(dec_upper_Rt)),
        col = add.alpha(mycol[1], alpha = 0.5),
        border = F)

lines(dates, mean_Rt, col = mycol[2])
lines(dates, dec_mean_Rt, col = mycol[1])
lines(dates, inc_mean_Rt, col = mycol[3])

abline(h = 1, lty = 2, col = "gray")


#########################################################
## Estimating the Transmission Dynamics from the Model ##
#########################################################

locations_adj = c("Shanghai_estimated", "Shanghai_10dec", "Shanghai_10inc")
populations = c(24860000, 24860000, 24860000)

t0names_month = format(as.Date(t0), "%b")

percent_scale = 100
scale_factor = 10000
quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975)
first = 1

simulations = 1000

## Read in Data and assign to variables 
cum_ili = list()
weekly_ili = list()

infection_attack_rate = list()
AR_quantiles = list()

weekly_inc = list()
mean_weekly_inc = list()
quantiles_weekly_inc = list()

peak_week_incidence = list()
peak_week = list()
PWInc_quantiles = list()
PW_quantiles = list()

eff_Rt = list()
weekly_Rt = list()
mean_weekly_Rt = list()
quantiles_weekly_Rt = list()

for(l in 1:length(locations_adj)){
  cum_ili[[l]] = list()
  weekly_ili[[l]] = list()
  
  infection_attack_rate[[l]] = list()
  AR_quantiles[[l]] = list()
  
  weekly_inc[[l]] = list()
  mean_weekly_inc[[l]] = list()
  quantiles_weekly_inc[[l]] = list()
  
  peak_week_incidence[[l]] = list()
  peak_week[[l]] = list()
  PWInc_quantiles[[l]] = list()
  PW_quantiles[[l]] = list()
  
  eff_Rt[[l]] = list()
  weekly_Rt[[l]] = list()
  mean_weekly_Rt[[l]] = list()
  quantiles_weekly_Rt[[l]] = list()
  
  for(t in 1:length(t0)){
    
    #### Read in Files ####
    cum_ili[[l]][[t]] <- read.csv(paste0("modeling_analysis/files/", locations_adj[l], "/CSV/", 
                                         t0names_month[t], "/cumulative_ili.csv"))[, -first]
    
    cum_ili[[l]][[t]] = cum_ili[[l]][[t]][, -ncol(cum_ili[[l]][[t]])]
    
    weekly_ili[[l]][[t]] <- read.csv(paste0("modeling_analysis/files/", locations_adj[l], "/CSV/", 
                                            t0names_month[t], "/weekly_ili.csv"))[, -first]
    
    weekly_ili[[l]][[t]] = weekly_ili[[l]][[t]][, -ncol(weekly_ili[[l]][[t]])]
    
    eff_Rt[[l]][[t]] = read.csv(paste0("modeling_analysis/files/", locations_adj[l], "/CSV/", 
                                       t0names_month[t], 
                                       "/effective_reproduction_number.csv"))[, -first]
    
    eff_Rt[[l]][[t]] = eff_Rt[[l]][[t]][, -ncol(eff_Rt[[l]][[t]])]
    
    
    #### Calculate Final Infection Attack Rate ####
    infection_attack_rate[[l]][[t]] = (cum_ili[[l]][[t]][, ncol(cum_ili[[l]][[t]])] 
                                       / populations[l])*percent_scale
    
    AR_quantiles[[l]][[t]] = quantile(infection_attack_rate[[l]][[t]],
                                      probs = quantiles)
    
    #### Calculate Weekly Incidence ####
    weekly_inc[[l]][[t]] = (weekly_ili[[l]][[t]] / populations[l])*scale_factor
    
    mean_weekly_inc[[l]][[t]] = apply(weekly_inc[[l]][[t]], 
                                      MARGIN = 2, FUN = mean)
    
    quantiles_weekly_inc[[l]][[t]] = apply(weekly_inc[[l]][[t]],
                                           MARGIN = 2, function(x) 
                                             quantile(x, c(0.025, 0.975)))
    
    #### Identify Peak Week and Peak Week Incidence ####
    peak_week_incidence[[l]][[t]] = NA
    peak_week[[l]][[t]] = NA
    
    for(s in 1:simulations){
      peak_week_incidence[[l]][[t]][s] = max(weekly_inc[[l]][[t]][s, ])
      
      peak_week[[l]][[t]][s] = which(weekly_inc[[l]][[t]][s, ] == 
                                       max(weekly_inc[[l]][[t]][s, ]))
    }
    
    peak_week[[l]][[t]] = as.numeric(format(as.Date(t0[t]) +
                                              peak_week[[l]][[t]]*7,
                                            "%W"))
    
    PWInc_quantiles[[l]][[t]] = quantile(peak_week_incidence[[l]][[t]],
                                         probs = quantiles)
    
    PW_quantiles[[l]][[t]] = quantile(peak_week[[l]][[t]], probs = quantiles)
    
    #### Effective Rt ####
    weekly_Rt[[l]][[t]] = matrix(NA, nrow = nrow(eff_Rt[[l]][[t]]), 
                                 ncol = round(ncol(eff_Rt[[l]][[t]])/7))
    
    for(i in 1:nrow(eff_Rt[[l]][[t]])){
      weekly_Rt[[l]][[t]][i, ] = daily_to_weekly_conversion(as.numeric(eff_Rt[[l]][[t]][i, ]),
                                                            n_weeks = ncol(weekly_Rt[[l]][[t]]),
                                                            fn = mean)
    }
    
    mean_weekly_Rt[[l]][[t]] = apply(weekly_Rt[[l]][[t]],
                                     MARGIN = 2, FUN = mean)
    
    quantiles_weekly_Rt[[l]][[t]] = apply(weekly_Rt[[l]][[t]],
                                          MARGIN = 2, function(x)
                                            quantile(x, c(0.025, 0.975)))
    
  }
}


placement = seq(1, (length(t0)*3), by = 3)

mycol = c("#2A91A9", "#89C97F", "#F2C55A")

idx_shanghai = c(1:3)

idx_all = c(1, 4, 5)

dates = seq(as.Date("2018-01-01"), as.Date("2018-05-01"), by = "1 month")
dates_at = format(dates, "%W")
dates_labels = format(dates, "%b")


## Figure 3B 

## Infection Attack Rates ##
plot(c(0, length(t0)*3), c(0, 60), type = "n", 
     ylab = "Infection Attack Rate (%)", xlab = "", axes = F)
axis(2)
for(t in 1:length(t0)){
  myboxplot(AR_quantiles[[1]][[t]], iteration = placement[t], Color = mycol[1])
  myboxplot(AR_quantiles[[2]][[t]], iteration = placement[t] + 1, Color = mycol[2])
  myboxplot(AR_quantiles[[3]][[t]], iteration = placement[t] + 2, Color = mycol[3])
  
  axis(1, at = placement[t] + 1, labels = t0names[t], las = 2)
}
legend("topright", legend = locations[idx_shanghai], fill = mycol)

## Peak Week Incidence ##
plot(c(0, length(t0)*3), c(0, 6), type = "n", 
     ylab = "Peak Week Incidence (per 10,000)", xlab = "", axes = F)
axis(2)
for(t in 1:length(t0)){
  myboxplot(PWInc_quantiles[[1]][[t]], iteration = placement[t], Color = mycol[1])
  myboxplot(PWInc_quantiles[[2]][[t]], iteration = placement[t] + 1, Color = mycol[2])
  myboxplot(PWInc_quantiles[[3]][[t]], iteration = placement[t] + 2, Color = mycol[3])
  
  axis(1, at = placement[t] + 1, labels = t0names[t], las = 2)
}
legend("topright", legend = locations[idx_shanghai], fill = mycol)

## Peak Week ##
plot(c(0, length(t0)*3), c(0, 16), type = "n", 
     ylab = "Peak Week", xlab = "", axes = F)
axis(2, at = dates_at, labels = dates_labels, las = 2)
for(t in 1:length(t0)){
  myboxplot(PW_quantiles[[1]][[t]], iteration = placement[t], Color = mycol[1])
  myboxplot(PW_quantiles[[2]][[t]], iteration = placement[t] + 1, Color = mycol[2])
  myboxplot(PW_quantiles[[3]][[t]], iteration = placement[t] + 2, Color = mycol[3])
  
  axis(1, at = placement[t] + 1, labels = t0names[t], las = 2)
}
legend("topleft", legend = locations[idx_shanghai], fill = mycol)







