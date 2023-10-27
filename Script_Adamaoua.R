library(tidyverse)
library(rstan)
library(shinystan)
library(gridExtra)
library("bayesplot")
library("tidybayes")
library(matrixStats)
library("loo")

rstan_options (auto_write = TRUE)   # telling stan not to recompile code that has already been compiled
options (mc.cores = parallel::detectCores() )
set.seed(3)


path = "/Cameroon_Data_Fits/Regional_Fits/Adamaoua/"
setwd(path)
source("Fit_Summary.R")   # Loads the function that store the outputs of the samplying in a folder


TotalPop = 1309666

L_h = 3080 
beta_h = 0.022  
beta_v = 0.48 
mu_h = 1/(60.3*12)  
mu_v = 0.033*30   
sigma_h  = 0.1*30 
gamma_h = 0.0035*30 


# Organizing the cases data by region
data_path = "/Cameroon_Data_Fits/CameroonCasesData/RegionalData/"
Raw_Cases = read.csv(paste(data_path, "Region Adamaoua.csv", sep="") )
# less than 5 years
Raw_Cases_L5 <- Raw_Cases %>%
	filter(Age_group == "Age_L5")

Raw_Cases_G5 <- Raw_Cases %>%
	filter(Age_group == "Age_G5")

ReportCases <- Raw_Cases_L5$Cases + Raw_Cases_G5$Cases
Dates = as.Date(unique(Raw_Cases$Date) )
MonthlyCases <- ReportCases

# initial condition
X0 = rep(0, 7)  # initializing the initial population vector
X0[1] = TotalPop  # initial susceptible population
X0[2] = 1 
X0[3] = 1 
X0[4] = 0 
X0[5] = 5e3 
X0[6] = 0

initial_infected = sum(X0[2:3])

state = c(X0[2], X0[3])/sum(X0[2:3])
X0[2] = state[1]   # Eh
X0[3] = state[2]   # Ih

y0_vars = X0

# Preparing the parameters for stan
t0 = 0								# Initial time
n_months = length(ReportCases)     # number of weeks
ts = 1:n_months

data_seir = list(
	L_h = L_h, 
	beta_h = beta_h,
	beta_v = beta_v,
	mu_h = mu_h,
	sigma_h  = sigma_h,
	gamma_h = gamma_h,
	mu_v = mu_v,
	PIE = pi,
	TotalPop = TotalPop,
	
	p = p,
	t0 = t0,
	n_months = n_months,
	ts = ts,
	inference = 0,
	doprint = 0,
	
	# Data to fit
	Reported_cases = MonthlyCases,
	
	# Priors   
	p_i0 = c(log(500), 1),   
	p_epsilon = c(log(0.9), 0.2), #
	p_pp = c(0.5, 0.5), #
	p_b0 = c(log(100), 0.5),#
	p_Iv0 = c(log(5000), 0.5),#
	p_w = c(log(1/12), 0.1),  #  1/12
	p_wh = c(log(0.2), 0.5),  #
	p_delta = c(log(0.02), 0.5), #
	p_Lv = c(log(5000), 0.5), #
	p_tc = c(-2.5, 0.1), #
	p_phi = c(5),                                # parameters for the prior for phi 
	
	
	# initial condition
	y0_vars = y0_vars
	
)


MalariaModel_Full <- stan_model(paste(path, "Adamaoua.stan", sep="") )
fit_seir <- sampling(MalariaModel_Full, data = data_seir, iter = 5,  chains = 1)


# checking inference
pars = c( "i0", "epsilon", "p", "b0", "Iv0", 'w', 'delta', 'L_v', 'tc', 'w_h', 'phi') # specifying parameters of interest
print(fit_seir, pars = pars)
stan_dens(fit_seir, pars = pars, separate_chains = TRUE)
traceplot(fit_seir, pars = pars)
pairs(fit_seir, pars = pars)


# All age groups
smr_pred <- cbind(as.data.frame(summary(fit_seir, pars = "pred_cases", probs = c(0.05, 0.25, 0.5, 0.75, 0.95))$summary), ts, MonthlyCases,  days = as.Date(Dates))
colnames(smr_pred) <- make.names(colnames(smr_pred)) # to remove % in the col names

Reg1_CI <- ggplot(smr_pred, mapping = aes(x = days)) +
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill ="steelblue4",  alpha = 0.35) +
  geom_ribbon(aes(ymin = X25., ymax = X75.), fill ="steelblue4",  alpha = 0.45) +
  geom_line(mapping = aes(x = days, y = X50.), col = "navy", size = 2.5) +
  geom_point(mapping = aes(y = MonthlyCases, x=days), size = 5.5) +
  labs(x = "Time (Months)", y = "Reported cases") +  theme_bw() + theme(text = element_text(size = 45)) +
  ggtitle("Adamawa") + 
  theme(plot.title=element_text(hjust=0.08, vjust=-0.95, margin=margin(t=50,b=-30), size=45)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2, fill=NA), 
        axis.line = element_line(colour = "black")) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1.0))
Reg1_CI

# saving the projecttion plot in a .png file
png(filename= "Adamaoua_MalariaFit_F.png",  width = 700, height = 600)
plot(Reg1_CI)
dev.off()


