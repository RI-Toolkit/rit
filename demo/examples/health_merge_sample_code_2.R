# set working directory and load packages
setwd("C:/Users/z5041702/GitHub/rit")

library(devtools)
load_all(export_all=FALSE)

param_file_5=US_HRS_5
param_file_3=US_HRS

model_type='F'

init_age=65
female=0
wave_index=8 # wave index
latent=0 # initial value of latent factor

trans_probs_3state=get_trans_probs(n_states=3, model_type, param_file=param_file_3, init_age, female, year = 2012, wave_index = 8, latent = 0)
trans_probs_5state=get_trans_probs(n_states=5, model_type, param_file=param_file_5, init_age, female, year = 2012, wave_index = 8, latent = 0)

l3=create_life_table(trans_probs_3state, init_age, init_state = 0, cohort = 100000)
l5=create_life_table(trans_probs_5state, init_age, init_state = 0, cohort = 100000)

l3_F=simulate_life_table(n_states=3, model_type, param_file=param_file_3, init_age, female, year = 2012, init_state = 0, wave_index = 8,latent=0,n_sim=100,cohort=100000,mean=TRUE)
l5_F=simulate_life_table(n_states=5, model_type, param_file=param_file_5, init_age, female, year = 2012, init_state = 0, wave_index = 8,latent=0,n_sim=100,cohort=100000,mean=TRUE)

simulated_path_3 <- simulate_health_state_paths(trans_probs_3state, init_age=65, init_state = 0, cohort = 10000)
simulated_path_5 <- simulate_health_state_paths(trans_probs_5state, init_age=65, init_state = 0, cohort = 10000)

# health5_stats(init_state=0, init_age=65, trans_probs=trans_probs_5state)

health_stats(n_states=3, init_age=65, init_state=0, trans_probs=trans_probs_3state)
health_stats(n_states=5, init_age=65, init_state=0, trans_probs=trans_probs_5state)

prob_plots(init_age=65, init_state = 0, trans_probs_3state)
prob_plots(init_age=65, init_state = 0, trans_probs_5state)



