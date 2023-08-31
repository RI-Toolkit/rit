#' function to get the first time leaving or entering different states for a number of individuals
#'
#' @param model_type
#' S for static model, T for trend model, F for frailty model
#'
#' @param state
#' 0 for first time leaving H state, only useful when initial state is 0
#' 1 for first time entering M state
#' 2 for first time entering D state
#' 3 for first time entering MD state
#' -1 for first time entering the dead state
#'
#' @param init_age
#' integer between 65 and 110 denoting initial age of individual. This needs to be same
#' initial age used in generation of `trans_probs` or `simulated_path`
#'
#' @param init_state
#' 0 for healthy, 1 for disabled
#'
#' @param trans_probs
#' a list of transition probability matrices, preferably generated from
#' \code{get_trans_probs}, only needed for static and trend models.
#'
#' @param simulated_path
#' the simulated path of individuals from the function \code{simulate_health_state_paths}
#'
#' @param female
#' 0 for male, 1 for female, compulsory variable for frailty model
#'
#' @param year
#' integer indicating current year, compulsory variable for frailty model
#'
#' @param wave_index
#' the wave index = (interview year - 1998)/2 + 1, compulsory variable for frailty model
#'
#' @param latent
#' initial value of latent factor, normally take the value 0, compulsory variable for frailty model
#'
#' @param param_file
#' string for file path of parameter file OR a tibble/dataframe of parameters, compulsory variable for frailty model
#'
#' @param n
#' integer denoting number of unique latent factor simulations
#'
#' @return
#' a column that consists the first time leaving or entering the state for a number of individuals
#'
#' @export
#'
#' @examples first_time_leave_H=health5_first_time_stats(health5_simulated_path_example, 0)
health5_first_time_stats=function(model_type, state, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female = NULL, year = NULL, wave_index = NULL, latent = NULL, param_file = NULL, n = 1000){

if (model_type=='S' | model_type=='T'){
    ##
    if (is.null(trans_probs) & is.null(simulated_path)) {
        stop('no transition probability matrices or simulated paths were provided')
    }
    if (!is.null(simulated_path) & !is.null(trans_probs)) {
        if (ncol(simulated_path) != 111-init_age +1) {
            stop('initial age does not correspond with size of simulated path')
        }
        # simulate path
        simulated_path <- simulated_path
    } else if (is.null(trans_probs)) {
        if (ncol(simulated_path) != 111-init_age +1) {
            stop('initial age does not correspond with size of simulated path')
        }

        simulated_path <- simulated_path
    } else {
        if (length(trans_probs) != 111-init_age) {
            stop('initial age does not correspond with number of transition probability matrices')
        }
        simulated_path <- simulate_health_state_paths(trans_probs, init_age, init_state, cohort = 10000)
    }
    ##
    first_time=matrix(nrow=nrow(simulated_path),ncol=1)
    if (state==0){
        for (i in 1:nrow(simulated_path)){
            first_time[i,]= which(simulated_path[i,] != 0)[1]-1.5
        }
        return(first_time)
    }
    if (state==1){
        for (i in 1:nrow(simulated_path)){
            first_time[i,]= which(simulated_path[i,] == 1)[1]-1.5
        }
    }
    if (state==2){
        for (i in 1:nrow(simulated_path)){
            first_time[i,]= which(simulated_path[i,] == 2)[1]-1.5
        }
    }
    if (state==3){
        for (i in 1:nrow(simulated_path)){
            first_time[i,]= which(simulated_path[i,] == 3)[1]-1.5
        }
    }
    if (state==-1){
        for (i in 1:nrow(simulated_path)){
            first_time[i,]= which(simulated_path[i,] == -1)[1]-1.5
        }
    }
}
if (model_type=='F'){
    first_time=matrix(nrow=n*10000,ncol=1)
    for (x in 1:n) {
        # simulate new frailty path for each iteration
        trans_probs <- get_trans_probs(n_states=5, model_type, param_file, init_age, female, year, wave_index, latent)
        simulated_path <- simulate_health_state_paths(trans_probs, init_age, init_state, cohort = 10000)
        for (i in 1:nrow(simulated_path)) {
            ##
            if (state==0){
                    first_time[(x-1)*10000+i,]= which(simulated_path[i,] != 0)[1]-1.5
            }
            if (state==1){
                    first_time[(x-1)*10000+i,]= which(simulated_path[i,] == 1)[1]-1.5
            }
            if (state==2){
                    first_time[(x-1)*10000+i,]= which(simulated_path[i,] == 2)[1]-1.5
            }
            if (state==3){
                    first_time[(x-1)*10000+i,]= which(simulated_path[i,] == 3)[1]-1.5
            }
            if (state==-1){
                    first_time[(x-1)*10000+i,]= which(simulated_path[i,] == -1)[1]-1.5
            }
            ##
        }
    }
}
return(first_time)
}

#' function to get the total time in different states for a number of individuals
#'
#' @param model_type
#' S for static model, T for trend model, F for frailty model
#'
#' @param state
#' 0 for total time in H state
#' 1 for total time in M state
#' 2 for total time in D state
#' 3 for total time in MD state
#' -1 for total time in dead state
#' 4 for total time alive or not in dead state
#'
#' @param init_age
#' integer between 65 and 110 denoting initial age of individual. This needs to be same
#' initial age used in generation of `trans_probs` or `simulated_path`
#'
#' @param init_state
#' 0 for healthy, 1 for disabled
#'
#' @param trans_probs
#' a list of transition probability matrices, preferably generated from
#' \code{get_trans_probs}, only needed for static and trend models.
#'
#' @param simulated_path
#' the simulated path of individuals from the function \code{simulate_health_state_paths}
#'
#' @param female
#' 0 for male, 1 for female, compulsory variable for frailty model
#'
#' @param year
#' integer indicating current year, compulsory variable for frailty model
#'
#' @param wave_index
#' the wave index = (interview year - 1998)/2 + 1, compulsory variable for frailty model
#'
#' @param latent
#' initial value of latent factor, normally take the value 0, compulsory variable for frailty model
#'
#' @param param_file
#' string for file path of parameter file OR a tibble/dataframe of parameters, compulsory variable for frailty model
#'
#' @param n
#' integer denoting number of unique latent factor simulations
#'
#' @return
#' a column that consists the total time in different states for a number of individuals
#'
#' @export
#'
#' @examples total_time_alive=health5_total_time_stats(health5_simulated_path_example, 4)
health5_total_time_stats=function(model_type, state, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female = NULL, year = NULL, wave_index = NULL, latent = NULL, param_file = NULL, n = 1000){
if (model_type=='S' | model_type=='T'){
    ##
    if (is.null(trans_probs) & is.null(simulated_path)) {
        stop('no transition probability matrices or simulated paths were provided')
    }
    if (!is.null(simulated_path) & !is.null(trans_probs)) {
        if (ncol(simulated_path) != 111-init_age +1) {
            stop('initial age does not correspond with size of simulated path')
        }
        # simulate path
        simulated_path <- simulated_path
    } else if (is.null(trans_probs)) {
        if (ncol(simulated_path) != 111-init_age +1) {
            stop('initial age does not correspond with size of simulated path')
        }
        simulated_path <- simulated_path
    } else {
        if (length(trans_probs) != 111-init_age) {
            stop('initial age does not correspond with number of transition probability matrices')
        }
        simulated_path <- simulate_health_state_paths(trans_probs, init_age, init_state, cohort = 10000)
    }
    ##
    total_time=matrix(nrow=nrow(simulated_path),ncol=1)
    if (state==0){
        for (i in 1:nrow(simulated_path)){
            if (simulated_path[i,1]==0){
                total_time[i,]= length(which(simulated_path[i,] == 0))-0.5
            }else{
                total_time[i,]= length(which(simulated_path[i,] == 0))
            }
        }
    }
    if(state==1){
        for (i in 1:nrow(simulated_path)){
            if (simulated_path[i,1]==1){
                total_time[i,]= length(which(simulated_path[i,] == 1))-0.5
            }else{
                total_time[i,]= length(which(simulated_path[i,] == 1))
            }
        }
    }
    if(state==2){
        for (i in 1:nrow(simulated_path)){
            if (simulated_path[i,1]==2){
                total_time[i,]= length(which(simulated_path[i,] == 2))-0.5
            }else{
                total_time[i,]= length(which(simulated_path[i,] == 2))
            }
        }
    }
    if(state==3){
        for (i in 1:nrow(simulated_path)){
            if (simulated_path[i,1]==3){
                total_time[i,]= length(which(simulated_path[i,] == 3))-0.5
            }else{
                total_time[i,]= length(which(simulated_path[i,] == 3))
            }
        }
    }
    if(state==-1){
        for (i in 1:nrow(simulated_path)){
            if (simulated_path[i,1]==-1){
                total_time[i,]= length(which(simulated_path[i,] == -1))-1
            }else{
                total_time[i,]= length(which(simulated_path[i,] == -1))-0.5
            }
        }
    }
    if(state==4){
        for (i in 1:nrow(simulated_path)){
            if (simulated_path[i,1]==-1){
                total_time[i,]= length(which(simulated_path[i,] != -1))
            }else{
                total_time[i,]= length(which(simulated_path[i,] != -1))-0.5
            }
        }
    }
}
    if (model_type=='F'){
        total_time=matrix(nrow=n*10000,ncol=1)
        for (x in 1:n) {
            # simulate new frailty path for each iteration
            trans_probs <- get_trans_probs(n_states=5, model_type, param_file, init_age, female, year, wave_index, latent)
            simulated_path <- simulate_health_state_paths(trans_probs, init_age, init_state, cohort = 10000)
            for (i in 1:nrow(simulated_path)) {
                ##
                if (state==0){
                        if (simulated_path[i,1]==0){
                            total_time[(x-1)*10000+i,]= length(which(simulated_path[i,] == 0))-0.5
                        }else{
                            total_time[(x-1)*10000+i,]= length(which(simulated_path[i,] == 0))
                        }
                }
                if(state==1){
                        if (simulated_path[i,1]==1){
                            total_time[(x-1)*10000+i,]= length(which(simulated_path[i,] == 1))-0.5
                        }else{
                            total_time[(x-1)*10000+i,]= length(which(simulated_path[i,] == 1))
                        }
                }
                if(state==2){
                        if (simulated_path[i,1]==2){
                            total_time[(x-1)*10000+i,]= length(which(simulated_path[i,] == 2))-0.5
                        }else{
                            total_time[(x-1)*10000+i,]= length(which(simulated_path[i,] == 2))
                        }
                }
                if(state==3){
                        if (simulated_path[i,1]==3){
                            total_time[(x-1)*10000+i,]= length(which(simulated_path[i,] == 3))-0.5
                        }else{
                            total_time[(x-1)*10000+i,]= length(which(simulated_path[i,] == 3))
                        }
                }
                if(state==-1){
                        if (simulated_path[i,1]==-1){
                            total_time[(x-1)*10000+i,]= length(which(simulated_path[i,] == -1))-1
                        }else{
                            total_time[(x-1)*10000+i,]= length(which(simulated_path[i,] == -1))-0.5
                        }
                }
                if(state==4){
                        if (simulated_path[i,1]==-1){
                            total_time[(x-1)*10000+i,]= length(which(simulated_path[i,] != -1))
                        }else{
                            total_time[(x-1)*10000+i,]= length(which(simulated_path[i,] != -1))-0.5
                        }
                }
                ##
            }
        }
    }
    return(total_time)
}


#' function to produce the mean and variance of a list of values
#'
#' @param input
#' the list of values to be studied
#' @return
#' mean and variance of the list of values
#'
#' @noRd
#'
#' @examples stats=health5_stats_produce(health5_first_time_stats(health5_simulated_path_example, 0))
health5_stats_produce=function(input){
    output=matrix(nrow=1, ncol = 2)
    colnames(output) <- c('expected_value', 'st_dev')
    output[1]=mean(input, na.rm = TRUE)
    output[2]=stats::sd(input, na.rm = TRUE)
    return(output)
}

#' Survival Statistics
#'
#' Produces statistics for 5-state model.
#'
#' @param model_type
#' S for static model, T for trend model, F for frailty model
#'
#' @param state
#' 0 for first time leaving H state, only useful when initial state is 0
#' 1 for first time entering M state
#' 2 for first time entering D state
#' 3 for first time entering MD state
#' -1 for first time entering the dead state
#'
#' @param init_age
#' integer between 65 and 110 denoting initial age of individual. This needs to be same
#' initial age used in generation of `trans_probs` or `simulated_path`
#'
#' @param init_state
#' 0 for healthy, 1 for disabled
#'
#' @param trans_probs
#' a list of transition probability matrices, preferably generated from
#' \code{get_trans_probs}, only needed for static and trend models.
#'
#' @param simulated_path
#' the simulated path of individuals from the function \code{simulate_health_state_paths}
#'
#' @param female
#' 0 for male, 1 for female, compulsory variable for frailty model
#'
#' @param year
#' integer indicating current year, compulsory variable for frailty model
#'
#' @param wave_index
#' the wave index = (interview year - 1998)/2 + 1, compulsory variable for frailty model
#'
#' @param latent
#' initial value of latent factor, normally take the value 0, compulsory variable for frailty model
#'
#' @param param_file
#' string for file path of parameter file OR a tibble/dataframe of parameters, compulsory variable for frailty model
#'
#' @param n
#' integer denoting number of unique latent factor simulations
#'
#' @return
#' dataframe output containing mean and standard deviation of different statistics
#'
#' @noRd
#'
#' @examples example
health5_stats <- function (model_type, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female = NULL, year = NULL, wave_index = NULL, latent = NULL, param_file = NULL, n = 1000) {
if (model_type=='S' | model_type=='T'){
        ##
    if (is.null(trans_probs) & is.null(simulated_path)) {
        stop('no transition probability matrices or simulated paths were provided')
    }
    if (!is.null(simulated_path) & !is.null(trans_probs)) {
        if (ncol(simulated_path) != 111-init_age +1) {
            stop('initial age does not correspond with size of simulated path')
        }
        # simulate path
        simulated_path <- simulated_path
    } else if (is.null(trans_probs)) {
        if (ncol(simulated_path) != 111-init_age +1) {
            stop('initial age does not correspond with size of simulated path')
        }
        simulated_path <- simulated_path
    } else {
        if (length(trans_probs) != 111-init_age) {
            stop('initial age does not correspond with number of transition probability matrices')
        }
        simulated_path <- simulate_health_state_paths(trans_probs, init_age, init_state, cohort = 10000)
    }
        ##

if (init_state == 0) {
    total_life=health5_total_time_stats(model_type, state=4, init_age, init_state, trans_probs, simulated_path)
    years_H=health5_total_time_stats(model_type, state=0, init_age, init_state, trans_probs, simulated_path)
    years_M=health5_total_time_stats(model_type, state=1, init_age, init_state, trans_probs, simulated_path)
    years_D=health5_total_time_stats(model_type, state=2, init_age, init_state, trans_probs, simulated_path)
    years_MD=health5_total_time_stats(model_type, state=3, init_age, init_state, trans_probs, simulated_path)
    years_disability=years_D+years_MD
    years_illness=years_M+years_MD
    first_H=health5_first_time_stats(model_type, state=0, init_age, init_state, trans_probs, simulated_path)
    first_M=health5_first_time_stats(model_type, state=1, init_age, init_state, trans_probs, simulated_path)
    first_D=health5_first_time_stats(model_type, state=2, init_age, init_state, trans_probs, simulated_path)
    first_MD=health5_first_time_stats(model_type, state=3, init_age, init_state, trans_probs, simulated_path)

    means <- c(mean(total_life), mean(years_H), mean(years_M),
               mean(years_D),mean(years_MD),mean(years_disability),mean(years_illness),mean(first_H, na.rm = TRUE),mean(first_M, na.rm = TRUE),mean(first_D, na.rm = TRUE),mean(first_MD, na.rm = TRUE))
    sds <- c(stats::sd(total_life), stats::sd(years_H), stats::sd(years_M),
             stats::sd(years_D),stats::sd(years_MD),stats::sd(years_disability),stats::sd(years_illness),
             stats::sd(first_H, na.rm = TRUE),stats::sd(first_M, na.rm = TRUE),stats::sd(first_D, na.rm = TRUE),stats::sd(first_MD, na.rm = TRUE))
    stats_df <- data.frame(
        'stats' = c('Mean years of life', 'Mean years in state H', 'Mean years in state M','Mean years in state D','Mean years in state MD',
                    'Mean years with disability','Mean years with illness','First time leaving state H','First time entering state M',
                    'First time entering state D', 'First time entering state MD'),
        'mean' = means,
        's.dev' = sds
    )
    return(stats_df)
}
##
if (init_state == 1) {
    total_life=health5_total_time_stats(model_type, state=4, init_age, init_state, trans_probs, simulated_path)
    years_H=health5_total_time_stats(model_type, state=0, init_age, init_state, trans_probs, simulated_path)
    years_M=health5_total_time_stats(model_type, state=1, init_age, init_state, trans_probs, simulated_path)
    years_D=health5_total_time_stats(model_type, state=2, init_age, init_state, trans_probs, simulated_path)
    years_MD=health5_total_time_stats(model_type, state=3, init_age, init_state, trans_probs, simulated_path)
    years_disability=years_D+years_MD
    years_illness=years_M+years_MD
    first_MD=health5_first_time_stats(model_type, state=3, init_age, init_state, trans_probs, simulated_path)

    means <- c(mean(total_life), mean(years_H), mean(years_M),
               mean(years_D),mean(years_MD),mean(years_disability),mean(years_illness),mean(first_MD, na.rm = TRUE))
    sds <- c(stats::sd(total_life), stats::sd(years_H), stats::sd(years_M),
             stats::sd(years_D),stats::sd(years_MD),stats::sd(years_disability),stats::sd(years_illness),
             stats::sd(first_MD, na.rm = TRUE))
    stats_df <- data.frame(
        'stats' = c('Mean years of life', 'Mean years in state H', 'Mean years in state M','Mean years in state D','Mean years in state MD',
                    'Mean years with disability','Mean years with illness','First time entering state MD'),
        'mean' = means,
        's.dev' = sds
    )
    return(stats_df)
}
##
if (init_state == 2) {
    total_life=health5_total_time_stats(model_type, state=4, init_age, init_state, trans_probs, simulated_path)
    years_H=health5_total_time_stats(model_type, state=0, init_age, init_state, trans_probs, simulated_path)
    years_M=health5_total_time_stats(model_type, state=1, init_age, init_state, trans_probs, simulated_path)
    years_D=health5_total_time_stats(model_type, state=2, init_age, init_state, trans_probs, simulated_path)
    years_MD=health5_total_time_stats(model_type, state=3, init_age, init_state, trans_probs, simulated_path)
    years_disability=years_D+years_MD
    years_illness=years_M+years_MD
    first_M=health5_first_time_stats(model_type, state=1, init_age, init_state, trans_probs, simulated_path)
    first_MD=health5_first_time_stats(model_type, state=3, init_age, init_state, trans_probs, simulated_path)

    means <- c(mean(total_life), mean(years_H), mean(years_M),
               mean(years_D),mean(years_MD),mean(years_disability),mean(years_illness),mean(first_M, na.rm = TRUE),mean(first_MD, na.rm = TRUE))
    sds <- c(stats::sd(total_life), stats::sd(years_H), stats::sd(years_M),
             stats::sd(years_D),stats::sd(years_MD),stats::sd(years_disability),stats::sd(years_illness),
             stats::sd(first_M, na.rm = TRUE),stats::sd(first_MD, na.rm = TRUE))
    stats_df <- data.frame(
        'stats' = c('Mean years of life', 'Mean years in state H', 'Mean years in state M','Mean years in state D','Mean years in state MD',
                    'Mean years with disability','Mean years with illness','First time entering state M',
                    'First time entering state MD'),
        'mean' = means,
        's.dev' = sds
    )
    return(stats_df)
}
##
if (init_state == 3) {
    total_life=health5_total_time_stats(model_type, state=4, init_age, init_state, trans_probs, simulated_path)
    years_H=health5_total_time_stats(model_type, state=0, init_age, init_state, trans_probs, simulated_path)
    years_M=health5_total_time_stats(model_type, state=1, init_age, init_state, trans_probs, simulated_path)
    years_D=health5_total_time_stats(model_type, state=2, init_age, init_state, trans_probs, simulated_path)
    years_MD=health5_total_time_stats(model_type, state=3, init_age, init_state, trans_probs, simulated_path)
    years_disability=years_D+years_MD
    years_illness=years_M+years_MD
    first_M=health5_first_time_stats(model_type, state=1, init_age, init_state, trans_probs, simulated_path)

    means <- c(mean(total_life), mean(years_H), mean(years_M),
               mean(years_D),mean(years_MD),mean(years_disability),mean(years_illness),mean(first_M, na.rm = TRUE))
    sds <- c(stats::sd(total_life), stats::sd(years_H), stats::sd(years_M),
             stats::sd(years_D),stats::sd(years_MD),stats::sd(years_disability),stats::sd(years_illness),
             stats::sd(first_M, na.rm = TRUE))
    stats_df <- data.frame(
        'stats' = c('Mean years of life', 'Mean years in state H', 'Mean years in state M','Mean years in state D','Mean years in state MD',
                    'Mean years with disability','Mean years with illness','First time entering state M'),
        'mean' = means,
        's.dev' = sds
    )
    return(stats_df)
}
}
if (model_type=='F'){
    if (init_state == 0) {
        total_life=health5_total_time_stats(model_type, state=4, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)
        years_H=health5_total_time_stats(model_type, state=0, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)
        years_M=health5_total_time_stats(model_type, state=1, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)
        years_D=health5_total_time_stats(model_type, state=2, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)
        years_MD=health5_total_time_stats(model_type, state=3, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)
        years_disability=years_D+years_MD
        years_illness=years_M+years_MD
        first_H=health5_first_time_stats(model_type, state=0, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)
        first_M=health5_first_time_stats(model_type, state=1, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)
        first_D=health5_first_time_stats(model_type, state=2, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)
        first_MD=health5_first_time_stats(model_type, state=3, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)

        means <- c(mean(total_life), mean(years_H), mean(years_M),
                   mean(years_D),mean(years_MD),mean(years_disability),mean(years_illness),mean(first_H, na.rm = TRUE),mean(first_M, na.rm = TRUE),mean(first_D, na.rm = TRUE),mean(first_MD, na.rm = TRUE))
        sds <- c(stats::sd(total_life), stats::sd(years_H), stats::sd(years_M),
                 stats::sd(years_D),stats::sd(years_MD),stats::sd(years_disability),stats::sd(years_illness),
                 stats::sd(first_H, na.rm = TRUE),stats::sd(first_M, na.rm = TRUE),stats::sd(first_D, na.rm = TRUE),stats::sd(first_MD, na.rm = TRUE))
        stats_df <- data.frame(
            'stats' = c('Mean years of life', 'Mean years in state H', 'Mean years in state M','Mean years in state D','Mean years in state MD',
                        'Mean years with disability','Mean years with illness','First time leaving state H','First time entering state M',
                        'First time entering state D', 'First time entering state MD'),
            'mean' = means,
            's.dev' = sds
        )
        return(stats_df)
    }
    ##
    if (init_state == 1) {
        total_life=health5_total_time_stats(model_type, state=4, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)
        years_H=health5_total_time_stats(model_type, state=0, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)
        years_M=health5_total_time_stats(model_type, state=1, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)
        years_D=health5_total_time_stats(model_type, state=2, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)
        years_MD=health5_total_time_stats(model_type, state=3, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)
        years_disability=years_D+years_MD
        years_illness=years_M+years_MD
        first_MD=health5_first_time_stats(model_type, state=3, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)

        means <- c(mean(total_life), mean(years_H), mean(years_M),
                   mean(years_D),mean(years_MD),mean(years_disability),mean(years_illness),mean(first_MD, na.rm = TRUE))
        sds <- c(stats::sd(total_life), stats::sd(years_H), stats::sd(years_M),
                 stats::sd(years_D),stats::sd(years_MD),stats::sd(years_disability),stats::sd(years_illness),
                 stats::sd(first_MD, na.rm = TRUE))
        stats_df <- data.frame(
            'stats' = c('Mean years of life', 'Mean years in state H', 'Mean years in state M','Mean years in state D','Mean years in state MD',
                        'Mean years with disability','Mean years with illness','First time entering state MD'),
            'mean' = means,
            's.dev' = sds
        )
        return(stats_df)
    }
    ##
    if (init_state == 2) {
        total_life=health5_total_time_stats(model_type, state=4, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)
        years_H=health5_total_time_stats(model_type, state=0, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)
        years_M=health5_total_time_stats(model_type, state=1, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)
        years_D=health5_total_time_stats(model_type, state=2, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)
        years_MD=health5_total_time_stats(model_type, state=3, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)
        years_disability=years_D+years_MD
        years_illness=years_M+years_MD
        first_M=health5_first_time_stats(model_type, state=1, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)
        first_MD=health5_first_time_stats(model_type, state=3, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)

        means <- c(mean(total_life), mean(years_H), mean(years_M),
                   mean(years_D),mean(years_MD),mean(years_disability),mean(years_illness),mean(first_M, na.rm = TRUE),mean(first_MD, na.rm = TRUE))
        sds <- c(stats::sd(total_life), stats::sd(years_H), stats::sd(years_M),
                 stats::sd(years_D),stats::sd(years_MD),stats::sd(years_disability),stats::sd(years_illness),
                 stats::sd(first_M, na.rm = TRUE),stats::sd(first_MD, na.rm = TRUE))
        stats_df <- data.frame(
            'stats' = c('Mean years of life', 'Mean years in state H', 'Mean years in state M','Mean years in state D','Mean years in state MD',
                        'Mean years with disability','Mean years with illness','First time entering state M',
                        'First time entering state MD'),
            'mean' = means,
            's.dev' = sds
        )
        return(stats_df)
    }
    ##
    if (init_state == 3) {
        total_life=health5_total_time_stats(model_type, state=4, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)
        years_H=health5_total_time_stats(model_type, state=0, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)
        years_M=health5_total_time_stats(model_type, state=1, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)
        years_D=health5_total_time_stats(model_type, state=2, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)
        years_MD=health5_total_time_stats(model_type, state=3, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)
        years_disability=years_D+years_MD
        years_illness=years_M+years_MD
        first_M=health5_first_time_stats(model_type, state=1, init_age, init_state, trans_probs = NULL, simulated_path = NULL, female, year, wave_index, latent, param_file, n)

        means <- c(mean(total_life), mean(years_H), mean(years_M),
                   mean(years_D),mean(years_MD),mean(years_disability),mean(years_illness),mean(first_M, na.rm = TRUE))
        sds <- c(stats::sd(total_life), stats::sd(years_H), stats::sd(years_M),
                 stats::sd(years_D),stats::sd(years_MD),stats::sd(years_disability),stats::sd(years_illness),
                 stats::sd(first_M, na.rm = TRUE))
        stats_df <- data.frame(
            'stats' = c('Mean years of life', 'Mean years in state H', 'Mean years in state M','Mean years in state D','Mean years in state MD',
                        'Mean years with disability','Mean years with illness','First time entering state M'),
            'mean' = means,
            's.dev' = sds
        )
        return(stats_df)
    }
}
}

