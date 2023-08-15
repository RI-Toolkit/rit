#' function to get the first time leaving or entering different states for a number of individuals
#' @param simulated_path
#' the simulated path of individuals from the function \code{simulate_health_state_paths}
#' @param state
#' 0 for first time leaving H state, only useful when initial state is 0
#' 1 for first time entering M state
#' 2 for first time entering D state
#' 3 for first time entering MD state
#' -1 for first time entering the dead state
#'
#' @return
#' a column that consists the first time leaving or entering the state for a number of individuals
#' @export health5_first_time_stats
#'
#' @examples first_time_leave_H=health5_first_time_stats(health5_simulated_path_example, 0)
health5_first_time_stats=function(simulated_path, state){
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
    return(first_time)
}

#' function to get the total time in different states for a number of individuals
#' @param simulated_path
#' the simulated path of individuals from the function \code{simulate_health_state_paths}
#' @param state
#' 0 for total time in H state
#' 1 for total time in M state
#' 2 for total time in D state
#' 3 for total time in MD state
#' -1 for total time in dead state
#' 4 for total time alive or not in dead state
#'
#' @return
#' a column that consists the total time in different states for a number of individuals
#' @export health5_total_time_stats
#'
#' @examples total_time_alive=health5_total_time_stats(health5_simulated_path_example, 4)
health5_total_time_stats=function(simulated_path, state){
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
        total_time=matrix(nrow=nrow(simulated_path),ncol=1)
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
    return(total_time)
}


#' function to produce the mean and variance of a list of values
#'
#' @param input
#' the list of values to be studied
#' @return
#' mean and variance of the list of values
#' @export health5_stats_produce
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
#' @param init_age
#' integer between 65 and 110 denoting initial age of individual
#' @param init_state
#' 0 for H state, 1 for M state, 2 for D state, 3 for MD state
#' @param trans_probs
#' list of transition probability matrices, generated from \code{health5_get_trans_probs}.
#' @param simulated_path
#' matrix containing lifetime path simulations from \code{health5_simulate_paths} function.
#'
#' @return
#' dataframe output containing mean and standard deviation of different statistics
#'
#' @export
#'
#' @examples example
health5_stats <- function (init_age, init_state, trans_probs) {

    simulated_path=health5_simulate_paths(trans_probs, init_age, init_state, cohort=50000)

if (init_state == 0) {
    total_life=health5_total_time_stats(simulated_path, 4)
    years_H=health5_total_time_stats(simulated_path, 0)
    years_M=health5_total_time_stats(simulated_path, 1)
    years_D=health5_total_time_stats(simulated_path, 2)
    years_MD=health5_total_time_stats(simulated_path, 3)
    years_disability=years_D+years_MD
    years_illness=years_M+years_MD
    first_H=health5_first_time_stats(simulated_path, 0)
    first_M=health5_first_time_stats(simulated_path, 1)
    first_D=health5_first_time_stats(simulated_path, 2)
    first_MD=health5_first_time_stats(simulated_path, 3)

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
    total_life=health5_total_time_stats(simulated_path, 4)
    years_H=health5_total_time_stats(simulated_path, 0)
    years_M=health5_total_time_stats(simulated_path, 1)
    years_D=health5_total_time_stats(simulated_path, 2)
    years_MD=health5_total_time_stats(simulated_path, 3)
    years_disability=years_D+years_MD
    years_illness=years_M+years_MD
    first_MD=health5_first_time_stats(simulated_path, 3)

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
    total_life=health5_total_time_stats(simulated_path, 4)
    years_H=health5_total_time_stats(simulated_path, 0)
    years_M=health5_total_time_stats(simulated_path, 1)
    years_D=health5_total_time_stats(simulated_path, 2)
    years_MD=health5_total_time_stats(simulated_path, 3)
    years_disability=years_D+years_MD
    years_illness=years_M+years_MD
    first_M=health5_first_time_stats(simulated_path, 1)
    first_MD=health5_first_time_stats(simulated_path, 3)

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
    total_life=health5_total_time_stats(simulated_path, 4)
    years_H=health5_total_time_stats(simulated_path, 0)
    years_M=health5_total_time_stats(simulated_path, 1)
    years_D=health5_total_time_stats(simulated_path, 2)
    years_MD=health5_total_time_stats(simulated_path, 3)
    years_disability=years_D+years_MD
    years_illness=years_M+years_MD
    first_M=health5_first_time_stats(simulated_path, 1)

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

