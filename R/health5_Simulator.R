#' function to get the path of set number of individuals
#'
#' @param list_trans_probs
#' a list of transition probability matrices, preferably generated from \code{get_trans_probs}.
#' @param init_age
#' the initial age of the path
#' @param closure_age
#' maximum life span
#' @param init_state
#' the initial state of all individuals (0 for H, 1 for M, 2 for D, 3 for MD)
#' @param cohort
#' the number of simulated cohorts
#'
#' @return
#' rows of individual paths in the states, 0 for H, 1 for M, 2 for D, 3 for MD, -1 for Dead
#' for each row it starts from the initial age as an input, and end at closure_age
#'
#' @noRd
#'
health5_simulate_paths <- function(list_trans_probs, init_age, closure_age, init_state, cohort) {
    # Dynamically determine the number of simulation steps from the transition matrices
    num_transitions <- length(list_trans_probs)
    num_periods <- num_transitions + 1

    # create empty matrix to contain simulated population
    simulated_pop <- matrix(0, nrow = cohort, ncol = num_periods)

    # initialise all individuals
    simulated_pop[, 1] <- init_state

    for (i in 2:num_periods) {

        # Loop through states 0 to 3 (H, M, D, MD)
        for (state in 0:3) {
            n_state <- sum(simulated_pop[, i-1] == state)

            # transition individuals only if there are any left in the state
            if (n_state > 0) {
                simulated_pop[simulated_pop[, i-1] == state, i] <- sample(
                    c(0, 1, 2, 3, -1),
                    n_state,
                    replace = TRUE,
                    prob = list_trans_probs[[i-1]][state + 1, ]
                )
            }
        }

        # death is an absorbing state
        simulated_pop[simulated_pop[, i-1] == -1, i] <- -1
    }

    return(simulated_pop)
}


#' the function to create lifetable for Static or Trend models
#'
#' @param list_trans_probs
#' a list of transition probability matrices, preferably generated from \code{{get_trans_probs}}.
#' @param init_age
#' the initial age of the path
#' @param closure_age
#' maximum life span
#' @param init_state
#' the initial state of all individuals (0=H, 1=M, 2=D, 3=MD)
#' @param cohort
#' the number of simulated cohorts
#'
#' @return
#' a life table data.frame
#' the row represents the age from the input initial age to closure age, and the columns are states H M D MD Dead
#'
#' @noRd
#'
health5_create_life_table <- function(list_trans_probs, init_age, closure_age, init_state, cohort) {
    num_transitions <- length(list_trans_probs)
    num_periods <- num_transitions + 1

    # Calculate fractional year step dynamically
    step <- (closure_age - init_age) / num_transitions

    state_status <- matrix(0, nrow = num_periods, ncol = 20)
    colnames(state_status) <- c(
        "Age", "Alive", "H", "M", "D", "MD", "Dead",
        "H_M", "H_D", "H_MD", "H_Dead", "M_MD", "M_Dead",
        "D_H", "D_M", "D_MD", "D_Dead", "MD_M", "MD_Dead", "H.M.D.MD_Dead"
    )

    # Initialize first row based on init_state
    state_status[1, "Age"] <- init_age
    state_status[1, "Alive"] <- 1

    if (init_state == 0) state_status[1, "H"] <- 1
    else if (init_state == 1) state_status[1, "M"] <- 1
    else if (init_state == 2) state_status[1, "D"] <- 1
    else if (init_state == 3) state_status[1, "MD"] <- 1

    for (i in 1:num_transitions) {
        # Matrix variables for current iteration
        P <- list_trans_probs[[i]]
        prev_states <- state_status[i, 3:7]

        # Use fast matrix multiplication to calculate expected states for the next period
        state_status[i+1, 3:7] <- prev_states %*% P

        state_status[i+1, "Alive"] <- 1 - state_status[i+1, "Dead"]
        state_status[i+1, "Age"] <- init_age + i * step

        # Detailed transition pathways
        state_status[i+1, "H_M"] <- prev_states[1] * P[1, 2]
        state_status[i+1, "H_D"] <- prev_states[1] * P[1, 3]
        state_status[i+1, "H_MD"] <- prev_states[1] * P[1, 4]
        state_status[i+1, "H_Dead"] <- prev_states[1] * P[1, 5]

        state_status[i+1, "M_MD"] <- prev_states[2] * P[2, 4]
        state_status[i+1, "M_Dead"] <- prev_states[2] * P[2, 5]

        state_status[i+1, "D_H"] <- prev_states[3] * P[3, 1]
        state_status[i+1, "D_M"] <- prev_states[3] * P[3, 2]
        state_status[i+1, "D_MD"] <- prev_states[3] * P[3, 4]
        state_status[i+1, "D_Dead"] <- prev_states[3] * P[3, 5]

        state_status[i+1, "MD_M"] <- prev_states[4] * P[4, 2]
        state_status[i+1, "MD_Dead"] <- prev_states[4] * P[4, 5]

        state_status[i+1, "H.M.D.MD_Dead"] <- state_status[i+1, "H_Dead"] +
            state_status[i+1, "M_Dead"] +
            state_status[i+1, "D_Dead"] +
            state_status[i+1, "MD_Dead"]
    }

    # Scale probabilities by cohort size
    state_status[, 2:20] <- state_status[, 2:20] * cohort
    return(as.data.frame(state_status))
}


#' the function to get n_sim number of simulated lifetables for Frailty model
#'
#' @param model_type
#' choose F for Frailty model
#' @param param_file
#' matrix of estimated parameters to construct the five state model.
#' @param female
#' female 1 if female, 0 if male
#' @param wave_index
#' the wave index = (interview year - 1998)/2 + 1
#' @param latent
#' initial value of latent factor, normally take the value 0
#' @param init_age
#' the initial age of the life table
#' @param closure_age
#' maximum life span
#' @param init_state
#' 0 for H state, 1 for M state, 2 for D state, 3 for MD state
#' @param n_sim
#' the number of simulations
#' @param cohort
#' number of people at the beginning of the life table
#' @param mean
#' TRUE to return expected life table, FALSE to return all simulated life tables
#' @param freq
#' integer denoting the number of transition steps per year
#'
#' @return
#' a list of n_sim number of life table matrices when mean=FALSE
#' or the mean life table when mean=TRUE
#'
#' @noRd
#'
health5_simulate_life_table <- function(model_type, param_file, female, wave_index, latent, init_age, closure_age, init_state, n_sim, cohort, mean, freq) {

    if (model_type != 'F') {
        stop('use frailty model to simulate lifetables')
    }

    # Pre-allocate list to prevent massive memory re-allocation in R
    state_status_full <- vector("list", n_sim)

    for (i in seq_len(n_sim)) {
        list_trans_probs <- health5_get_trans_probs(model_type, param_file, init_age, closure_age, female, wave_index, latent, freq)
        state_status_full[[i]] <- health5_create_life_table(list_trans_probs, init_age, closure_age, init_state, cohort)
    }

    if (mean) {
        # Reduce handles the dataframe addition cleanly
        return(Reduce('+', state_status_full) / n_sim)
    } else {
        return(state_status_full)
    }
}
