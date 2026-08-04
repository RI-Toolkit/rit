# simulation functions

#' Simulate cohort life path
#'
#' Simulates the path each life takes in an initial cohort using transition probabilities.
#' Supports multi-frequency transition matrices.
#'
#' @param trans_probs
#' a list of transition probability matrices, preferably generated from \code{health3_get_trans_probs}.
#' @param init_age
#' integer between 0 and closure age denoting current age.
#' @param closure_age
#' maximum life span
#' @param init_state
#' 0 for healthy, 1 for disabled
#' @param cohort
#' integer (default 10000) denoting number of people in the simulation
#'
#' @return
#' a matrix where each row represents a new individual, and the columns represent
#' the individual's movement through each state over each time step.
#'
#' -1 (death) is absorbing, so if an individual enters that state, the rest of the row will be -1.
#'
#' @noRd
#'
health3_simulate_paths <- function(trans_probs, init_age, closure_age, init_state, cohort = 10000) {
    # screening for errors
    if (init_state != 0 & init_state != 1) {
        stop('invalid state, 0 for healthy and 1 for disabled')
    }

    if (init_age < 0 | init_age > closure_age) {
        stop('invalid age')
    }

    if (as.integer(init_age) != init_age) {
        stop('initial age must be an integer')
    }

    if (cohort < 0) {
        stop('cohort needs to be a positive integer')
    }

    if (as.integer(cohort) != cohort) {
        stop('cohort needs to be an integer')
    }

    # Dynamically determine the number of simulation steps from the length of trans_probs
    num_transitions <- length(trans_probs)
    num_periods <- num_transitions + 1

    # create empty matrix to contain simulated population
    simulated_pop <- matrix(0, nrow = cohort, ncol = num_periods)

    # initialise all individuals
    simulated_pop[, 1] <- init_state

    for (i in 2:num_periods) {

        # count how many people are in each state
        n_healthy <- sum(simulated_pop[, i-1] == 0)
        n_disabled <- sum(simulated_pop[, i-1] == 1)

        # transition healthy individuals (only if there are any left)
        if (n_healthy > 0) {
            simulated_pop[simulated_pop[, i-1] == 0, i] <- sample(
                c(0, 1, -1),
                n_healthy,
                replace = TRUE,
                prob = trans_probs[[i-1]][1, ]
            )
        }

        # transition disabled individuals (only if there are any left)
        if (n_disabled > 0) {
            simulated_pop[simulated_pop[, i-1] == 1, i] <- sample(
                c(0, 1, -1),
                n_disabled,
                replace = TRUE,
                prob = trans_probs[[i-1]][2, ]
            )
        }

        # death is an absorbing state
        simulated_pop[simulated_pop[, i-1] == -1, i] <- -1
    }

    return(simulated_pop)
}
