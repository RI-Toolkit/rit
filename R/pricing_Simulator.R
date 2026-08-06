###############################################################################
###### MAIN SIMULATION WRAPPER

#' Policy Cashflow Simulator
#'
#' Simulate cash flows using Monte-Carlo methods for various policies
#'
#' @name simulate_cf
#' @param policy
#' Policy type to simulate:
#' `policy object created using a create_policy function
#' @param init_age
#' Initial age of the policyholder in years
#' @param seed
#' Seed for random generator
#' @param n
#' Number of paths to simulate (Monte-Carlo method)
#' @param state
#' Simulated state matrix via Health-State / Aggregate Mortality
#' @param econ_var
#' Simulated economic variables via Economic Scenario Generator
#' @param cohort_death_probs
#' Simulated death probabilities for cohort (for Pooled Annuity)
#' @param frequency
#' string selecting the simulation frequency: "year", "quarter", or "month". Default is "year".
#' @return
#' Matrix of cash flow vectors for each simulated path
#' @export simulate_cf
simulate_cf <- function(policy, init_age = 65, seed = NULL, n = 1000, state = NULL, econ_var = NULL, cohort_death_probs = NULL) {

    frequency = policy$frequency[1]

    # Set cash flow function based on input policy
    cf_func <- switch(policy$name[1],
                      "AP" = cf_account_based_pension,
                      "RM" = cf_reverse_mortgage,
                      "VA" = cf_variable_annuity,
                      "PA" = cf_pooled_annuity,
                      "CA" = cf_care_annuity,
                      "LA" = cf_life_annuity)

    # If not provided, get states for each path (matrix)
    if (is.null(state)) {
        state <- get_state_simulation(policy, age = init_age, female = 1, seed, n, frequency = frequency)
    }

    # Validate formatting of mortality state data
    if (nrow(state) != n) {
        stop("Error: State matrix does not fit number of paths requested")
    }

    period <- ncol(state)

    # If not provided, get economic data for each path (list of matrices)
    if (is.null(econ_var)) {
        econ_var <- get_econ_simulation(state, n, seed, frequency = frequency)

        econ_var <- list(stock=econ_var$market_index, infla=econ_var$inflation_index,
                         zcp3m=econ_var$zcp3m_yield, house=econ_var$home_index,
                         sdf=econ_var$discount_factors)
    }

    # Validate formatting of economic data
    for (i in names(econ_var)) {
        if (ncol(econ_var[[i]]) < ncol(state)) {
            error_message <- paste("Error: Reduced dimension size in economic simulation -> $", i,
                                   " (expected ", nrow(state), "x", ncol(state),
                                   ", got ", nrow(econ_var[[i]]), "x", ncol(econ_var[[i]]), ")",
                                   sep="")
            stop(error_message)
        }
    }

    # Get matrix of economic variables for each path (passing frequency down)
    data <- get_policy_scenario(policy, age = init_age, female = 1, seed, n, period, econ_var, cohort_death_probs, frequency = frequency)

    # Initialize output matrix
    cf <- matrix(nrow = n, ncol = ncol(state))

    # Generate cash flows for each state vector
    for (i in seq(1, n)) {
        cf[i,] <- cf_func(policy, state[i,], data[[i]])
    }

    # Round cashflows to cents
    cf <- round(cf, 2)

    result <- list(cf = cf, sdf = unname(econ_var$sdf))

    return(result)
}


###############################################################################
###### POLICY SCENARIO FUNCTION

#' Scenario Generator
#'
#' Generates and encapsulates scenario data for a given policy
#'
#' @param policy
#' Policy object containing necessary parameters (see create_policy_ )
#' @param age
#' Initial age of policyholder in years
#' @param female
#' Gender of policyholder, 0 for male, 1 for female
#' @param seed
#' Seed for random generator
#' @param n
#' Number of paths to simulate (Monte-Carlo method)
#' @param period
#' Number of periods to simulate
#' @param econ_var
#' Simulated economic variables via Economic Scenario Generator
#' @param death_probs
#' Simulated death probabilities for cohort
#' @param frequency
#' string selecting the simulation frequency: "year", "quarter", or "month".
#'
#' @return
#' Data frame containing all variables generated using other modules
get_policy_scenario <- function(policy, age, female, seed, n, period, econ_var, death_probs, frequency) {

    var_sim <- econ_var

    if (policy$name[1] == "AP") {

        # Get all relevant economic variables
        infla <- get_inflation_rate(var_sim)
        stock <- get_stock_return(var_sim)

        # Organise economic inputs into a data.frame for each path
        data <- list()
        for (i in seq(1, n)) {
            temp <- data.frame(infla = infla[i, ],
                               stock = stock[i, ])
            data <- append(data, list(temp))
        }

    } else if (policy$name[1] == "CA" | policy$name[1] == "LA") {

        # Get all relevant economic variables
        #infla <- get_inflation_rate(var_sim)

        # Organise economic inputs into a data.frame for each path
        data <- list()
        for (i in seq(1, n)) {
            #temp <- data.frame(infla = infla[i, ])
            temp <- data.frame()
            data <- append(data, list(temp))
        }

    } else if (policy$name[1] == "PA") {

        # Get all relevant health variables for pool (passing frequency down)
        pool_r <- get_pool_realised(age, female, seed, n, policy$size, death_probs, frequency = frequency)
        pool_e <- get_pool_expected(age, female, seed, policy$size, death_probs, frequency = frequency)

        # Get all relevant economic variables
        stock <- get_stock_return(var_sim)
        stock <- stock[, 1:length(pool_e)]

        # Organise economic inputs into a data.frame for each path
        data <- list()
        for (i in seq(1, n)) {
            temp <- data.frame(pool_r = pool_r[i, ],
                               pool_e = pool_e,
                               stock = stock[i, ])
            data <- append(data, list(temp))
        }

    } else if (policy$name[1] == "RM") {

        # Get all relevant economic variables
        zcp3m <- get_zcp3m_yield(var_sim)
        house <- get_house_return(var_sim)

        # Organise economic inputs into a data.frame for each path
        data <- list()
        for (i in seq(1, n)) {
            temp <- data.frame(house = house[i, ],
                               zcp3m = zcp3m[i, ])
            data <- append(data, list(temp))
        }

    } else if (policy$name[1] == "VA") {

        # Get all relevant economic / health variables
        stock <- get_stock_return(var_sim)

        # Organise economic inputs into a data.frame for each path
        data <- list()
        for (i in seq(1, n)) {
            temp <- data.frame(stock = stock[i, ])
            data <- append(data, list(temp))
        }

    } else {
        stop("Error: invalid policy object")
    }

    return(data)

}

get_econ_simulation <- function(state, n, seed, frequency) {
    esg_names <- c("ASX200", "CPI", "home_index", "zcp3m_yield", "discount_factors")
    # Generalize naming to match generic pricing module inputs
    gen_names <- c("market_index", "inflation_index", "home_index", "zcp3m_yield", "discount_factors")
    simulated_vars <- esg_var_simulator(ncol(state), n, frequency = frequency, return_sdf = TRUE, seed = seed)
    filtered_vars <- simulated_vars[esg_names]
    names(filtered_vars) <- gen_names
    return(filtered_vars)
}

###############################################################################
###### STATE SIMULATION FUNCTION

get_state_simulation <- function(policy, age, female, seed, n, frequency) {
    if (policy$name[1] == "CA") {
        if (nrow(policy) == 2) {
            probs <- get_trans_probs(3, 'S', rit::US_HRS_3, init_age = age, closure_age = 110, female = 1, frequency = frequency)
        } else if (nrow(policy) == 4) {
            probs <- get_trans_probs(5, 'S', rit::US_HRS_5, init_age = age, closure_age = 110, female = 1, frequency = frequency)
        } else {
            stop("Error: CA policy object needs to have 2 or 4 rows")
        }
        return(simulate_health_state_paths(probs, init_age = age, closure_age = 110, cohort = n))
    } else if (policy$name[1] == "RM") {
        probs <- get_trans_probs(3, 'S', rit::US_HRS_3, init_age = age, closure_age = 110, female = 1, frequency = frequency)
        return(simulate_health_state_paths(probs, init_age = age, closure_age = 110, cohort = n))
    } else {
        # Pass frequency down to get_aggregate_mortality
        return(get_aggregate_mortality(age, female, seed, n, frequency = frequency))
    }
}

###############################################################################
###### PLACEHOLDER FUNCTIONS


# ------------------------------------------------------------------------
# ---- Health State Module

# get_health_state_3 <- function(age, female = 1, seed = 0, n = 1000) {
#     trans_probs <-  get_trans_probs(3, 'T', rit::US_HRS_3, age, closure_age = 110, (female = 1), year = 2022)
#     return(simulate_health_state_paths(trans_probs, age, 0, closure_age = 110, n))
# }
#
# get_health_state_5 <- function(age, female = 1, seed = 0, n = 1000) {
#     trans_probs <-  get_trans_probs(5, 'T', rit::US_HRS_5, age, closure_age = 110, (female = 1), year = 2022, latent = 0)
#     return(simulate_health_state_paths(trans_probs, age, 0, closure_age = 110, n))
# }

# ------------------------------------------------------------------------
# ---- Aggregate Mortality Module



get_aggregate_mortality <- function(age, female = 1, seed = 0, n = 1000, frequency) {
    utils::capture.output(suppressWarnings(
        mortality <- sim_indiv_path(init_age = age, female = female, death_probs = NULL, closure_age = 110, n_sim = n, seed = seed, frequency = frequency)
    ))
    return(mortality)
}

get_pool_realised <- function(age, female = 1, seed = 0, n = 1000, cohort = 1000, death_probs = NULL, frequency) {

    closure_age <- 110
    if (!is.null(death_probs)) {
        closure_age <- age + length(death_probs) - 1
    }

    utils::capture.output(suppressWarnings(
        pool <- sim_cohort_path_realised(init_age = age, female = female, death_probs = death_probs, closure_age = closure_age, cohort = cohort, n_sim = n, seed = seed, frequency = frequency)
    ))

    return(pool)
}

get_pool_expected <- function(age, female = 1, seed = 0, cohort = 1000, death_probs = NULL, frequency) {

    closure_age <- 110
    if (!is.null(death_probs)) {
        closure_age <- age + length(death_probs) - 1
    }

    utils::capture.output(suppressWarnings(
        pool <- sim_cohort_path_expected(init_age = age, female = female, death_probs = death_probs, closure_age = closure_age, cohort = cohort, frequency = frequency)
    ))

    return(pool)
}

# ------------------------------------------------------------------------
# ---- Economic Scenario Generator Module

get_perc_change <- function(df) {
    result <- df
    for (i in seq(1, NCOL(df) - 1)) {
        result[,i] <- (df[,i + 1]/df[,i]) - 1
    }
    result[, ncol(df)] <- result[, ncol(df) - 1]
    return(result)
}

get_zcp3m_yield <- function(var_sim) {
    return((unname(var_sim$zcp3m)))
}

get_inflation_rate <- function(var_sim) {
    cpi <- (unname(var_sim$infla))
    return(get_perc_change(cpi))
}

get_house_return <- function(var_sim) {
    home_index <- (unname(var_sim$house))
    return(get_perc_change(home_index))
}

get_stock_return <- function(var_sim) {

    # Check the names available in the provided var_sim list
    # and extract the correct one dynamically
    if ("stock" %in% names(var_sim)) {
        asx <- unname(var_sim$stock)
    } else if ("ASX" %in% names(var_sim)) {
        asx <- unname(var_sim$ASX)
    } else if ("ASX200" %in% names(var_sim)) {
        asx <- unname(var_sim$ASX200)
    } else if ("market_index" %in% names(var_sim)) {
        asx <- unname(var_sim$market_index)
    } else {
        stop("Error: Could not locate stock data. Ensure econ_var contains 'stock' or 'ASX'.")
    }

    return(get_perc_change(asx))
}
