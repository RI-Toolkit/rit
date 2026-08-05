#' function to get the first time leaving or entering different states for a number of individuals
#'
#' @param model_type S for static model, T for trend model, F for frailty model
#' @param state 0 for leaving H, 1 for entering M, 2 for entering D, 3 for entering MD, -1 for entering Dead
#' @param init_age integer between 0 and 110 denoting initial age of individual
#' @param closure_age maximum life span
#' @param init_state 0 for healthy, 1 for disabled
#' @param trans_probs a list of transition probability matrices
#' @param simulated_path the simulated path of individuals
#' @param female 0 for male, 1 for female
#' @param year integer indicating current year
#' @param wave_index the wave index = (interview year - 1998)/2 + 1
#' @param latent initial value of latent factor
#' @param param_file parameter file
#' @param n integer denoting number of unique latent factor simulations
#' @param frequency string selecting the simulation frequency: "year", "quarter", or "month".
#'
#' @return a column that consists the first time leaving or entering the state
#'
#' @export
health5_first_time_stats <- function(model_type, state, init_age, closure_age, init_state, trans_probs = NULL, simulated_path = NULL, female = NULL, year = NULL, wave_index = NULL, latent = NULL, param_file = NULL, n = 1000, frequency) {

    # Map the string frequency to numeric steps per year
    if (frequency == "year") {
        freq <- 1
    } else if (frequency == "quarter") {
        freq <- 4
    } else if (frequency == "month") {
        freq <- 12
    } else {
        stop("Frequency must be one of 'year', 'quarter', and 'month'.")
    }

    # Internal vectorized helper
    calc_first_time <- function(SP) {
        if (state == 0) {
            logical_mat <- (SP != 0)
        } else {
            logical_mat <- (SP == state)
        }

        has_state <- rowSums(logical_mat) > 0
        res <- rep(NA_real_, nrow(SP))
        if (any(has_state)) {
            idx <- max.col(logical_mat[has_state, , drop = FALSE], ties.method = "first")
            res[has_state] <- pmax(0, (idx - 1.5) / freq)
        }
        return(res)
    }

    if (model_type %in% c('S', 'T')) {
        if (is.null(trans_probs) & is.null(simulated_path)) stop('missing matrices or path')

        if (is.null(simulated_path)) {
            SP <- health5_simulate_paths(trans_probs, init_age, closure_age, init_state, cohort = 10000)
        } else {
            SP <- simulated_path
        }
        return(matrix(calc_first_time(SP), ncol = 1))

    } else if (model_type == 'F') {
        first_time <- numeric(n * 10000)
        for (x in 1:n) {
            TP <- health5_get_trans_probs(model_type, param_file, init_age, closure_age, female, wave_index, latent, freq)
            SP <- health5_simulate_paths(TP, init_age, closure_age, init_state, cohort = 10000)
            first_time[((x - 1) * 10000 + 1):(x * 10000)] <- calc_first_time(SP)
        }
        return(matrix(first_time, ncol = 1))
    }
}


#' function to get the total time in different states for a number of individuals
#'
#' @param state 0=H, 1=M, 2=D, 3=MD, -1=Dead, 4=Alive (not Dead)
#' @export
health5_total_time_stats <- function(model_type, state, init_age, closure_age, init_state, trans_probs = NULL, simulated_path = NULL, female = NULL, year = NULL, wave_index = NULL, latent = NULL, param_file = NULL, n = 1000, frequency) {

    # Map the string frequency to numeric steps per year
    if (frequency == "year") {
        freq <- 1
    } else if (frequency == "quarter") {
        freq <- 4
    } else if (frequency == "month") {
        freq <- 12
    } else {
        stop("Frequency must be one of 'year', 'quarter', and 'month'.")
    }

    # Internal vectorized helper
    calc_total_time <- function(SP) {
        if (state == 4) {
            # Total time alive (with right-censoring safety)
            time <- rep((ncol(SP) - 1) / freq, nrow(SP))
            has_died <- rowSums(SP == -1) > 0
            if (any(has_died)) {
                death_idx <- max.col(SP == -1, ties.method = "first")
                time[has_died] <- pmax(0, (death_idx[has_died] - 1.5) / freq)
            }
            return(time)
        } else if (state == -1) {
            # Total time dead
            time <- rowSums(SP == -1) / freq
            has_died <- time > 0
            adj <- ifelse(SP[, 1] == -1, 1.0 / freq, 0.5 / freq)
            time[has_died] <- time[has_died] - adj[has_died]
            return(time)
        } else {
            # Total time in impaired/healthy states
            time <- rowSums(SP == state) / freq
            adj <- ifelse(SP[, 1] == state, 0.5 / freq, 0)
            return(time - adj)
        }
    }

    if (model_type %in% c('S', 'T')) {
        if (is.null(trans_probs) & is.null(simulated_path)) stop('missing matrices or path')

        if (is.null(simulated_path)) {
            SP <- health5_simulate_paths(trans_probs, init_age, closure_age, init_state, cohort = 10000)
        } else {
            SP <- simulated_path
        }
        return(matrix(calc_total_time(SP), ncol = 1))

    } else if (model_type == 'F') {
        total_time <- numeric(n * 10000)
        for (x in 1:n) {
            TP <- health5_get_trans_probs(model_type, param_file, init_age, closure_age, female, wave_index, latent, freq)
            SP <- health5_simulate_paths(TP, init_age, closure_age, init_state, cohort = 10000)
            total_time[((x - 1) * 10000 + 1):(x * 10000)] <- calc_total_time(SP)
        }
        return(matrix(total_time, ncol = 1))
    }
}


#' function to produce the mean and variance of a list of values
#'
#' @noRd
health5_stats_produce <- function(input) {
    output <- matrix(nrow = 1, ncol = 2)
    colnames(output) <- c('expected_value', 'st_dev')
    output[1] <- mean(input, na.rm = TRUE)
    # Divide standard deviation by sqrt(n) for valid non-NA values
    n_valid <- sum(!is.na(input))
    output[2] <- stats::sd(input, na.rm = TRUE) / sqrt(n_valid)
    return(output)
}


#' Survival Statistics
#'
#' Produces statistics for 5-state model.
#'
#' @noRd
health5_stats <- function(model_type, init_age, closure_age, init_state, trans_probs = NULL, simulated_path = NULL, female = NULL, year = NULL, wave_index = NULL, latent = NULL, param_file = NULL, n = 1000, frequency) {

    # Map the string frequency to numeric steps per year
    if (frequency == "year") {
        freq <- 1
    } else if (frequency == "quarter") {
        freq <- 4
    } else if (frequency == "month") {
        freq <- 12
    } else {
        stop("Frequency must be one of 'year', 'quarter', and 'month'.")
    }

    # Inline evaluation logic to prevent re-simulating the paths over and over
    calc_first <- function(SP, st) {
        logical_mat <- if (st == 0) (SP != 0) else (SP == st)
        has_state <- rowSums(logical_mat) > 0
        res <- rep(NA_real_, nrow(SP))
        if (any(has_state)) {
            idx <- max.col(logical_mat[has_state, , drop = FALSE], ties.method = "first")
            res[has_state] <- pmax(0, (idx - 1.5) / freq)
        }
        return(res)
    }

    calc_total <- function(SP, st) {
        if (st == 4) {
            time <- rep((ncol(SP) - 1) / freq, nrow(SP))
            has_died <- rowSums(SP == -1) > 0
            if (any(has_died)) {
                death_idx <- max.col(SP == -1, ties.method = "first")
                time[has_died] <- pmax(0, (death_idx[has_died] - 1.5) / freq)
            }
            return(time)
        } else if (st == -1) {
            time <- rowSums(SP == -1) / freq
            has_died <- time > 0
            adj <- ifelse(SP[, 1] == -1, 1.0 / freq, 0.5 / freq)
            time[has_died] <- time[has_died] - adj[has_died]
            return(time)
        } else {
            time <- rowSums(SP == st) / freq
            adj <- ifelse(SP[, 1] == st, 0.5 / freq, 0)
            return(time - adj)
        }
    }

    if (model_type %in% c('S', 'T')) {
        if (is.null(trans_probs) & is.null(simulated_path)) stop('missing matrices or path')
        SP <- if (is.null(simulated_path)) health5_simulate_paths(trans_probs, init_age, closure_age, init_state, cohort = 10000) else simulated_path

        total_life <- calc_total(SP, 4)
        years_H <- calc_total(SP, 0)
        years_M <- calc_total(SP, 1)
        years_D <- calc_total(SP, 2)
        years_MD <- calc_total(SP, 3)
        first_H <- calc_first(SP, 0)
        first_M <- calc_first(SP, 1)
        first_D <- calc_first(SP, 2)
        first_MD <- calc_first(SP, 3)

    } else if (model_type == 'F') {
        N_tot <- n * 10000
        total_life <- numeric(N_tot)
        years_H <- numeric(N_tot); years_M <- numeric(N_tot)
        years_D <- numeric(N_tot); years_MD <- numeric(N_tot)
        first_H <- numeric(N_tot); first_M <- numeric(N_tot)
        first_D <- numeric(N_tot); first_MD <- numeric(N_tot)

        for (x in 1:n) {
            TP <- health5_get_trans_probs(model_type, param_file, init_age, closure_age, female, wave_index, latent, freq)
            SP <- health5_simulate_paths(TP, init_age, closure_age, init_state, cohort = 10000)

            idx <- ((x - 1) * 10000 + 1):(x * 10000)
            total_life[idx] <- calc_total(SP, 4)
            years_H[idx] <- calc_total(SP, 0)
            years_M[idx] <- calc_total(SP, 1)
            years_D[idx] <- calc_total(SP, 2)
            years_MD[idx] <- calc_total(SP, 3)
            first_H[idx] <- calc_first(SP, 0)
            first_M[idx] <- calc_first(SP, 1)
            first_D[idx] <- calc_first(SP, 2)
            first_MD[idx] <- calc_first(SP, 3)
        }
    }

    years_disability <- years_D + years_MD
    years_illness <- years_M + years_MD

    n_tot <- length(total_life)

    # Format the return structure based on the initial state and calculate SE
    if (init_state == 0) {
        means <- c(mean(total_life), mean(years_H), mean(years_M), mean(years_D), mean(years_MD), mean(years_disability), mean(years_illness), mean(first_H, na.rm=TRUE), mean(first_M, na.rm=TRUE), mean(first_D, na.rm=TRUE), mean(first_MD, na.rm=TRUE))
        sds <- c(
            stats::sd(total_life) / sqrt(n_tot),
            stats::sd(years_H) / sqrt(n_tot),
            stats::sd(years_M) / sqrt(n_tot),
            stats::sd(years_D) / sqrt(n_tot),
            stats::sd(years_MD) / sqrt(n_tot),
            stats::sd(years_disability) / sqrt(n_tot),
            stats::sd(years_illness) / sqrt(n_tot),
            stats::sd(first_H, na.rm=TRUE) / sqrt(sum(!is.na(first_H))),
            stats::sd(first_M, na.rm=TRUE) / sqrt(sum(!is.na(first_M))),
            stats::sd(first_D, na.rm=TRUE) / sqrt(sum(!is.na(first_D))),
            stats::sd(first_MD, na.rm=TRUE) / sqrt(sum(!is.na(first_MD)))
        )
        stats_df <- data.frame(
            'stats' = c('Mean years of life', 'Mean years in state H', 'Mean years in state M','Mean years in state D','Mean years in state MD', 'Mean years with disability','Mean years with illness','First time leaving state H','First time entering state M', 'First time entering state D', 'First time entering state MD'),
            'mean' = means, 's.dev' = sds
        )
    } else if (init_state == 1) {
        means <- c(mean(total_life), mean(years_H), mean(years_M), mean(years_D), mean(years_MD), mean(years_disability), mean(years_illness), mean(first_MD, na.rm=TRUE))
        sds <- c(
            stats::sd(total_life) / sqrt(n_tot),
            stats::sd(years_H) / sqrt(n_tot),
            stats::sd(years_M) / sqrt(n_tot),
            stats::sd(years_D) / sqrt(n_tot),
            stats::sd(years_MD) / sqrt(n_tot),
            stats::sd(years_disability) / sqrt(n_tot),
            stats::sd(years_illness) / sqrt(n_tot),
            stats::sd(first_MD, na.rm=TRUE) / sqrt(sum(!is.na(first_MD)))
        )
        stats_df <- data.frame(
            'stats' = c('Mean years of life', 'Mean years in state H', 'Mean years in state M','Mean years in state D','Mean years in state MD', 'Mean years with disability','Mean years with illness','First time entering state MD'),
            'mean' = means, 's.dev' = sds
        )
    } else if (init_state == 2) {
        means <- c(mean(total_life), mean(years_H), mean(years_M), mean(years_D), mean(years_MD), mean(years_disability), mean(years_illness), mean(first_M, na.rm=TRUE), mean(first_MD, na.rm=TRUE))
        sds <- c(
            stats::sd(total_life) / sqrt(n_tot),
            stats::sd(years_H) / sqrt(n_tot),
            stats::sd(years_M) / sqrt(n_tot),
            stats::sd(years_D) / sqrt(n_tot),
            stats::sd(years_MD) / sqrt(n_tot),
            stats::sd(years_disability) / sqrt(n_tot),
            stats::sd(years_illness) / sqrt(n_tot),
            stats::sd(first_M, na.rm=TRUE) / sqrt(sum(!is.na(first_M))),
            stats::sd(first_MD, na.rm=TRUE) / sqrt(sum(!is.na(first_MD)))
        )
        stats_df <- data.frame(
            'stats' = c('Mean years of life', 'Mean years in state H', 'Mean years in state M','Mean years in state D','Mean years in state MD', 'Mean years with disability','Mean years with illness','First time entering state M', 'First time entering state MD'),
            'mean' = means, 's.dev' = sds
        )
    } else if (init_state == 3) {
        means <- c(mean(total_life), mean(years_H), mean(years_M), mean(years_D), mean(years_MD), mean(years_disability), mean(years_illness), mean(first_M, na.rm=TRUE))
        sds <- c(
            stats::sd(total_life) / sqrt(n_tot),
            stats::sd(years_H) / sqrt(n_tot),
            stats::sd(years_M) / sqrt(n_tot),
            stats::sd(years_D) / sqrt(n_tot),
            stats::sd(years_MD) / sqrt(n_tot),
            stats::sd(years_disability) / sqrt(n_tot),
            stats::sd(years_illness) / sqrt(n_tot),
            stats::sd(first_M, na.rm=TRUE) / sqrt(sum(!is.na(first_M)))
        )
        stats_df <- data.frame(
            'stats' = c('Mean years of life', 'Mean years in state H', 'Mean years in state M','Mean years in state D','Mean years in state MD', 'Mean years with disability','Mean years with illness','First time entering state M'),
            'mean' = means, 's.dev' = sds
        )
    }

    return(stats_df)
}
