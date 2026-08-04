# statistics functions

#' Average future lifetime
#'
#' Calculates the average future life time and its standard deviation
#' given initial state and age of an individual by simulating life time paths.
#'
#' @param model_type S for static model, T for trend model, F for frailty model
#' @param init_age integer between 0 and closure age denoting initial age
#' @param closure_age maximum life span
#' @param init_state 0 for healthy, 1 for disabled
#' @param trans_probs list of transition probability matrices
#' @param simulated_path a matrix of lifetime paths simulations
#' @param female 0 for male, 1 for female (compulsory for frailty model)
#' @param year integer indicating current year (compulsory for frailty model)
#' @param param_file parameter file (compulsory for frailty model)
#' @param n integer denoting number of unique latent factor simulations
#' @param freq integer denoting the number of transition steps per year
#'
#' @return numeric output for average and standard deviation of future lifetime
#' @export
health3_afl <- function(model_type, init_age, closure_age = 110, init_state, trans_probs = NULL, simulated_path = NULL, female = NULL, year = NULL, param_file = NULL, n = 1000, freq) {

    if (!model_type %in% c('S', 'T', 'F')) {
        stop('invalid model type, use S for static, T for trend, and F for frailty model')
    }

    if (model_type %in% c('S', 'T')) {
        if (is.null(trans_probs) & is.null(simulated_path)) {
            stop('no transition probability matrices or simulated paths were provided')
        }

        if (is.null(simulated_path)) {
            SP <- health3_simulate_paths(trans_probs, init_age, closure_age, init_state, 10000)
        } else {
            SP <- simulated_path
        }

        # Vectorized calculation for time of death (with safety for survivors)
        has_died <- rowSums(SP == -1) > 0
        future_lifetimes <- rep((ncol(SP) - 1) / freq, nrow(SP))

        if (any(has_died)) {
            death_idx <- max.col(SP == -1, ties.method = "first")
            future_lifetimes[has_died] <- (death_idx[has_died] - 1.5) / freq
        }

    } else if (model_type == 'F') {
        future_lifetimes <- numeric(n * 10000)

        for (x in 1:n) {
            TP <- health3_get_trans_probs('F', param_file, init_age, closure_age, female, year, freq)
            SP <- health3_simulate_paths(TP, init_age, closure_age, init_state, cohort = 10000)

            has_died <- rowSums(SP == -1) > 0
            fl <- rep((ncol(SP) - 1) / freq, nrow(SP))

            if (any(has_died)) {
                death_idx <- max.col(SP == -1, ties.method = "first")
                fl[has_died] <- (death_idx[has_died] - 1.5) / freq
            }

            idx_range <- ((x - 1) * 10000 + 1):(x * 10000)
            future_lifetimes[idx_range] <- fl
        }
    }

    return(list('mean' = mean(future_lifetimes), 's.dev' = stats::sd(future_lifetimes)))
}

#' Healthy Future lifetime
#'
#' Calculates the expected future lifetime spent and its standard deviation
#' in the healthy state.
#'
#' @export
health3_hfl <- function(model_type, init_age, closure_age = 110, init_state, trans_probs = NULL, simulated_path = NULL, female = NULL, year = NULL, param_file = NULL, n = 1000, freq) {

    if (!model_type %in% c('S', 'T', 'F')) {
        stop('invalid model type, use S for static, T for trend, and F for frailty model')
    }

    if (model_type %in% c('S', 'T')) {
        if (is.null(trans_probs) & is.null(simulated_path)) {
            stop('no transition probability matrices or simulated paths were provided')
        }

        if (is.null(simulated_path)) {
            SP <- health3_simulate_paths(trans_probs, init_age, closure_age, init_state, 10000)
        } else {
            SP <- simulated_path
        }

        # Vectorized calculation
        healthy_lifetimes <- rowSums(SP == 0) / freq
        if (init_state == 0) healthy_lifetimes <- healthy_lifetimes - (0.5 / freq)

    } else if (model_type == 'F') {
        healthy_lifetimes <- numeric(n * 10000)

        for (x in 1:n) {
            TP <- health3_get_trans_probs('F', param_file, init_age, closure_age, female, year, freq)
            SP <- health3_simulate_paths(TP, init_age, closure_age, init_state, cohort = 10000)

            hl <- rowSums(SP == 0) / freq
            if (init_state == 0) hl <- hl - (0.5 / freq)

            idx_range <- ((x - 1) * 10000 + 1):(x * 10000)
            healthy_lifetimes[idx_range] <- hl
        }
    }
    return(list('mean' = mean(healthy_lifetimes), 's.dev' = stats::sd(healthy_lifetimes)))
}

#' Disabled Future Lifetime
#'
#' Calculates the average future lifetime spent in disabled state and the standard
#' deviation by simulating life time paths.
#'
#' @export
health3_dfl <- function(model_type, init_age, closure_age = 110, init_state, trans_probs = NULL, simulated_path = NULL, female = NULL, year = NULL, param_file = NULL, n = 1000, freq) {

    if (!model_type %in% c('S', 'T', 'F')) {
        stop('invalid model type')
    }

    if (model_type %in% c('S', 'T')) {
        if (is.null(trans_probs) & is.null(simulated_path)) {
            stop('no transition probability matrices or simulated paths were provided')
        }

        if (is.null(simulated_path)) {
            SP <- health3_simulate_paths(trans_probs, init_age, closure_age, init_state, 10000)
        } else {
            SP <- simulated_path
        }

        disabled_lifetimes <- rowSums(SP == 1) / freq
        if (init_state == 1) disabled_lifetimes <- disabled_lifetimes - (0.5 / freq)

    } else if (model_type == 'F') {
        disabled_lifetimes <- numeric(n * 10000)

        for (x in 1:n) {
            TP <- health3_get_trans_probs('F', param_file, init_age, closure_age, female, year, freq)
            SP <- health3_simulate_paths(TP, init_age, closure_age, init_state, cohort = 10000)

            dl <- rowSums(SP == 1) / freq
            if (init_state == 1) dl <- dl - (0.5 / freq)

            idx_range <- ((x - 1) * 10000 + 1):(x * 10000)
            disabled_lifetimes[idx_range] <- dl
        }
    }
    return(list('mean' = mean(disabled_lifetimes), 's.dev' = stats::sd(disabled_lifetimes)))
}

#' Time until onset of disability (conditional on being disabled)
#'
#' Uses simulation to produce an average time until a healthy individual becomes disabled
#' during their life time.
#'
#' @export
health3_time_to_disabled <- function(model_type, init_age, closure_age = 110, init_state, trans_probs = NULL, simulated_path = NULL, female = NULL, year = NULL, param_file = NULL, n = 1000, freq) {

    if (init_state != 0) stop('initial state needs to be 0')

    if (model_type %in% c('S', 'T')) {
        if (is.null(trans_probs) & is.null(simulated_path)) stop('missing trans_probs or path')

        if (is.null(simulated_path)) {
            SP <- health3_simulate_paths(trans_probs, init_age, closure_age, 0, 10000)
        } else {
            SP <- simulated_path
        }

        has_disabled <- rowSums(SP == 1) > 0
        SP_dis <- SP[has_disabled, , drop = FALSE]

        if (nrow(SP_dis) > 0) {
            fd_idx <- max.col(SP_dis == 1, ties.method = "first")
            first_time <- (fd_idx - 1.5) / freq
        } else {
            first_time <- numeric(0)
        }

    } else if (model_type == 'F') {
        first_time <- numeric(0) # dynamic length since it's conditional

        for (x in 1:n) {
            TP <- health3_get_trans_probs('F', param_file, init_age, closure_age, female, year, freq)
            SP <- health3_simulate_paths(TP, init_age, closure_age, 0, cohort = 10000)

            has_disabled <- rowSums(SP == 1) > 0
            SP_dis <- SP[has_disabled, , drop = FALSE]

            if (nrow(SP_dis) > 0) {
                fd_idx <- max.col(SP_dis == 1, ties.method = "first")
                first_time <- c(first_time, (fd_idx - 1.5) / freq)
            }
        }
    }
    return(list('mean' = mean(first_time), 's.dev' = stats::sd(first_time)))
}

#' Survival Statistics
#'
#' Produces statistics including: total expected lifetime, healthy lifetime,
#' disabled lifetime, onset of disability (if initial state is healthy).
#'
#' @noRd
health3_survival_stats <- function(model_type, init_age, closure_age, init_state, trans_probs = NULL, simulated_path = NULL, female = NULL, year = NULL, param_file = NULL, n = 1000, freq) {

    if (model_type %in% c('S', 'T')) {

        if (is.null(simulated_path)) {
            SP <- health3_simulate_paths(trans_probs, init_age, closure_age, init_state, 10000)
        } else {
            SP <- simulated_path
        }

        # Highly optimized vectorized calculations
        has_died <- rowSums(SP == -1) > 0
        total_lifetime <- rep((ncol(SP) - 1) / freq, nrow(SP))

        if (any(has_died)) {
            death_idx <- max.col(SP == -1, ties.method = "first")
            total_lifetime[has_died] <- (death_idx[has_died] - 1.5) / freq
        }

        healthy_lifetime <- rowSums(SP == 0) / freq
        disabled_lifetime <- rowSums(SP == 1) / freq

        if (init_state == 0) {
            healthy_lifetime <- healthy_lifetime - (0.5 / freq)

            has_disabled <- rowSums(SP == 1) > 0
            first_disabled <- rep(NA, nrow(SP))
            if (any(has_disabled)) {
                fd_idx <- max.col(SP[has_disabled, , drop = FALSE] == 1, ties.method = "first")
                first_disabled[has_disabled] <- (fd_idx - 1.5) / freq
            }
        } else {
            disabled_lifetime <- disabled_lifetime - (0.5 / freq)
        }

    } else if (model_type == 'F') {

        total_lifetime <- numeric(n * 10000)
        healthy_lifetime <- numeric(n * 10000)
        disabled_lifetime <- numeric(n * 10000)
        first_disabled <- rep(NA, n * 10000)

        for (x in 1:n) {
            TP <- health3_get_trans_probs('F', param_file, init_age, closure_age, female, year, freq)
            SP <- health3_simulate_paths(TP, init_age, closure_age, init_state, cohort = 10000)

            idx_range <- ((x - 1) * 10000 + 1):(x * 10000)

            has_died <- rowSums(SP == -1) > 0
            tl <- rep((ncol(SP) - 1) / freq, nrow(SP))

            if (any(has_died)) {
                death_idx <- max.col(SP == -1, ties.method = "first")
                tl[has_died] <- (death_idx[has_died] - 1.5) / freq
            }
            total_lifetime[idx_range] <- tl

            hl <- rowSums(SP == 0) / freq
            dl <- rowSums(SP == 1) / freq

            if (init_state == 0) {
                hl <- hl - (0.5 / freq)

                has_disabled <- rowSums(SP == 1) > 0
                if (any(has_disabled)) {
                    fd_idx <- max.col(SP[has_disabled, , drop = FALSE] == 1, ties.method = "first")
                    abs_idx <- idx_range[has_disabled]
                    first_disabled[abs_idx] <- (fd_idx - 1.5) / freq
                }
            } else {
                dl <- dl - (0.5 / freq)
            }

            healthy_lifetime[idx_range] <- hl
            disabled_lifetime[idx_range] <- dl
        }
    }

    # Format final return dataframe
    if (init_state == 0) {
        means <- c(mean(total_lifetime), mean(healthy_lifetime), mean(disabled_lifetime), mean(first_disabled, na.rm = TRUE))
        sds <- c(stats::sd(total_lifetime), stats::sd(healthy_lifetime), stats::sd(disabled_lifetime), stats::sd(first_disabled, na.rm = TRUE))
        stats_df <- data.frame(
            'stats' = c('Mean years of life', 'Mean years in state H', 'Mean years in state F', 'First time entering state F'),
            'mean' = means,
            's.dev' = sds
        )
    } else {
        means <- c(mean(total_lifetime), mean(healthy_lifetime), mean(disabled_lifetime))
        sds <- c(stats::sd(total_lifetime), stats::sd(healthy_lifetime), stats::sd(disabled_lifetime))
        stats_df <- data.frame(
            'stats' = c('Mean years of life', 'Mean years in state H', 'Mean years in state F'),
            'mean' = means,
            's.dev' = sds
        )
    }

    return(stats_df)
}
