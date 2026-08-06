# ---------------------------------------------------------------------------- #
# --------------------------- Account Based Pension -------------------------- #
# ---------------------------------------------------------------------------- #

#' Account Based Pension - Cashflow Simulator
#'
#' Simulate cash flows for Account Based Pension
#' @param policy
#' Policy object containing necessary parameters (see create_policy_AP)
#' @param state
#' State vector containing state values for entire duration
#' @param data
#' Data frame containing all variables generated using other modules
#'
#' @return
#' Vector of cashflows for at each time point
cf_account_based_pension <- function(policy, state, data) {

    # Extract relevant policy variables
    balance <- policy$bal[1]
    frequency <- policy$frequency[1]
    rate <- policy$rate

    # Map string frequency to numeric periods per year
    if (frequency == "year") {
        freq_num <- 1
    } else if (frequency == "quarter") {
        freq_num <- 4
    } else if (frequency == "month") {
        freq_num <- 12
    } else {
        stop("Frequency must be one of 'year', 'quarter', and 'month'.")
    }

    # Initialize output vector
    cf <- rep(0, times = length(state))

    i <- 1
    while (state[i] != -1 & i < length(state)) {    # while PH is not dead

        # Determine which year of the policy we are in to get the correct annual rate
        year_idx <- floor((i - 1) / freq_num) + 1

        # Safety check: prevent out-of-bounds if state sequence outlasts rate array
        if (year_idx > length(rate)) {
            year_idx <- length(rate)
        }

        # Apply market returns for the fractional period
        balance <- balance * (1 + data$stock[i])

        # Record cashflow to output vector (fraction of the annual rate)
        cf[i] <- balance * (rate[year_idx] / freq_num)

        # Update balance
        balance <- balance - cf[i]

        i <- i + 1
    }

    # Withdraw remaining balance after death (payout to family)
    cf[i] <- balance

    return(cf)
}


# ---------------------------------------------------------------------------- #
# ------------------------------- Care Annuity ------------------------------- #
# ---------------------------------------------------------------------------- #


#' Care Annuity - Cashflow Simulator (Matrix Optimized)
#'
#' @param policy Policy object containing necessary parameters (benefit, increase, min, defer)
#' @param state State matrix or vector (0=H, 1=M, 2=D, 3=MD, -1=Dead)
#' @param data Data frame containing all variables generated using other modules
#'
#' @return Matrix or Vector of cashflows at each time point
cf_care_annuity <- function(policy, state, data = NULL) {

    # Extract relevant policy variables
    increase <- policy$increase
    benefit <- policy$benefit
    minimum <- policy$min

    # Extract deferral/waiting period. Default to 0 if not explicitly defined.
    defer <- policy$defer
    if (is.null(defer)) defer <- 0

    # Check if simulation passed a 2D matrix (all paths at once) or a 1D vector
    if (is.matrix(state)) {
        cf <- matrix(0, nrow = nrow(state), ncol = ncol(state))

        # Extract column indices to represent periods (t = 1, 2, 3...)
        t_mat <- col(state)

        # Growth factor compounded per period
        growth_factor <- (1 + increase)^(t_mat - 1)

        is_alive <- (state != -1)

        # Check if the current period is strictly after the waiting period
        after_waiting <- (t_mat > defer)

        # The guarantee period starts after the deferral period ends.
        in_guarantee <- (t_mat <= (defer + minimum[1]))

        # 1. Base benefit
        eligible_for_base <- (is_alive | in_guarantee) & after_waiting
        if (any(eligible_for_base)) {
            cf[eligible_for_base] <- benefit[1] * growth_factor[eligible_for_base]
        }

        # 2. Additional impaired benefits
        impaired <- is_alive & (state > 0) & after_waiting
        if (any(impaired)) {
            cf[impaired] <- cf[impaired] + benefit[state[impaired] + 1] * growth_factor[impaired]
        }

        return(cf)

    } else {
        # Fallback for 1D vector (single path execution)
        t <- seq_along(state)

        # Growth factor compounded per period
        growth_factor <- (1 + increase)^(t - 1)
        cf <- numeric(length(state))

        is_alive <- (state != -1)

        after_waiting <- (t > defer)
        in_guarantee <- (t <= (defer + minimum[1]))

        eligible_for_base <- (is_alive | in_guarantee) & after_waiting
        if (any(eligible_for_base)) {
            cf[eligible_for_base] <- benefit[1] * growth_factor[eligible_for_base]
        }

        impaired <- is_alive & (state > 0) & after_waiting
        if (any(impaired)) {
            cf[impaired] <- cf[impaired] + benefit[state[impaired] + 1] * growth_factor[impaired]
        }

        return(cf)
    }
}

# ---------------------------------------------------------------------------- #
# ------------------------------- Life Annuity ------------------------------- #
# ---------------------------------------------------------------------------- #


#' Life Annuity - Cashflow Simulator
#'
#' @param policy
#' Policy object containing necessary parameters (see create_policy_LA)
#' @param state
#' State vector containing state values for entire duration (-1 = dead)
#' @param data
#' Data frame containing all variables generated using other modules
#'
#' @return
#' Vector of cashflows at each time point
cf_life_annuity <- function(policy, state, data) {

    # Extract relevant policy variables
    # (Per documentation, these are already in frequency units)
    increase <- policy$increase[1]
    benefit  <- policy$benefit[1]
    d        <- policy$defer[1]

    # Initialize output vector
    cf <- numeric(length(state))

    # Identify all periods where the policyholder is alive
    alive <- state != -1

    if (!any(alive)) {
        return(cf)
    }

    # Create a vector of step indices (1, 2, 3... length(state))
    i <- seq_along(state)

    # Calculate the compounded benefit for all alive periods.
    # This exactly mirrors your original logic: benefit * (1+increase) per step
    cf[alive] <- benefit * (1 + increase)^i[alive]

    # Enforce the deferment period by zeroing out early cashflows
    cf[i <= d] <- 0

    return(cf)
}


# ---------------------------------------------------------------------------- #
# ------------------------------ Pooled Annuity ------------------------------ #
# ---------------------------------------------------------------------------- #


#' Pooled Annuity - Cashflow Simulator
#'
#' @param policy
#' Policy object containing necessary parameters (see create_policy_PA)
#' @param state
#' State vector containing state values for entire duration
#' @param data
#' Data frame containing all variables generated using other modules
#'
#' @return
#' Vector of cashflows for at each time point
cf_pooled_annuity <- function(policy, state, data) {

    # Extract relevant policy variables
    size <- policy$size
    benefit <- policy$benefit
    interest <- policy$interest

    # Initialize output vector
    cf <- rep(0, times = length(state))

    i <- 1
    while (state[i] != -1 & i < length(state)) {     # while PH is not dead

        # Get benefit if alive
        cf[i] <- benefit

        # Skip rate calculations for final year
        if (i + 1 > length(state)) break

        # Calculate expected and realized survivorship rates
        suv_e <- data$pool_e[i + 1] / data$pool_e[i]
        suv_r <- data$pool_r[i + 1] / data$pool_r[i]

        if (suv_r == 0) break

        # Mortality experience adjustment factor
        mea <- suv_e / suv_r

        # Interest rate adjustment factor
        ira <- (1 + data$stock[i]) / (1 + interest)

        # Scale benefit from period t to t + 1
        benefit <- benefit * mea * ira

        i <- i + 1
    }

    return(cf)
}


# ---------------------------------------------------------------------------- #
# ----------------------------- Reverse Mortgage ----------------------------- #
# ---------------------------------------------------------------------------- #

#' Reverse Mortgage - Cashflow Simulator
#'
#' @param policy
#' Policy object containing necessary parameters (see create_policy_RM)
#' @param state
#' State vector containing state values for entire duration
#' @param data
#' Data frame containing all variables generated using other modules
#'
#' @return
#' Vector of cashflows for at each time point
cf_reverse_mortgage <- function(policy, state, data) {

    # Extract relevant policy variables
    LVR <- policy$LVR
    cost <- policy$trans_cost
    value <- policy$value
    margin <- policy$margin

    # Initialize output vector
    cf <- rep(0, times = length(state))

    # Get loan amount for policyholder
    loan <- LVR * value
    #cf[1] <- loan      # only value NNEG!

    i <- 1
    while (state[i] == 0 & i < length(state)) {     # while PH is healthy

        # Compound loan value over 1 year period (excluding lending margin)
        loan <- loan * (1 + data$zcp3m[i])

        # Update house value after 1 year period
        value <- value * (1 + data$house[i])

        i <- i + 1
    }

    # Add excess from lending margin to loan
    loan <- loan * exp((i - 1) * margin)

    # Calculate cashflow from sale (includes negative value)
    cf[i] <- max(loan - (1 - cost) * value, 0)

    return(cf)
}

# ---------------------------------------------------------------------------- #
# ----------------------------- Variable Annuity ----------------------------- #
# ---------------------------------------------------------------------------- #


#' Variable Annuity - Cashflow Simulator
#'
#' @param policy
#' Policy object containing necessary parameters (see create_policy_VA)
#' @param state
#' State vector containing state values for entire duration
#' @param data
#' Data frame containing all variables generated using other modules
#'
#' @return
#' Vector of cashflows for at each time point
cf_variable_annuity <- function(policy, state, data) {

    # Extract relevant policy variables
    value <- policy$value
    contract_length <- policy$length
    withdraw_prop <- policy$prop
    g_fee <- policy$g_fee

    # Initialize output vector
    cf <- rep(0, times = length(state))

    max_withdraw <- value * withdraw_prop
    account_value <- value

    i <- 1
    while (i <= contract_length) {

        # Compound account value - expenses for withdraw guarantee
        account_value <- account_value * (1 + data$stock[i]) * exp(-g_fee)

        if(state[i] == -1) {
            cf[i] = account_value
            break
        } else {
            if (i < contract_length) {
                cf[i] <- max_withdraw
                account_value <- max(account_value - max_withdraw, 0)
            } else {
                cf[i] <- account_value
            }
        }
        i <- i + 1
    }

    return(cf)
}
