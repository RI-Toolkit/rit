# ---------------------------------------------------------------------------- #
# --------------------------- Account Based Pension -------------------------- #
# ---------------------------------------------------------------------------- #

#' Account Based Pension - Policy Object Generator
#'
#' Generates policy object for Account Based Pension policies
#'
#' @param balance
#' Initial balance of Pension account
#' @param init_age
#' Initial age of policyholder
#' @param drawdown_rate
#' Annual drawdown rate from initial age to 111 and the default value is the minimum
#' drawdown rate set by the Australian Taxation Office from age 65 to 111
#' @param frequency
#' string selecting the simulation frequency: "year", "quarter", or "month". Default is "year".
#' @return
#' Policy object
#' @export create_policy_AP
#' @examples
#' ap <- create_policy_AP(1000000, 50000)
create_policy_AP <- function(balance, init_age = 65, drawdown_rate = c(rep(0.05,10), rep(0.06,5), rep(0.07,5), rep(0.09,5), rep(0.11,5), rep(0.14,17)), frequency = "year") {

    if (!frequency %in% c("year", "quarter", "month")) stop("Frequency must be one of 'year', 'quarter', and 'month'.")
    if (length(balance) > 1)    stop("Invalid balance: length(balance) = 1")
    if (balance < 0)            stop("Invalid balance: balance > 0")
    if (length(drawdown_rate) != 111 - init_age + 1)  stop("Invalid drawdown rate: length of drawdown_rate should be 111 - init_age + 1")

    minimum_drawdown_rate <- c(rep(0.05,10), rep(0.06,5), rep(0.07,5), rep(0.09,5), rep(0.11,5), rep(0.14,17))
    ind <- (init_age - 65 + 1):length(minimum_drawdown_rate)

    if(any(drawdown_rate < minimum_drawdown_rate[ind])) stop("Invalid drawdown rate: value should be greater than the minimum drawdown rate")

    pol <- data.frame(name = c("AP"),
                      bal = c(balance),
                      age = c(init_age:111),
                      rate = c(drawdown_rate),
                      frequency = c(frequency))

    return(pol)
}

# ---------------------------------------------------------------------------- #
# ------------------------------- Care Annuity ------------------------------- #
# ---------------------------------------------------------------------------- #

#' Care Annuity - Policy Object Generator
#'
#' Generates policy object for Care Annuity policies
#'
#' @param benefit
#' Annual benefit for each policy (vector)
#' @param increase
#' Annual rate of increase for each policy (vector)
#' @param min
#' Minimum guaranteed period for each policy (vector)
#' @param defer
#' Deferment period of policy
#' @param frequency
#' string selecting the simulation frequency: "year", "quarter", or "month". Default is "year".
#' @return
#' Policy object
#' @export create_policy_CA
#' @examples
#' ca <- create_policy_CA(c(60000, 1200), 0.04, 5, 0)
create_policy_CA <- function(benefit, increase, min, defer, frequency = "year") {

    if (!frequency %in% c("year", "quarter", "month")) stop("Frequency must be one of 'year', 'quarter', and 'month'.")
    if (min < 0)       stop("Invalid min: min > 0")
    if (sum(benefit < 0))   stop("Invalid benefit: benefit[i] > 0")
    if (increase < 0)  stop("Invalid increase: increase > 0")
    if (defer < 0)     stop("Invalid defer: defer > 0")

    n_policies = length(benefit)

    if (length(benefit) == 2) {
        state <- c("H", "S")
    } else if (length(benefit) == 4) {
        state <- c("H", "M", "D", "MD")
    } else {
        stop("Invalid input: vector of benefit must be of length 2 or 4")
    }

    pol <- data.frame(name = c("CA"),
                      state = c(state),
                      increase = c(increase),
                      benefit = c(benefit),
                      min = c(min),
                      defer = c(defer),
                      frequency = c(frequency))

    return(pol)
}


# ---------------------------------------------------------------------------- #
# ------------------------------- Life Annuity ------------------------------- #
# ---------------------------------------------------------------------------- #

#' Life Annuity - Policy Object Generator
#'
#' Generates policy object for Life Annuity  policies
#'
#' @param benefit
#' Annual benefit for policy
#' @param defer
#' (Optional) Deferment period of policy
#' @param increase
#' (Optional) Annual rate of increase for policy
#' @param frequency
#' string selecting the simulation frequency: "year", "quarter", or "month". Default is "year".
#' @return
#' Policy object
#' @export create_policy_LA
#' @examples
#' la <- create_policy_LA(60000, 5, 0)
create_policy_LA <- function(benefit, defer = 0, increase = 0, frequency = "year") {

    if (!frequency %in% c("year", "quarter", "month")) stop("Frequency must be one of 'year', 'quarter', and 'month'.")
    if (length(defer) > 1)      stop("Invalid defer: length(defer) = 1")
    if (length(benefit) > 1)    stop("Invalid benefit: length(benefit) = 1")
    if (length(increase) > 1)   stop("Invalid increase: length(increase) = 1")

    if (defer < 0)              stop("Invalid defer: defer > 0")
    if (benefit < 0)            stop("Invalid benefit: benefit > 0")
    if (increase < 0)           stop("Invalid increase: increase > 0")

    pol <- data.frame(name = c("LA"),
                      benefit = c(benefit),
                      defer = c(defer),
                      increase = c(increase),
                      frequency = c(frequency))

    return(pol)
}

# ---------------------------------------------------------------------------- #
# ------------------------------ Pooled Annuity ------------------------------ #
# ---------------------------------------------------------------------------- #

#' Pooled Annuity - Policy Object Generator
#'
#' Generates policy object for Pooled Annuity policies
#'
#' @param benefit
#' Initial annual benefit for policy
#' @param size
#' Size of pool
#' @param interest
#' Assumed investment earnings rate (flat)
#' @param frequency
#' string selecting the simulation frequency: "year", "quarter", or "month". Default is "year".
#' @return
#' Policy object
#' @export create_policy_PA
#' @examples
#' pa <- create_policy_PA(60000, 10000, 0.05)
create_policy_PA <- function(benefit, size, interest, frequency = "year") {

    if (!frequency %in% c("year", "quarter", "month")) stop("Frequency must be one of 'year', 'quarter', and 'month'.")
    if (length(size) > 1)       stop("Invalid size: length(size) = 1")
    if (length(benefit) > 1)    stop("Invalid benefit: length(benefit) = 1")
    if (length(interest) > 1)   stop("Invalid interest: length(interest) = 1")

    if (size < 0)               stop("Invalid size: size > 0")
    if (benefit < 0)            stop("Invalid benefit: benefit > 0")
    if (interest < 0)           stop("Invalid interest: interest > 0")

    if (!all.equal(size, as.integer(size)))
        stop("Invalid size: size needs to be integer")

    pol <- data.frame(name = c("PA"),
                      benefit = c(benefit),
                      size = c(size),
                      interest = c(interest),
                      frequency = c(frequency))

    return(pol)
}

# ---------------------------------------------------------------------------- #
# ----------------------------- Reverse Mortgage ----------------------------- #
# ---------------------------------------------------------------------------- #

#' Reverse Mortgage - Policy Object Generator
#'
#' Generates policy object for Reverse Mortgage policies
#'
#' @param value
#' Initial value of property
#' @param margin
#' Lending margin
#' @param LVR
#' Loan to Value ratio for PH
#' @param trans_cost
#' Transaction cost associated with sale
#' @param frequency
#' string selecting the simulation frequency: "year", "quarter", or "month". Default is "year".
#' @return
#' Policy object
#' @export create_policy_RM
#' @examples
#' rm <- create_policy_RM(600000, 0.64, 0.01, 0.04)
create_policy_RM <- function(value, LVR, trans_cost, margin, frequency = "year") {

    if (!frequency %in% c("year", "quarter", "month")) stop("Frequency must be one of 'year', 'quarter', and 'month'.")
    if (length(LVR) > 1)        stop("Invalid LVR: length(LVR) = 1")
    if (length(value) > 1)      stop("Invalid value: length(value) = 1")
    if (length(margin) > 1)     stop("Invalid margin: length(margin) = 1")
    if (length(trans_cost) > 1) stop("Invalid trans_cost: length(t_cost) = 1")

    if (value < 0)              stop("Invalid value: value > 0")
    if (margin < 0)             stop("Invalid margin: margin > 0")
    if (LVR < 0 | LVR > 1)      stop("Invalid LVR: 0 <= LVR <= 1")
    if (trans_cost < 0)         stop("Invalid trans_cost: 0 < trans_cost < 1")
    if (trans_cost > 1)         stop("Invalid trans_cost: 0 < trans_cost < 1")

    pol <- data.frame(name = c("RM"),
                      value = c(value),
                      LVR = c(LVR),
                      trans_cost = c(trans_cost),
                      margin = c(margin),
                      frequency = c(frequency))

    return(pol)
}

# ---------------------------------------------------------------------------- #
# ----------------------------- Variable Annuity ----------------------------- #
# ---------------------------------------------------------------------------- #

#' Variable Annuity - Policy Object Generator
#'
#' Generates policy object for Variable Annuity policies
#'
#' @param value
#' Initial Account value (single upfront premium)
#' @param length
#' Length of contract (in years)
#' @param prop
#' Proportion of the premium that can be withdrawn annually
#' @param g_fee
#' Continuous guarantee fees for GMWB
#' @param frequency
#' string selecting the simulation frequency: "year", "quarter", or "month". Default is "year".
#'
#' @return
#' Policy object
#' @export create_policy_VA
#' @examples
#' va <- create_policy_VA(1000000, 30, 0.1, 0.01)
create_policy_VA <- function(value, length, prop, g_fee, frequency = "year") {

    if (!frequency %in% c("year", "quarter", "month")) stop("Frequency must be one of 'year', 'quarter', and 'month'.")
    if (length(prop) > 1)       stop("Invalid prop: length(prop) = 1")
    if (length(g_fee) > 1)      stop("Invalid g_fee: length(g_fee) = 1")
    if (length(length) > 1)     stop("Invalid length: length(length) = 1")
    if (length(prop) > 1)       stop("Invalid prop: length(prop) = 1")

    if (g_fee < 0)              stop("Invalid g_fee: g_fee > 0")
    if (value < 0)              stop("Invalid value: value > 0")
    if (length < 0)             stop("Invalid length: length > 0")
    if (prop < 0 | prop > 1)    stop("Invalid prop: 0 <= prop <= 1")

    pol <- data.frame(name = c("VA"),
                      value = c(value),
                      length = c(length),
                      prop = c(prop),
                      g_fee = c(g_fee),
                      frequency = c(frequency))

    return(pol)
}
