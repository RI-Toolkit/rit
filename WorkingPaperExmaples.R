# Open rit.Rprog

library("devtools")
load_all(export_all = FALSE)


#######################
# Example 1
la <- create_policy_LA(benefit = 100, defer = 0, increase = 0.01)
cf_la <- simulate_cf(la)
val_la <- value_policy(la, cf_la)


#######################
# Exmaple 2
construct_state_matrix <- function(ages_at_death, max_years) {
    n_paths <- length(ages_at_death)

    # Create matrix with all entries = -1 (i.e. PH is dead)
    state <- matrix(rep(-1, n_paths*max_years), nrow=n_paths,
                    ncol=max_years, byrow=TRUE)

    for (i in seq(1, n_paths)) {
        # Set all entries in i'th row prior to death as 0 (i.e. PH is alive)
        death_yr <- ages_at_death[i]
        state[i, 1:death_yr] <- rep(0, death_yr)
    }
    return(state)
}

get_ages_at_death <- function(surv_probs, max_years, n_paths) {
    ages_at_death <- c()
    for (i in seq(1, n_paths)) {
        p <- stats::runif(length(surv_probs))
        lv <- min(which((p > surv_probs) == TRUE))
        ages_at_death <- c(ages_at_death, lv)
    }
    return(ages_at_death)
}

get_px_from_lx <- function(lx) {
    px <- lx
    px[1] <- 1
    for (i in seq(2, length(lx) - 1)) {
        px[i] <- lx[i+1]/lx[i]
    }
    px[length(lx)] <- 0
    return(px)
}

library(lifecontingencies); data(soaLt);
initial_age <- 65
soa08Act <- with(soaLt, new("actuarialtable",interest=0.06,
                            x = x,lx = Ix,name = "SOA2008"))
lx08 <- soa08Act@lx[initial_age + 1:length(soa08Act@lx)]

# Extract survival probabilities from soa0xAct
surv <- get_px_from_lx(lx08)

# Simulate ages at death based on survival probabilities
death_ages <- get_ages_at_death(surv, max_years = 100, n_paths = 100)

# Construct state matrix
state <- construct_state_matrix(death_ages, max_years = 100)

la <- create_policy_LA(benefit = 100, defer = 0, increase = 0.01)
cf_la <- simulate_cf(la, n = 100, state = state)
val_la <- value_policy(la, cf_la)


#######################
# Example 3
sdf_ex3 <- list(sdf = t(matrix(rep((1+0.03)^-1, 100*100), ncol  = 100)))
la <- create_policy_LA(benefit = 100, defer = 0, increase = 0.01)
cf_la <- simulate_cf(la, n = 100, econ_var = sdf_ex3)
val_la <- value_policy(la, cf_la)



#######################
# Example 4

calculatefee <- function(prop, length, value, seed) {

    state <- sim_indiv_path(init_age = 65, female = 1, seed = seed, n = N)
    econ_var <- esg_var_simulator(num_years = NCOL(state), num_paths = N, frequency = "year", seed = seed)
    econ_var$sdf <- econ_var$discount_factors

    # Define the objective function
    f <- function(fee) {
        GMWB <- create_policy_VA(value = value, length = length, prop = prop, g_fee = fee)
        policy_cf <- simulate_cf(GMWB, n = N, state = state, econ_var = econ_var)

        # Turn off displaying figures in the function value_policy()
        pdf(file = tempfile())
        on.exit(dev.off(), add = TRUE)
        invisible(capture.output(policy_payoff <- value_policy(policy = GMWB, cashflows = policy_cf)))
        error <- policy_payoff[["stats"]][["mean"]] - value
        return(error)
    }

    # The Bisection root-finding method:
    # Check if the root is bracketed within the interval
    x0 <- 0
    x1 <- 0.2

    if (f(x0) * f(x1) > 0) {
        message("Enter valid interval !!!")
        return(NA) # Return NA (Not Available) to indicate failure
    } else {
        x2 <- (x0 + x1) / 2
        err <- abs(x0 - x1)

        # Initialize root to ensure it exists if the loop doesn't run (rare case)
        root <- x2

        while (err > 1e-4) {
            if (f(x0) * f(x2) < 0) {
                x1 <- x2
            } else {
                x0 <- x2
            }
            x2 <- (x0 + x1) / 2
            err <- abs(x2 - x1)
            root <- x2
        }
        return(root)
    }
}

# Product details
prop <- 0.02
length <- 15
value <- 100
seed <- 2026
N <- 1e3

fair_fee <- calculatefee(prop = prop, length = length, value = value, seed = seed)





#######################
# Example 5
# 1. Create the Care Annuity Policy Object
# Benefits: 100 (Healthy), 150 (Disabled)
# Increase: 2% p.a., Guarantee: 5 years
ca_policy <- create_policy_CA(benefit = c(100, 150), increase = 0.02, min = 5)

# 2. Simulate Cashflows
# Note: When 'state' is not provided, simulate_cf() automatically uses the default Health State module parameters (the US HRS 3-state model)
cf_ca <- simulate_cf(ca_policy, init_age = 65, n = 500)

# 3. Value the Policy
val_ca <- value_policy(ca_policy, cf_ca)
