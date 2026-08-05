# First open rit.Rprog

library("devtools")
load_all(export_all = FALSE)


#######################
# Example 1
la <- create_policy_LA(benefit = 1, defer = 0, increase = 0, frequency = "year")
cf_la <- simulate_cf(la, n = 10000, seed = 2026)
val_la <- value_policy(la, cf_la)


#######################
# Example 2

construct_state_matrix <- function(ages_at_death, max_years) {
    n_paths <- length(ages_at_death)

    # Create matrix with all entries = -1 (i.e. Policyholder is dead)
    state <- matrix(-1, nrow = n_paths, ncol = max_years)

    for (i in seq_len(n_paths)) {
        death_yr <- ages_at_death[i]

        # If they survived at least 1 full year, mark those completed years as 0 (Alive)
        if (death_yr > 1) {
            state[i, 1:(death_yr - 1)] <- 0
        }
    }
    return(state)
}

get_ages_at_death <- function(surv_probs, max_years, n_paths) {
    ages_at_death <- numeric(n_paths)
    for (i in seq(1, n_paths)) {
        p <- stats::runif(length(surv_probs))
        lv <- min(which((p > surv_probs) == TRUE))
        ages_at_death[i] <- lv
    }
    return(ages_at_death)
}

get_px_from_lx <- function(lx) {
    px <- numeric(length(lx))

    # Correctly calculate the 1-year survival prob starting from the very first year
    for (i in seq(1, length(lx) - 1)) {
        px[i] <- lx[i+1] / lx[i]
    }
    px[length(lx)] <- 0
    return(px)
}

library(lifecontingencies)
data(soaLt)

initial_age <- 65

# Initialize actuarial table
soa08Act <- with(soaLt, new("actuarialtable", interest=0.06,
                            x = x, lx = Ix, name = "SOA2008"))
lx08 <- soa08Act@lx[soa08Act@x >= initial_age]

# Extract survival probabilities from soa0xAct
surv <- get_px_from_lx(lx08)

# Simulate ages at death based on survival probabilities
N <- 10000
death_ages <- get_ages_at_death(surv, max_years = 100, n_paths = N)

# Construct state and SDF matrix
state <- construct_state_matrix(death_ages, max_years = 100)
sdf_ex2 <- list(sdf = t(matrix(rep((1+0.03)^-1, N*100), ncol  = N)))

# Policy valuation
la <- create_policy_LA(benefit = 1, defer = 0, increase = 0, frequency = "year")
cf_la <- simulate_cf(la, state = state, econ_var = sdf_ex2, n = N)
val_la <- value_policy(la, cf_la)

# Exact valuation
exact_val <- axn(actuarialtable = soa08Act, x = 65, payment = "arrears", i = 0.03)

# Check the difference
simulated_mean <- val_la$stats$mean
cat("--- Alignment Check ---\n")
cat("Exact Annuity-Immediate APV:      ", exact_val, "\n")
cat("Simulated Annuity-Immediate Mean: ", simulated_mean, "\n")
cat("Difference:                       ", abs(exact_val - simulated_mean), "\n")


#######################
# Example 3
sdf_ex3 <- list(sdf = t(matrix(rep((1+0.03)^-1, 1000*100), ncol  = 1000)))
la <- create_policy_LA(benefit = 100, defer = 0, increase = 0.01)
cf_la <- simulate_cf(la, econ_var = sdf_ex3)
val_la <- value_policy(la, cf_la)



#######################
# Example 4
calculatefee_VA <- function(prop, length, value, seed) {

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

fair_fee <- calculatefee_VA(prop = prop, length = length, value = value, seed = seed)





#######################
# Example 5
# 1. Create the Care Annuity Policy Object
# Benefits: 100 (Healthy), 150 (Disabled)
# Increase: 2% p.a., Guarantee: 5 years
ca_policy <- create_policy_CA(benefit = c(1000, 3000), increase = 0, min = 0, defer = 3)

# 2. Simulate Cashflows
# Note: When 'state' is not provided, simulate_cf() automatically uses the default Health State module parameters (the US HRS 3-state model)
cf_ca <- simulate_cf(ca_policy, init_age = 65, freq = 12, n = 1000)

# 3. Value the Policy
val_ca <- value_policy(ca_policy, cf_ca)


#######################
# Example 6
# 1. Create the Reverse Mortgage Policy Object
# Initial house value: $600,000; LVR = 100%; Transaction cost = 2.5%

calculatefee_RM <- function(value, LVR, trans_cost, seed) {

    state <- sim_indiv_path(init_age = 65, female = 1, seed = seed, n = N)
    econ_var <- esg_var_simulator(num_years = NCOL(state), num_paths = N, frequency = "year", seed = seed)
    econ_var$sdf <- econ_var$discount_factors
    econ_var$zcp3m <- econ_var$zcp3m_yield
    econ_var$house <- econ_var$home_index

    # Define the objective function
    f <- function(fee) {
        RM <- create_policy_RM(value = value, LVR = LVR, trans_cost = gamma, margin = fee)
        policy_cf <- simulate_cf(RM, n = N, state = state, econ_var = econ_var)
        L0 <- value * LVR
        Lt <- L0 * t(apply(1 + fee + econ_var$zcp3m, 1, cumprod))
        premium <- mean(rowSums((state + 1) * policy_cf$sdf * Lt * fee))

        # Turn off displaying figures in the function value_policy()
        pdf(file = tempfile())
        on.exit(dev.off(), add = TRUE)
        invisible(capture.output(policy_payoff <- value_policy(policy = RM, cashflows = policy_cf)))
        error <- policy_payoff[["stats"]][["mean"]] - premium
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
value <- 600000
LVR <- 0.64
gamma <- 0.01
seed <- 123
N <- 1e3

fair_margin <- calculatefee_RM(value = value, LVR = LVR, trans_cost = gamma, seed = seed)

rm <- create_policy_RM(value = value, LVR = LVR, trans_cost = gamma, margin = fair_margin)
cf_rm <- simulate_cf(rm)
val_rm <- value_policy(rm, cf_rm)


#######################
# Example 7
# Transform real-world mortality to risk-neutral measures
# 1. Load data and extract cohort mortality for a 65-year-old male: Use the AUS male rates from the package data
rates_P <- mortality_AUS_data$rate$male
qx_P_incomplete <- rate2rate(rates_P, from = "central", to = "prob")
old_ages <- 91:110
ages <- mortality_AUS_data$age
qx_P <- complete_old_age(qx_P_incomplete, ages, old_ages, method = "DG", type = "prob", closure_age = 110)
qx_P <- period2cohort(qx_P, 0:110, init_age = 65)
qx_P <- qx_P[,1] # The cohort born in 1970
# We use the package's cohort simulator to get the P-measure survival curve
surv_P <- sim_cohort_path_expected(init_age = 65, female = 0, death_probs = qx_P, closure_age = 110)
surv_P <- surv_P/surv_P[1]

# 2. Apply the Wang Transform to obtain Q-measure survival probabilities: The 'lambda' parameter represents the market price of longevity risk
surv_Q <- survivalP2Q(surv_P, method = "wang", lambda = 0.1)

# 3. Compare Expected Curtate Future Lifetime (Ex)
ex_P <- sum(surv_P) - 0.5 # Approximation for complete life expectancy
ex_Q <- sum(surv_Q) - 0.5
cat("Real-world Life Expectancy (P):", round(ex_P, 2), "\n") # 13.22
cat("Risk-neutral Life Expectancy (Q):", round(ex_Q, 2), "\n") # 14.01

# 4. Visualise the distortion: A higher survival curve in Q-measure implies a higher annuity price
plot(cumprod(surv_Q), type = "l", col = "red", lwd = 2, ylab = "Survival Probability", xlab = "Time (in years)")
lines(cumprod(surv_P), col = "blue", lwd = 2, lty = 2)
legend("topright", legend = c("Risk-Neutral (Q)", "Real-World (P)"), col = c("red", "blue"), lty = c(1, 2))



#######################
# Example 8 (Sherris and Wei (2021))
# 1. Create the Care Annuity Policy Object
# Benefits: 1000 (Healthy), 3000 (Disabled)
ca_policy <- create_policy_CA(benefit = c(1000, 0, 3000, 3000), increase = 0, min = 0, defer = 3, frequency = 'month')
sdf_ex8 <- list(sdf = t(matrix(rep((1+0.03)^-(1/12), 10000*12*100), ncol  = 10000)))
trans_probs <- get_trans_probs(n_states=5, model_type='S', param_file=US_HRS_5, init_age=65, female=1, year = 2022, latent = 0, frequency = 'month')
simulated_path <- simulate_health_state_paths(trans_probs, init_age=65, init_state = 0, cohort = 100)

# 2. Simulate Cashflows
# Note: When 'state' is not provided, simulate_cf() automatically uses the default Health State module parameters (the US HRS 3-state model)
cf_ca <- simulate_cf(ca_policy, init_age = 65, econ_var = sdf_ex8, n = 100, state = simulated_path)

# 3. Value the Policy
val_ca <- value_policy(ca_policy, cf_ca)


#############
trans_probs_5 <- get_trans_probs(n_states = 5, model_type = 'T', param_file = US_HRS_5, init_age = 87, female = 0, year = 2022, latent = 0, frequency = 'year')
lifetable_5 <- create_life_table(trans_probs_5, init_age = 87, init_state = 0, cohort = 100000)
head(lifetable_5,3)
simulated_path_5 <- simulate_health_state_paths(trans_probs_5, init_age = 87, init_state = 0, cohort = 10000)
prob_plots(init_age = 87, init_state = 0, trans_probs = trans_probs_5, frequency = 'year')
health_stats(n_states = 5, model_type = 'T', init_age = 87, init_state = 0, trans_probs = trans_probs_5, frequency = 'year')




