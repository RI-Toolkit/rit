#' function to calculate transition rates at a certain age over a fractional year
#'
#' @param model_type
#' S for static model, T for trend model, F for frailty model
#' @param param_file
#' matrix of estimated parameters to construct the five state model. Generally,
#' use US_HRS_5 for 5 state model.
#' @param age
#' age of the individual
#' @param female
#' female 1 if female, 0 if male
#' @param wave_index
#' the wave index = (interview year - 1998)/2 + 1
#' @param latent
#' initial value of latent factor, normally take the value 0
#' @param frequency
#' string selecting the simulation frequency: "year", "quarter", or "month".
#' @return
#' 12 times 1 vector of transition rates for the 12 types of transitions
#'
#' @noRd
#'
health5_get_trans_rates = function(model_type, param_file, age, female, wave_index, latent, frequency){

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

    # Calculate fractional year step internally
    step <- 1 / freq

    # Subset parameters based on model type
    if (model_type == 'S') {
        param_file <- param_file[1:3, 3:14]
    } else if (model_type == 'T') {
        param_file <- param_file[6:9, 3:14]
    } else if (model_type == 'F') {
        param_file <- param_file[11:15, 3:14]
    }

    # Construct column vectors and advance by the fractional `step`
    if (model_type == 'S'){
        vari_x = matrix(c(1, age, female), ncol=1)
        vari_x1 = vari_x + cbind(c(0, step, 0))
    } else if (model_type == 'T'){
        vari_x = matrix(c(1, age, female, wave_index), ncol=1)
        vari_x1 = vari_x + cbind(c(0, step, 0, 0))
    } else if (model_type == 'F'){
        vari_x = matrix(c(1, age, female, wave_index, latent), ncol=1)
        vari_x1 = vari_x + cbind(c(0, step, 0, 0, 0))
    }

    # Matrix calculation to get a column vector of the ln transition rates
    ln_trans_rate_x = t(param_file) %*% vari_x
    ln_trans_rate_x1 = t(param_file) %*% vari_x1

    # Exact integration to get transition rates over the [age, age + step] interval
    # param_file[2,] is the age coefficient
    trans_rate = t(exp(ln_trans_rate_x1) / param_file[2,] - exp(ln_trans_rate_x) / param_file[2,])

    return(trans_rate)
}


#' function to get a matrix of transition probabilities at a certain age over a fractional year
#'
#' @param model_type
#' S for static model, T for trend model, F for frailty model
#' @param param_file
#' matrix of estimated parameters to construct the five state model.
#' @param age
#' age of the individual
#' @param female
#' female 1 if female, 0 if male
#' @param wave_index
#' the wave index = (interview year - 1998)/2 + 1
#' @param latent
#' initial value of latent factor, normally take the value 0
#' @param frequency
#' string selecting the simulation frequency: "year", "quarter", or "month".
#'
#' @return
#' 5 times 5 matrix of transitions probabilities, the states are H M D MD Dead
#' on the rows and columns
#' @import expm
#'
#' @noRd
#'
health5_get_trans_probs_at_age = function(model_type, param_file, age, female, wave_index, latent, frequency){

    # Pass freq down to the rates function
    trans_rate = health5_get_trans_rates(model_type, param_file, age, female, wave_index, latent, frequency)

    trans_rate_matrix = rbind(
        c(-sum(trans_rate[1:4]), trans_rate[1], trans_rate[2], trans_rate[3], trans_rate[4]),
        c(0, -sum(trans_rate[5:6]), 0, trans_rate[5], trans_rate[6]),
        c(trans_rate[7], trans_rate[8], -sum(trans_rate[c(7,8,9,10)]), trans_rate[9], trans_rate[10]),
        c(0, trans_rate[11], 0, -sum(trans_rate[11:12]), trans_rate[12]),
        c(0, 0, 0, 0, 0)
    )

    trans_prob_matrix = expm::expm(trans_rate_matrix)
    return(trans_prob_matrix)
}


#' the function to get a full list of transition probability matrices from the
#' initial age to closure age at a specified frequency
#'
#' @param model_type
#' S for static model, T for trend model, F for frailty model
#' @param param_file
#' matrix of estimated parameters to construct the five state model. Generally,
#' use US_HRS_5 for 5 state model.
#' @param init_age
#' the initial age of the transition probability matrices
#' @param closure_age
#' maximum life span
#' @param female
#' female 1 if female, 0 if male
#' @param wave_index
#' the wave index
#' @param latent
#' initial value of latent factor, normally take the value 0
#' @param frequency
#' string selecting the simulation frequency: "year", "quarter", or "month".
#'
#' @return a list of 5 times 5 transition probability matrices
#' @import readxl expm
#'
#' @noRd
#'
health5_get_trans_probs = function(model_type, param_file, init_age, closure_age, female, wave_index, latent, frequency){

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

    step <- 1 / freq
    ages <- seq(init_age, closure_age - step, by = step)
    n_steps <- length(ages)

    # Pre-allocate list of transition probability matrices for performance
    trans_prob_matrix = vector("list", n_steps)

    for (i in seq_len(n_steps)){
        a <- ages[i]

        # Wave index scales fractionally. The HRS wave index is scaled by 2 years in the data.
        current_wave_index <- wave_index + (a - init_age) / 2

        # calculate transition probability matrix for each fractional age step, passing freq down
        trans_prob_matrix[[i]] = health5_get_trans_probs_at_age(
            model_type, param_file, a, female, current_wave_index, latent, frequency
        )

        if (model_type == 'F'){
            # simulate the latent factor random walk, scaling the variance by the fractional step (1/freq)
            latent = latent + stats::rnorm(1, 0, sqrt(0.5 / freq))
        }
    }

    return(trans_prob_matrix)
}
