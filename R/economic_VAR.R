esg_var_simulator = function (num_years = 5, num_paths = 1000, frequency = "quarter", perc_change = FALSE, return_sdf = TRUE, seed = NULL) {

    # ################
    # # error messages
    is.wholenumber = function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
    if (num_years <= 0 | num_paths <= 0 | !is.wholenumber(num_years)
        | !is.wholenumber(num_paths)) {
        stop("Number of years and paths to simulate must be positivie integers. ")

    } else if (!frequency %in% c("year", "quarter", "month")) {
        stop ("Frequency must be one of 'year', 'quarter', and 'month'. ")

    } else if (!is.logical(perc_change) | !is.logical(return_sdf)) {
        stop ("perc_change and return_sdf must be logical. ")

    }

    # ##########################################################
    # # VAR(2) calibrated coefficients (for stationary series) #
    # ##########################################################

    # variable names
    var_names = c("zcp3m_yield", "zcp10y_spread", "home_index", "GDP", "CPI", "ASX200", "AUD")
    sim_var_names = c(var_names, "mortgage_rate", "unemployment_rate")

    VAR = var_model()
    {
        intercept = VAR$intercept
        coef = as.matrix(VAR$coef)
        covres = as.matrix(VAR$covres)
        init_stat_2025q4 = as.numeric(VAR$init_stat_2025q4)
        init_orig = VAR$init_orig
        mortgage_rate = VAR$mortgage_rate
        unemployment_rate = VAR$unemployment_rate
        init_st = VAR$init_st
    }

    coef = as.matrix(coef)

    # ##################
    # # initialisation #
    # ##################

    init_qtr = as.Date("2026-1-1")
    num_pred = 4 * num_years
    time_index = seq(from = init_qtr, length.out = num_pred + 1, by = "quarter")
    path_index = paste("trajectory_", 1:num_paths, sep = "")
    progression = floor(num_paths / 5)

    # ############################
    # # step-by-step simulations #
    # ############################

    # Pre-generate all random normal shocks at once for speed
    set.seed(seed)
    all_noise = matrix(data = stats::rnorm(length(intercept) * num_pred * num_paths, 0, 1),
                       nrow = length(intercept))

    # Optimized single-path simulation engine using pre-allocated matrix multiplication
    var_path = function (num_pred, path_idx) {
        path = matrix(NA, nrow = num_pred, ncol = length(var_names))
        row.names(path) = as.character(time_index[-1])
        colnames(path) = var_names

        new_init = init_stat_2025q4
        chol_cov = t(chol(covres))

        # Extract the specific noise slice for this path
        col_start = (path_idx - 1) * num_pred + 1
        col_end = path_idx * num_pred
        path_noise = all_noise[, col_start:col_end, drop = FALSE]

        for (i in 1:num_pred) {
            e = path_noise[, i]
            zt = intercept + as.vector(coef %*% new_init) + as.matrix(chol_cov) %*% e
            path[i, ] = zt
            new_init = zt
        }
        return (as.data.frame(path))
    }

    # ########################################
    # # simulation for the stationary series #
    # ########################################

    prog_ind = 1; cat("Progress: 0% \n")
    var_sim_stationary = function (num_pred, num_paths) {
        v_path = vector("list", num_paths)
        for (x in 1:num_paths) {
            if (x == progression * prog_ind && prog_ind < 5) {
                cat(paste(20 * prog_ind, "%\n", sep = ""))
                prog_ind <<- prog_ind + 1
            }
            v_path[[x]] = var_path(num_pred, x)
        }
        return (v_path)
    }

    stat = var_sim_stationary(num_pred, num_paths)
    stat = lapply(stat, function (x) {cbind(x[,1:2], x[,3:7])})

    # ################################################
    # # convert forecast variables -> original units #
    # ################################################

    index2grow_inv = function (x, init) {
        Reduce (function (init, x) {init * exp(x)}, c(init, x), accumulate = TRUE)
    }

    sim = replicate(n = length(var_names),
                    expr = {data.frame(matrix(NA, nrow = num_pred, ncol = num_paths))},
                    simplify = F)
    sim = lapply(1:7, function (y) { lapply(1:num_paths, function (x) {stat[[x]][,y]}) })
    sim = lapply(sim, function (x) {as.data.frame(x)})
    sim[[1]] = rbind(init_orig[1], as.data.frame(sim[[1]])) # zcp3m
    sim[[2]] = rbind(init_orig[2], as.data.frame(sim[[2]])) # zcp10y_spread
    sim[[3]] = apply(sim[[3]], 2, function (x) {index2grow_inv(x, init_orig[3])}) # home_index
    sim[[4]] = apply(sim[[4]], 2, function (x) {index2grow_inv(x, init_orig[4])}) # GDP
    sim[[5]] = apply(sim[[5]], 2, function (x) {index2grow_inv(x, init_orig[5])}) # CPI
    sim[[6]] = apply(sim[[6]], 2, function (x) {index2grow_inv(x, init_orig[6])}) # ASX200
    sim[[7]] = apply(sim[[7]], 2, function (x) {index2grow_inv(x, init_orig[7])}) # AUD
    sim[[8]] = sim[[1]] + mortgage_rate # mortgage_rate
    sim[[9]] = sim[[2]] + unemployment_rate # unemployment_rate
    sim = lapply(sim, function(x){row.names(x) = time_index; colnames(x) = path_index; return (x)})
    names(sim) = sim_var_names

    # ###############################
    # # stochastic discount factors #
    # ###############################

    if (isTRUE(return_sdf)) {
        # Deterministic risk-free discounting calculation (fully vectorized loop)
        st = as.data.frame(matrix(NA, nrow = num_pred, ncol = num_paths))
        for (x in 1:num_paths) {
            st[, x] = exp(- (stat[[x]][, 1] / 4))
        }
        st = rbind(init_st, st)
        row.names(st) = as.character(time_index)
        colnames(st) = path_index

        sim[[length(sim_var_names) + 1]] = st
        names(sim)[length(sim_var_names) + 1] = "discount_factors"
    }

    # #################
    # # Adj frequency #
    # #################

    output = list()
    if (frequency == "month") {
        time_index_month = seq(from = init_qtr, length.out = num_years * 12 + 1, by = "month")

        qtr2month = function (x) {
            qtr_data = zoo::zoo (x, time_index)
            month_data = zoo::zoo (NA, time_index_month)
            data = merge (qtr_data, month_data)
            data$month_data = zoo::na.approx(data$qtr_data, rule=12)
            return (as.vector(data$month_data))
        }
        output = lapply(sim, function (x) apply(x, 2, qtr2month))
        output = lapply(output, function(x) {
            row.names(x) = as.character(time_index_month)
            return (x[-nrow(x), ])
        })

    } else if (frequency == "quarter") {
        output = lapply(sim, function(x) {x = x[-nrow(x), ]})

    } else if (frequency == "year") {
        time_index_year = seq(from = init_qtr, length.out = num_years, by = "year")
        output = lapply(sim, function (x) apply(x[-nrow(x), ], 2, function (y) {colMeans(matrix(y, nrow=4))} ))
        output = lapply(output, function(x) {row.names(x) = as.character(time_index_year); x})
    }

    output = lapply(output, function(x){x = t(as.data.frame(x))})

    if (isTRUE(return_sdf)) {
        if (frequency == "month") {
            output[["discount_factors"]] <- output[["discount_factors"]] ^ (1/3)
        } else if (frequency == "year") {
            sdf_qtr <- sim[["discount_factors"]][-nrow(sim[["discount_factors"]]), ]
            sdf_year <- apply(sdf_qtr, 2, function(y) {
                apply(matrix(y, nrow=4), 2, prod)
            })
            row.names(sdf_year) <- as.character(time_index_year)
            output[["discount_factors"]] <- t(as.data.frame(sdf_year))
        }
    }

    # #############
    # # Adj units #
    # #############

    if (isTRUE(perc_change)) {
        ref_level = lapply(output, function (x) {x = as.data.frame(x[,1]); colnames(x) = paste("ref_level", time_index[1]);x})

        output = lapply(names(output), function (name) {
            x <- output[[name]]
            if (name == "discount_factors") {
                return(x)
            } else {
                return((x[,-1] - x[,-ncol(x)]) / x[,-ncol(x)])
            }
        })

        output = lapply(1:length(output), function (x) {cbind(ref_level[[x]],output[[x]])})
        names(output) = names(sim)
    }

    cat("100% \n")
    return (output)
}
