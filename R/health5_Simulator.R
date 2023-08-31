#' function to get the path of set number of individuals
#'
#' @param list_trans_probs
#' a list of transition probability matrices, preferably generated from \code{get_trans_probs}.
#' @param init_age
#' the initial age of the path
#' @param init_state
#' the initial state of all individuals
#' @param cohort
#' the number of simulated cohorts
#'
#' @return
#' rows of individual paths in the states, 0 for H, 1 for M, 2 for D, 3 for MD, -1 for Dead
#' for each row it starts from the initial age as an input, and end at age 110
#'
#' @examples example
health5_simulate_paths <- function(list_trans_probs, init_age, init_state, cohort) {
    # init_state 0 for H, 1 for M, 2 for D, 3 for MD, -1 for Dead

    # create empty matrix to contain simulated population
    simulated_pop <- matrix(0, nrow = cohort, ncol = 110-init_age+2)

    # initialise all individuals
    simulated_pop[, 1] <- init_state

    for (i in 2:ncol(simulated_pop)) {
        simulated_pop[simulated_pop[,i-1] == 0, i] <- sample(c(0, 1, 2, 3, -1),
                                                             sum(simulated_pop[, i-1] == 0),
                                                             replace = TRUE,
                                                             prob = list_trans_probs[[i-1]][1, ])


        simulated_pop[simulated_pop[,i-1] == 1, i] <- sample(c(0, 1, 2, 3, -1),
                                                             sum(simulated_pop[, i-1] == 1),
                                                             replace = TRUE,
                                                             prob = list_trans_probs[[i-1]][2, ])

        simulated_pop[simulated_pop[,i-1] == 2, i] <- sample(c(0, 1, 2, 3, -1),
                                                             sum(simulated_pop[, i-1] == 2),
                                                             replace = TRUE,
                                                             prob = list_trans_probs[[i-1]][3, ])


        simulated_pop[simulated_pop[,i-1] == 3, i] <- sample(c(0, 1, 2, 3, -1),
                                                             sum(simulated_pop[, i-1] == 3),
                                                             replace = TRUE,
                                                             prob = list_trans_probs[[i-1]][4, ])

        simulated_pop[simulated_pop[,i-1] == -1, i] <- -1
    }
    return(simulated_pop)
}

#' the function to create lifetable for Static or Trend models
#'
#' @param list_trans_probs
#' a list of transition probability matrices, preferably generated from \code{{get_trans_probs}}.
#' @param init_age
#' the initial age of the path
#' @param init_state
#' the initial state of all individuals
#' @param cohort
#' the number of simulated cohorts
#'
#' @return
#' a life table matrix
#' the row represents the age from the input initial age to 110, and the columns are states H M D MD Dead
#'
#'
#' @examples example
health5_create_life_table=function(list_trans_probs,init_age,init_state,cohort){
        # list of 46 matrices of transition probabilities for this simulation
        #list of lifetables
        state_status=matrix(nrow = 110-init_age+1, ncol = 20)
        colnames(state_status) <- c("Age", "Alive", "H", "M", "D", "MD", "Dead","H_M","H_D","H_MD","H_Dead","M_MD","M_Dead","D_H","D_M","D_MD","D_Dead","MD_M","MD_Dead","H.M.D.MD_Dead")
        if (init_state==0){
            # initial state status is 1 in the healthy state and 0 for the others
            state_status[1,]=c(init_age,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
        }
        if (init_state==1){
            # initial state status is 1 in the M state and 0 for the others
            state_status[1,]=c(init_age,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
        }
        if (init_state==2){
            # initial state status is 1 in the D state and 0 for the others
            state_status[1,]=c(init_age,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
        }
        if (init_state==3){
            # initial state status is 1 in the MD state and 0 for the others
            state_status[1,]=c(init_age,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
        }

        for (age in init_age:(110-1)){
            for (j in 3:7){
                state_status[age-init_age+2,j]=state_status[age-init_age+1,3]*list_trans_probs[[age-init_age+1]][1,j-2]+
                    state_status[age-init_age+1,4]*list_trans_probs[[age-init_age+1]][2,j-2]+
                    state_status[age-init_age+1,5]*list_trans_probs[[age-init_age+1]][3,j-2]+
                    state_status[age-init_age+1,6]*list_trans_probs[[age-init_age+1]][4,j-2]+
                    state_status[age-init_age+1,7]*list_trans_probs[[age-init_age+1]][5,j-2]
                # state_status is essentially the life table at each age
                # firstly times the probability/or number(if assume 10000 individuals just simply times 10000)
                # of the state at the beginning of this year with the probability of entering state j from this state
                # then take the sum of all possible initial states of this year
            }
            state_status[age-init_age+2,'Alive']=1-state_status[age-init_age+2,'Dead']
            state_status[age-init_age+2,'Age']=age+1

            state_status[age-init_age+2,'H_M']=state_status[age-init_age+1,'H']*list_trans_probs[[age-init_age+1]][1,2]
            state_status[age-init_age+2,"H_D"]=state_status[age-init_age+1,'H']*list_trans_probs[[age-init_age+1]][1,3]
            state_status[age-init_age+2,"H_MD"]=state_status[age-init_age+1,'H']*list_trans_probs[[age-init_age+1]][1,4]
            state_status[age-init_age+2,"H_Dead"]=state_status[age-init_age+1,'H']*list_trans_probs[[age-init_age+1]][1,5]
            state_status[age-init_age+2,"M_MD"]=state_status[age-init_age+1,'M']*list_trans_probs[[age-init_age+1]][2,4]
            state_status[age-init_age+2,"M_Dead"]=state_status[age-init_age+1,'M']*list_trans_probs[[age-init_age+1]][2,5]
            state_status[age-init_age+2,"D_H"]=state_status[age-init_age+1,'D']*list_trans_probs[[age-init_age+1]][3,1]
            state_status[age-init_age+2,"D_M"]=state_status[age-init_age+1,'D']*list_trans_probs[[age-init_age+1]][3,2]
            state_status[age-init_age+2,"D_MD"]=state_status[age-init_age+1,'D']*list_trans_probs[[age-init_age+1]][3,4]
            state_status[age-init_age+2,"D_Dead"]=state_status[age-init_age+1,'D']*list_trans_probs[[age-init_age+1]][3,5]
            state_status[age-init_age+2,"MD_M"]=state_status[age-init_age+1,'MD']*list_trans_probs[[age-init_age+1]][4,2]
            state_status[age-init_age+2,"MD_Dead"]=state_status[age-init_age+1,'MD']*list_trans_probs[[age-init_age+1]][4,5]
            state_status[age-init_age+2,"H.M.D.MD_Dead"]=state_status[age-init_age+2,"H_Dead"]+state_status[age-init_age+2,"M_Dead"]+state_status[age-init_age+2,"D_Dead"]+state_status[age-init_age+2,"MD_Dead"]

            expected_time_state=colSums(state_status) # the order is H M D MD Dead
        }
        # the size of the lists will be large
        state_status[,c(2:20)]=state_status[,c(2:20)]*cohort
        state_status=data.frame(state_status)
    return(state_status)
}

#' the function to get n_sim number of simulated lifetables for Frailty model
#'
#' @param model_type
#' choose F for Frailty model
#' @param param_file
#' matrix of estimated parameters to construct the five state model. Generally,
#' use US_HRS_5 for 5 state model.
#' @param female
#' female 1 if female, 0 if male
#' @param wave_index
#' the wave index = (interview year - 1998)/2 + 1
#' @param latent
#' initial value of latent factor, normally take the value 0
#' @param init_age
#' the initial age of the life table
#' @param init_state
#' 0 for H state, 1 for M state, 2 for D state, 3 for MD state
#' @param n_sim
#' the number of simulations
#' @param cohort
#' number of people at the beginning of the life table
#' @param mean
#' TRUE to return expected life table, FALSE to return all simulated life tables
#'
#' @return
#' a list of n_sim number of life table matrices when mean=FALSE
#' for each matrix, the row represents the age from the input initial age to 110, and the columns are states H M D MD Dead
#' or the mean life table when mean=TRUE
#'
#' @examples example
health5_simulate_life_table=function(model_type, param_file, female, wave_index,latent,init_age,init_state,n_sim, cohort, mean){
    if (model_type != 'F') {
        stop('use frailty model to simulate lifetables')
    }
    state_status_full=list()
    for (i in 1: n_sim){
    list_trans_probs=health5_get_trans_probs(model_type, param_file, init_age, female, wave_index, latent)
    state_status=health5_create_life_table(list_trans_probs,init_age,init_state, cohort)
    state_status_full=append(state_status_full,list(state_status))
    }
    if (mean != TRUE){
        return(state_status_full)
    }
    if(mean == TRUE){
        return(Reduce('+', state_status_full)/n_sim)
    }
}

