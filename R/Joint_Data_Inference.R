# joint Inference
incidenceMerge = function(Traj, gridsize){
  n = dim(Traj)[1]
  k = floor(n / gridsize)
  Traj2 = matrix(ncol = 2, nrow = k)
  Traj2[,1] = Traj[seq(gridsize,n,by = gridsize),1]
  for(i in 1:k){
    Traj2[i,2] = sum(Traj[((i-1)*gridsize+1) : (i * gridsize),2])
  }
  return(Traj2)
}


incid_Sampling = function(Traj, rho, phi){
  res = apply(Traj, 1, function(x){
    return(c(x[1], rnbinom(1, size = phi, mu = rho * x[2])))
  })
  return(t(res))
}


MCMC_setup_Joint_Inference = function(tree,raw_incidence_obs, times,t_correct,N,gridsize=50,niter = 1000,burn = 500,thin = 5, changetime,DEMS = c("E", "I"),
         prior=list(pop_pr=c(1,10), R0_pr=c(1,7), mu_pr = c(3,0.2), gamma_pr = c(3,0.2), ch_pr = 1,hyper_pr=c(0.01,0.01)),
         proposal = list(pop_prop = 1, R0_prop = c(0.01), mu_prop=0.1, gamma_prop = 0.2, ch_prop=0.05),
         control = list(), likelihood = "volz", model = "SIR", Index = c(0,1), nparam = 2,PCOV = NULL,REAL = F){

  Times = min(times)
  for(i in 2:length(times)){
    Times = c(Times,seq(times[i-1],times[i],length.out = gridsize+1)[-1])
  }
  gridset = seq(1,length(Times),by=gridsize)
  grid = times
  cut_times = max(grid)
  coal_obs = tree
  if(REAL){
    newtree = coal_lik_init_tree(tree, grid, t_correct[1], t_correct[2])
    cut_times = max(grid) - as.numeric(as.Date(t_correct[1]) - as.Date(t_correct[2]))/365
    tree = newtree[[1]]
    t_correct = newtree[[2]]
    coal_obs = summarize_phylo(tree)
  }else{
    cut_times = max(grid) - t_correct
  }
  c = 0
  for(i in 1:length(changetime)){
    if(changetime[i] <= cut_times){
        c = c + 1
    }else{
      break;
    }
  }
  Init = coal_lik_init(coal_obs$samp_times, coal_obs$n_sampled, coal_obs$coal_times, grid,t_correct)
  grid = grid[1:(sum(grid < cut_times) + 1)]
  mgrid = length(grid)
  n = dim(raw_incidence_obs)[1]
  i = 1
  k = 2

  l = numeric(n)
  l2 = numeric(n)
  l3 = numeric(n)
  l2[1] =  sum((raw_incidence_obs[1,1] >= grid) & (2 * raw_incidence_obs[1,1] - raw_incidence_obs[2,1] < grid))
  for(j in 2:n){
     l2[j] = sum((raw_incidence_obs[j,1] >= grid) & (raw_incidence_obs[j-1,1] < grid))
     l3[j] = sum((raw_incidence_obs[j,1] > grid) & (raw_incidence_obs[j-1,1] <= grid))
  }

  if(min(grid)<=2*raw_incidence_obs[1,1] - raw_incidence_obs[2,1]){
    Incidence_start = which.max(l2>0)
  }else{
    Incidence_start = which.max(l3>0)
  }
  Incidence_end = n - which.max(rev(l2>0)) + 1
  j = Incidence_start
  time_inc_start = 0
 # if(j == 1){
#    time_inc_start = which.max((raw_incidence_obs[1,1] >= grid) & (2 * raw_incidence_obs[1,1] - raw_incidence_obs[2,1] < grid))
#  }else{
   # time_inc_start = which.max((raw_incidence_obs[j,1] >= grid) & (raw_incidence_obs[j-1,1] < grid))
    time_inc_start = which(abs(grid - raw_incidence_obs[Incidence_start,1] + raw_incidence_obs[2,1] - raw_incidence_obs[1,1])<0.0000000001)
 # }

  #for(i in 1:(length(grid)-1)){
  #  l = l + ((raw_incidence_obs[,1] > grid[i]) & (raw_incidence_obs[,1] <= grid[i+1]))
    #if(time_inc_start == 0 && sum(l)>0)
#  }
  #mask = (l > 0)
  #l = l[mask]
  #Incidence_start = which.max(mask)
  #Incidence_end = n - which.max(rev(mask)) + 1
  ###
  # test

  ###
  time_inc_end = which(abs(grid - raw_incidence_obs[Incidence_end,1])<0.00000001)
 # if(raw_incidence_obs[Incidence_end,1]> times[time_inc_end]){
#    stop("error with the grid")
#  }
  Incidence_obs = raw_incidence_obs[Incidence_start:Incidence_end,]
  IncidenceList = list(Incidence_obs = Incidence_obs, gridrep = time_inc_start:time_inc_end - 1, idmx = time_inc_end - 1)
  if(is.null(t_correct)){
    t_correct = max(coal_obs$coal)
  }
  x_i = c(length(changetime),nparam,Index[1],Index[2], c,mgrid)
  print("time grids")
  print(times[1:max(time_inc_end, Init$ng + 1)])
  print("_________________")
  MCMC_setting = list(Init = Init, Incidence = IncidenceList, times = Times,t_correct = t_correct,x_r = c(N,changetime),
                      gridsize = gridsize,gridset = gridset, niter = niter,burn = burn,thin = thin,x_i = x_i,
                      prior = prior, proposal = proposal, control = control, p = length(DEMS),
                      reps = 1, likelihood = likelihood, model = model,PCOV = PCOV)
  cat("MCMC set up ready \n")

  return(MCMC_setting)
}

#####

MCMC_initialize_Integrate = function(MCMC_setting, enable = c(T,T)){

  logMultiNorm = NaN
  coalLog = NaN
  IncidLog = NaN
  ########
  par_probs = numeric(5)
  while(is.nan(logMultiNorm)||is.nan(coalLog) || is.nan(IncidLog)){
    state = numeric(MCMC_setting$p)
    state[1] = MCMC_setting$x_r[1]
    for(i in 2:MCMC_setting$p){

      if(is.null(MCMC_setting$control$alpha)){
        state[i] = rlnorm(1, MCMC_setting$prior$pop_pr[2*i-3],MCMC_setting$prior$pop_pr[2 * i-2])
      }else{
        state[i] = MCMC_setting$control$alpha[i-1] * MCMC_setting$x_r[1]
      }
      par_probs[i-1] = dnorm(log(state[i]), MCMC_setting$prior$pop_pr[2 * i - 3],MCMC_setting$prior$pop_pr[2 * i - 2],log = T)
    }

    if(is.null(MCMC_setting$control$R0)){
      R0 = runif(1,MCMC_setting$prior$R0_pr[1],MCMC_setting$prior$R0_pr[2])
    }else{
      R0 = MCMC_setting$control$R0
    }

    par_probs[2] = dunif(R0, MCMC_setting$prior$R0_pr[1], MCMC_setting$prior$R0_pr[2], log = T)

    if(MCMC_setting$model == "SEIR" || MCMC_setting$model == "SEIR2"){

      if(is.null(MCMC_setting$control$mu)){
        mu = exp(rnorm(1,MCMC_setting$prior$mu_pr[1],MCMC_setting$prior$mu_pr[2]))
      }else{
        mu = MCMC_setting$control$mu
      }
      par_probs[4] = dnorm(log(mu),MCMC_setting$prior$mu_pr[1],MCMC_setting$prior$mu_pr[2],log = T)
    }

    if(is.null(MCMC_setting$control$gamma)){
      gamma = exp(rnorm(1,MCMC_setting$prior$gamma_pr[1], MCMC_setting$prior$gamma_pr[2]))
    }else{
      gamma = MCMC_setting$control$gamma
    }

    par_probs[3] = dnorm(log(gamma),MCMC_setting$prior$gamma_pr[1],MCMC_setting$prior$gamma_pr[2],log = T)

    ch=c()

    if(is.null(MCMC_setting$control$hyper)){
      hyper = rgamma(1, MCMC_setting$prior$hyper_pr[1],MCMC_setting$prior$hyper_pr[2])
    }else{
      hyper = MCMC_setting$control$hyper
    }
    incidence_par = NULL
    if(!is.null(MCMC_setting$control$incidence_par)){
      incidence_par = MCMC_setting$control$incidence_par
    }
    par_probs[5] = dgamma(hyper, MCMC_setting$prior$hyper_pr[1],MCMC_setting$prior$hyper_pr[2], log = T)

    if(length(MCMC_setting$x_r) > 1){
      if(is.null(MCMC_setting$control$ch)){
        ch = rlnorm(MCMC_setting$x_i[1],0,1/hyper)
      }else{
        ch = MCMC_setting$control$ch
      }
    }
    ##########
    if(MCMC_setting$x_i[5] < MCMC_setting$x_i[1]){
      ch[(MCMC_setting$x_i[5] + 1):MCMC_setting$x_i[1]] = 1
    }
    ########
    if(MCMC_setting$model == "SEIR" || MCMC_setting$model == "SEIR2"){
      param = c(R0,mu,gamma,ch,hyper)
    }else if(MCMC_setting$model == "SIR"){
      param = c(R0, gamma,ch,hyper)
    }

    paramlist = New_Param_List(param, state, MCMC_setting$gridsize, MCMC_setting$times, MCMC_setting$x_r,
                               MCMC_setting$x_i, model = MCMC_setting$model)

    # Ode_Traj_thin = ODE_rk45(state, MCMC_setting$times, param, MCMC_setting$x_r, MCMC_setting$x_i,model = MCMC_setting$model)

    Ode_Traj_coarse = paramlist$Ode
    FT = paramlist$FT
    betaN = paramlist$betaN
    # FT = KF_param_chol(Ode_Traj_thin,param,MCMC_setting$gridsize, MCMC_setting$x_r,MCMC_setting$x_i,model = MCMC_setting$model)

    #betaN = beta
    if(is.null(MCMC_setting$control$traj)){
      Latent = Traj_sim_general_noncentral(Ode_Traj_coarse,FT,MCMC_setting$t_correct)
      #LatentTraj = Latent$SimuTraj
      LatentTraj = Latent$SimuTraj
      OriginTraj = Latent$OriginTraj
      logMultiNorm = Latent$logMultiNorm
      logOrigin = Latent$logOrigin
    }else{
      OriginTraj = MCMC_setting$control$traj
      LatentTraj = TransformTraj(Ode_Traj_coarse, OriginTraj,FT)
      logMultiNorm = log_like_traj_general_adjust(OriginTraj,Ode_Traj_coarse,FT,MCMC_setting$gridsize, MCMC_setting$t_correct)
      logOrigin = - sum(OriginTraj * OriginTraj) / 2
      if( sum(abs(LatentTraj[1,2:(MCMC_setting$p+1)] - state)) > 1){
        print("not consistent")
      }
      logMultiNorm = log_like_traj_general_adjust(LatentTraj,Ode_Traj_coarse,
                                                  FT,MCMC_setting$gridsize,MCMC_setting$t_correct)
    }

    #coalLog = Structural_Coal_lik(MCMC_setting$Init_Detail,LatentTraj,param, MCMC_setting$x_r, MCMC_setting$x_i,
    #                             model = "SEIR2")
    coalLog = volz_loglik_nh_adj(MCMC_setting$Init,LatentTraj,
                              betaN, MCMC_setting$t_correct, MCMC_setting$x_i[3:4],enable = enable[1])
    #IncidLog = log_incidence(LatentTraj, MCMC_setting$Incidence$Incidence_obs, Incid_Par = incidence_par)
    IncidLog = log_incidence(LatentTraj, MCMC_setting$Incidence, Incid_Par = incidence_par, Pois = F,enable = enable[2])
    print(paste("coalescent likelihood after initialization ", coalLog))

    if(!is.nan((coalLog))){
      coalLog = ifelse(coalLog <= - 10000000,NaN, coalLog)
    }
    if(!is.nan((IncidLog))){
      IncidLog = ifelse(IncidLog <= - 10000000,NaN, IncidLog)
    }
    plot(LatentTraj[,1],LatentTraj[,3],type="l",xlab = "time",col="blue", lwd = 2)
    #lines(LatentTraj[,1],LatentTraj[,3],col="red", lwd = 2)
  }

  chprobs = -0.5 * sum((log(param[(MCMC_setting$x_i[2] + 1):(MCMC_setting$x_i[1] + MCMC_setting$x_i[2])]) * hyper)^2)
  #LogGamma = dlnorm(gamma,MCMC_setting$prior$gamma_pr[1],MCMC_setting$prior$gamma_pr[2],log = T)
  #LogE = dlnorm(mu, MCMC_setting$prior$mu_pr[1], MCMC_setting$prior$mu_pr[2],log = T)
  plot(LatentTraj[,1],LatentTraj[,3],type="l",xlab = "time",col="blue", lwd = 2)

  paras =  c(state,param)
  MCMC_obj = list(par = paras,incid_par = incidence_par,LatentTraj = LatentTraj,logMultiNorm = logMultiNorm,p = MCMC_setting$p,
                  Ode_Traj_coarse = Ode_Traj_coarse, FT = FT, coalLog = coalLog, IncidLog = IncidLog, OriginTraj = OriginTraj,logOrigin = logOrigin,
                  par_probs = par_probs, chprobs = chprobs, betaN = betaN)
  ##########
  cat("Initialize MCMC \n")
  print(paste("population size = ", MCMC_setting$N))
  print(paste(paras))
  return(MCMC_obj)
}



sigmoid = function(x){
  return(1/(1+exp(-x)))
}

sigmoid_inv = function(y){
  return(log(y) - log(1-y))
}





Update_incidence_Pars = function(MCMC_setting, MCMC_obj,update=c(1,1),joint = F){
  npar = length(MCMC_obj$incid_par)

  if(joint){
    incid_par_new = MCMC_obj$incid_par
    t_old = sigmoid_inv(MCMC_obj$incid_par[1])
    t_new = t_old + rnorm(1,0,MCMC_setting$proposal$incid_par_prop[1])
    incid_par_new[1] = sigmoid(t_new)
    t_old2 = log(MCMC_obj$incid_par[2])
    t_new2 = t_old2 + rnorm(1,0,MCMC_setting$proposal$incid_par_prop[2])
    incid_par_new[2] = exp(t_new2)

    IncidLog_new = log_incidence(MCMC_obj$LatentTraj, MCMC_setting$Incidence, Incid_Par = incid_par_new,F,T)

    pr_diff = dnorm(t_new,MCMC_setting$prior$rho_pr[1],MCMC_setting$prior$rho_pr[2],log = T) -
      dnorm(t_old,MCMC_setting$prior$rho_pr[1],MCMC_setting$prior$rho_pr[2],log = T) +
      dnorm(t_new2,MCMC_setting$prior$phi_pr[1],MCMC_setting$prior$phi_pr[2],log = T) -
      dnorm(t_old2,MCMC_setting$prior$phi_pr[1],MCMC_setting$prior$phi_pr[2],log = T)

    u = runif(1,0,1)
    if(log(u) < IncidLog_new - MCMC_obj$IncidLog + pr_diff){
      MCMC_obj$incid_par = incid_par_new
      MCMC_obj$IncidLog = IncidLog_new
    }
  }else{
    if(update[1] == 1){
    # rho
    incid_par_new = MCMC_obj$incid_par
    t_old = sigmoid_inv(MCMC_obj$incid_par[1])
    t_new = t_old + rnorm(1,0,MCMC_setting$proposal$incid_par_prop[1])
    incid_par_new[1] = sigmoid(t_new)

    IncidLog_new = log_incidence(MCMC_obj$LatentTraj, MCMC_setting$Incidence, Incid_Par = incid_par_new,F,T)
    #pr_diff = dgamma(t_new,MCMC_setting$prior$rho_pr[1],MCMC_setting$prior$rho_pr[2],log = T) - dgamma(t_old,MCMC_setting$prior$rho_pr[1],MCMC_setting$prior$rho_pr[2],log = T)
    pr_diff = dnorm(t_new,MCMC_setting$prior$rho_pr[1],MCMC_setting$prior$rho_pr[2],log = T) - dnorm(t_old,MCMC_setting$prior$rho_pr[1],MCMC_setting$prior$rho_pr[2],log = T)
    u = runif(1,0,1)
    if(log(u) < IncidLog_new - MCMC_obj$IncidLog + pr_diff){
      MCMC_obj$incid_par = incid_par_new
      MCMC_obj$IncidLog = IncidLog_new
      }
    }
    if(update[2] == 1){

      incid_par_new = MCMC_obj$incid_par
      t_old = log(MCMC_obj$incid_par[2])
      t_new = t_old + rnorm(1,0,MCMC_setting$proposal$incid_par_prop[2])
      incid_par_new[2] = exp(t_new)

      IncidLog_new = log_incidence(MCMC_obj$LatentTraj, MCMC_setting$Incidence, Incid_Par = incid_par_new,F,T)

      pr_diff = dnorm(t_new,MCMC_setting$prior$phi_pr[1],MCMC_setting$prior$phi_pr[2],log = T) - dnorm(t_old,MCMC_setting$prior$phi_pr[1],MCMC_setting$prior$phi_pr[2],log = T)
      u = runif(1,0,1)
      if(log(u) < IncidLog_new - MCMC_obj$IncidLog + pr_diff){
        MCMC_obj$incid_par = incid_par_new
        MCMC_obj$IncidLog = IncidLog_new
        }
      }
    }
  return(MCMC_obj)
}

#' @title MCMC for SIR model with LNA and changepoint on R0
#'
#' @param coal_obj Observed coalescent data. A list contains the sufficient statistics for coalescent: coalescant time,
#' sampling time and number of active lineages sampled at each sampling time. For a phylogentic tree object, this is summarize_phylo(tree)
#'
#' @param times A fine time grid used for Ode integration
#' @param t_correct The time period from start point to the last sampling time, a short t_correct can be used for infering
#' infectious deseased at early stage based on truncated tree
#' @param gridsize size of grid to integrate LNA. If dt is the interval of ODE integration, then then t_i - t_i-1
#' @param niter total number of iterations for MCMC
#' @param burn Number of iterations to burn
#' @param thin Number of iterations for thinning
#' @param changetime a vector of times that allows changing R0.
#' @param DEMS The name of each dimension for infectious disease trajectory
#' @param prior list of vectors for the values of the parameter prior
#' @param proposal The step size for M-H proposal
#' @param control The list of parameters that allows changing
#' @param ESS_vec A 0-1 vector indicate which parameters to update using ESS
#' @param likelihood Using volz likelihood by default
#' @param model SIR or SEIR model
#' @param Index Index for the S and I
#' @param nparam number of parameters for infectious disease model. 2 for SIR model, 3 for SEIR model
#' @param method sequentially update parameters of jointly update parameters
#' @param opitons updating details for MCMC
#' @param verbose boolean, if print parameters values during MCMC iterations
#'
General_MCMC_Integrate_ESlice = function(coal_obs,incidence_obs,times,t_correct,N,gridsize=1000, niter = 1000, burn = 0, thin = 5, changetime, DEMS=c("S","I"),
                                    prior=list(pop_pr=c(1,1,1,10), R0_pr=c(0.7,0.2), mu_pr = c(3,0.2), gamma_pr = c(3,0.2), hyper_pr = c(0.001,0.001)),
                                    proposal = list(pop_prop = 0.5, R0_prop = c(0.01), mu_prop=0.1, gamma_prop = 0.2, hyper_prop=0.05),
                                    control = list(), ESS_vec = c(1,1,1,1,1,1,1),likelihood = "volz",model = "SIR",
                                    Index = c(0,2), nparam=2, method = "seq",options = list(joint = F, PCOV = NULL,beta = 0.05, burn1 = 5000,
                                                                                            parIdlist = NULL, priorIdlist = NULL,up = 2000, tune = 0.01, hyper = F), enable = c(T,T), REAL = F,verbose = T){

  MCMC_setting = MCMC_setup_Joint_Inference(coal_obs, incidence_obs, times,t_correct,N,gridsize,niter,burn,
                                    thin,changetime, DEMS,prior,proposal,
                                    control,likelihood,model,Index,nparam, options$PCOV,REAL = REAL)


  p = MCMC_setting$p
  nparam = sum(MCMC_setting$x_i[1:2]) + 1

  MCMC_obj = MCMC_initialize_Integrate(MCMC_setting)

  prlist = list(a = prior[[1]][1:2], b = prior[[2]][1:2],c = prior[[3]][1:2],d = prior[[5]][1:2])
  params_incid = NULL
  if(enable[2]){
    params_incid = matrix(nrow = niter, ncol = 2)
  }
  if (MCMC_setting$likelihood == "volz") { # volz likelihood model

    params = matrix(nrow = niter, ncol = nparam + MCMC_obj$p)
  } else if (MCMC_setting$likelihood == "structural") { # structured coalescent model

    params = matrix(nrow = niter, ncol =  nparam + MCMC_obj$p)

  }else{ # other models
    params = matrix(nrow = niter, ncol = sum(MCMC_setting$x_i) + MCMC_obj$p)
  }

  l = numeric(niter)
  l1 = l
  l3 = l
  #l2 = matrix(ncol = 5, nrow = niter)
  tjs = array(dim = c(dim(MCMC_obj$LatentTraj),niter/thin))
  #if(length(times) > MCMC_setting$x_i[6]){
    #PredMX = matrix(nrow = niter, ncol = length(times) - MCMC_setting$x_i[6])
 # }else{
  #  PredMX = NULL
 # }

  #' updateVec
  #' 1 parameter for intial state
  #' 2 R0
  #' 3 gamma
  #' 4 mu
  #' 5 hyper-parameter controls the smoothness of changepoints
  #' 6 changepoints
  #' 7 LNA noise in trajectory

  ESS_temp = ESS_vec
  ESS_temp[5] = 0
  for (i in 1:MCMC_setting$niter) {
    if(i %% 10000 == 0){
      print(i)
    }
    if (i %% 100 == 0 && verbose == T) {
      print(i)
      print(MCMC_obj$par)
      print(l1[i-1])
      print(l2[i-1])

      plot(MCMC_obj$LatentTraj[,1], MCMC_obj$LatentTraj[, MCMC_setting$x_i[4] + 2],type="l",xlab = "time", ylab = "Infected")
      lines(MCMC_obj$Ode_Traj_coarse[,1],MCMC_obj$Ode_Traj_coarse[, MCMC_setting$x_i[4] + 2],col="red",lty=2)
    }

    if(i < options$burn1){
      # update gamma,
      MCMC_obj = tryCatch({Update_Param_each_Norm(MCMC_obj, MCMC_setting, options$parIdlist,
                                                  options$isLoglist, options$priorIdlist,options$priorIdlist,hyper = F,joint = T, enable = enable)$MCMC_obj
      }, error = function(cond){
        message(cond)
        # Choose a return value in case of error
        return(MCMC_obj)
      })

      # update R0 chpts
      if(ESS_vec[6] != 0){
        MCMC_obj = tryCatch({update_Par_ESlice_combine2(MCMC_obj,MCMC_setting, priorList = prlist, ESS_vec = ESS_temp, i, enable = enable)},
                            error = function(cond){
                              message(cond)
                              # Choose a return value in case of error
                              return(MCMC_obj)
                            })
      }

      MCMC_obj = Update_incidence_Pars(MCMC_setting, MCMC_obj, update = c(1,1), F)
    }else{
      # Metropolis-Hasting step, update I_0 gamma
      MCMC_obj = tryCatch({Update_Param_each_Norm(MCMC_obj, MCMC_setting, options$parIdlist,
                                                  options$isLoglist, options$priorIdlist,options$priorIdlist,hyper = options$hyper,joint = T, enable = enable)$MCMC_obj
      }, error = function(cond){
        message(cond)
        # Choose a return value in case of error
        return(MCMC_obj)
      })

      if(ESS_vec[6] != 0){
      # update R_0, changepoints and hyperparameter
        MCMC_obj = tryCatch({update_Par_ESlice_combine2(MCMC_obj,MCMC_setting,prlist,ESS_vec,i, enable = enable)},
                            error = function(cond){
                              message(cond)
                              # Choose a return value in case of error
                              return(MCMC_obj)
                            })
      }
      # update trajectories
      MCMC_obj = tryCatch({updateTraj_general_NC(MCMC_obj,MCMC_setting,i,T, enable = enable)$MCMC_obj},
                            error = function(cond){
                              message(cond)
                              # Choose a return value in case of error
                              return(MCMC_obj)
                            })
      if(enable[2]){
        MCMC_obj = Update_incidence_Pars(MCMC_setting, MCMC_obj, update = c(1,1), F)
      }
    }
    if(i %% thin == 0){
      tjs[,,as.integer(i/thin)] = MCMC_obj$LatentTraj
    }
  #  if(length(times) > MCMC_setting$x_i[6]){
  #    PredMX[i,] = LNA_integrate_pred(LNA_traj = MCMC_obj$LatentTraj, incid_par = MCMC_obj$incid_par, MCMC_setting$x_i[6])
#  }
    params[i,] = MCMC_obj$par
    if(enable[2]){
      params_incid[i,] = MCMC_obj$incid_par
    }
    l[i] = MCMC_obj$logOrigin
    l1[i] = MCMC_obj$coalLog
    l3[i] = MCMC_obj$IncidLog
  }
  return(list(par = params, incid_par = params_incid, Trajectory = tjs,l=l,l1=l1,l3 = l3, MX = MCMC_setting$PCOV, MCMC_setting = MCMC_setting, MCMC_obj = MCMC_obj))
}

log_post = function(MCMC_res, pref = F, incid = F){
  MCMC_setting = MCMC_res$MCMC_setting
  res = MCMC_res$l + MCMC_res$l1 + MCMC_res$l3
  d = dim(MCMC_res$par)[2]
  logpar = sapply(MCMC_res$par[,3],function(x) return(dlnorm(x, MCMC_setting$prior$R0_pr[1], MCMC_setting$prior$R0_pr[2], log = T))) +
    sapply(MCMC_res$par[,4],function(x) return(dlnorm(x, MCMC_setting$prior$gamma_pr[1], MCMC_setting$prior$gamma_pr[2], log = T))) +
    sapply(MCMC_res$par[,d],function(x) return(dlnorm(x, MCMC_setting$prior$hyper_pr[1], MCMC_setting$prior$hyper_pr[2], log = T)))
  if(incid){
    logpar = logpar + sapply(MCMC_res$incid_par[,1],function(x) return(dnorm(sigmoid_inv(x), MCMC_setting$prior$rho_pr[1], MCMC_setting$prior$rho_pr[2], log = T))) +
      sapply(MCMC_res$incid_par[,2],function(x) return(dlnorm(x, MCMC_setting$prior$phi_pr[1], MCMC_setting$prior$phi_pr[2], log = T)))
  }
  if(pref){
    logpar = logpar + sapply(MCMC_res$pref_par[,1],function(x) return(dnorm(x, MCMC_setting$prior$a_pr[1], MCMC_setting$prior$a_pr[2], log = T))) +
      sapply(MCMC_res$pref_par[,2],function(x) return(dnorm(x, MCMC_setting$prior$b_pr[1], MCMC_setting$prior$b_pr[2], log = T)))
  }
  return(logpar + res)
}
####################

General_MCMC_Integrate_ODE = function(coal_obs,incidence_obs,times,t_correct,N,gridsize=1000, niter = 1000, burn = 0, thin = 5, changetime, DEMS=c("S","I"),
                                         prior=list(pop_pr=c(1,1,1,10), R0_pr=c(0.7,0.2), mu_pr = c(3,0.2), gamma_pr = c(3,0.2), hyper_pr = c(0.001,0.001)),
                                         proposal = list(pop_prop = 0.5, R0_prop = c(0.01), mu_prop=0.1, gamma_prop = 0.2, hyper_prop=0.05),
                                         control = list(), ESS_vec = c(1,1,1,1,1,1,1),likelihood = "volz",model = "SIR",
                                         Index = c(0,2), nparam=2, method = "seq",options = list(parIdlist = NULL, priorIdlist = NULL, hyper = F), enable = c(T,F), REAL = F,verbose = T){

  MCMC_setting = MCMC_setup_Joint_Inference(coal_obs, incidence_obs, times,t_correct,N,gridsize,niter,burn,
                                            thin,changetime, DEMS,prior,proposal,
                                            control,likelihood,model,Index,nparam, options$PCOV,REAL = REAL)


  p = MCMC_setting$p
  nparam = sum(MCMC_setting$x_i[1:2]) + 1

  MCMC_obj = MCMC_initialize_Integrate(MCMC_setting)

  prlist = list(a = prior[[1]][1:2], b = prior[[2]][1:2],c = prior[[3]][1:2],d = prior[[5]][1:2])
  #params_incid = matrix(nrow = niter, ncol = 2)

  if (MCMC_setting$likelihood == "volz") { # volz likelihood model

    params = matrix(nrow = niter, ncol = nparam + MCMC_obj$p)
  } else if (MCMC_setting$likelihood == "structural") { # structured coalescent model

    params = matrix(nrow = niter, ncol =  nparam + MCMC_obj$p)

  }else{ # other models
    params = matrix(nrow = niter, ncol = sum(MCMC_setting$x_i) + MCMC_obj$p)
  }

  l1 = numeric(niter)
  l3 = l1
  #l2 = matrix(ncol = 5, nrow = niter)
  tjs = array(dim = c(dim(MCMC_obj$LatentTraj),niter/thin))



  #' updateVec
  #' 1 parameter for intial state
  #' 2 R0
  #' 3 gamma
  #' 4 mu
  #' 5 hyper-parameter controls the smoothness of changepoints
  #' 6 changepoints
  #' 7 LNA noise in trajectory

  ESS_temp = ESS_vec
  #ESS_temp[5] = 0
  for (i in 1:MCMC_setting$niter) {
    if(i %% 10000 == 0){
      print(i)
    }
    if (i %% 100 == 0 && verbose == T) {
      print(i)
      print(MCMC_obj$par)
      print(l1[i-1])
      print(l2[i-1])

      plot(MCMC_obj$LatentTraj[,1], MCMC_obj$LatentTraj[, MCMC_setting$x_i[4] + 2],type="l",xlab = "time", ylab = "Infected")
      lines(MCMC_obj$Ode_Traj_coarse[,1],MCMC_obj$Ode_Traj_coarse[, MCMC_setting$x_i[4] + 2],col="red",lty=2)
    }

      # update gamma,
      MCMC_obj = tryCatch({Update_Param_each_Norm(MCMC_obj, MCMC_setting, options$parIdlist,
                                                  options$isLoglist, options$priorIdlist,options$priorIdlist,hyper = F,joint = T, enable = enable)$MCMC_obj
      }, error = function(cond){
        message(cond)
        # Choose a return value in case of error
        return(MCMC_obj)
      })

      # update R0 chpts
      MCMC_obj = tryCatch({update_Par_ESlice_combine2(MCMC_obj,MCMC_setting, priorList = prlist, ESS_vec = ESS_temp, i, enable = enable)},
                          error = function(cond){
                            message(cond)
                            # Choose a return value in case of error
                            return(MCMC_obj)
                          })

      #MCMC_obj = Update_incidence_Pars(MCMC_setting, MCMC_obj, update = c(1,1), F)

    if(i %% thin == 0){
      tjs[,,as.integer(i/thin)] = MCMC_obj$LatentTraj
    }
    params[i,] = MCMC_obj$par
    #params_incid[i,] = MCMC_obj$incid_par
    l1[i] = MCMC_obj$coalLog
    l3[i] = MCMC_obj$IncidLog
  }
  return(list(par = params, Trajectory = tjs,l1=l1,l3 = l3, MX = MCMC_setting$PCOV, MCMC_setting = MCMC_setting, MCMC_obj = MCMC_obj))
}

