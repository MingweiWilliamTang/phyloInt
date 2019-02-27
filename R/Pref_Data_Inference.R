Sample_byIncid = function(Traj, beta1, beta0,left=NULL,right = NULL){
  n = dim(Traj)[1]
  if(is.null(left)){
    left = Traj[1,1]
  }
  if(is.null(right)){
    right = Traj[n,1]
  }
  res = apply(Traj, 1, function(x){
    return(c(x[1], rpois(1, lambda = exp(beta0 + beta1 * log(x[2] + 1)))))
  })
  res = t(res)
  res[(res[,1]< left) |(res[,1]>right) ,2] = 0
  samples = res[(which(res[,2]>0)),]
  t_max = samples[dim(samples)[1],1]
  t_correct = max(Traj[,1]) - t_max
  return(list(samp_time = t_max - rev(samples[,1]), n_sampled = rev(samples[,2]), t_correct = t_correct))
}



Update_Param_pref_each_Norm = function(MCMC_obj, MCMC_setting, parIdlist, isLoglist, priorIdlist, proposeIdlist,hyper = T, joint = T, enable = c(T, T)){

  if(is.null((parIdlist))){
    return(list(MCMC_obj = MCMC_obj))
  }

  if(!check_list_eq(parIdlist, isLoglist)){
    stop("parameters do not match")
  }

  if(!check_list_eq(parIdlist, priorIdlist)){
    stop("priors do not match")
  }

  if(!check_list_eq(parIdlist, proposeIdlist)){
    stop("priors do not match")
  }

  d = length(parIdlist)
  p = MCMC_obj$p
  x_i = MCMC_setting$x_i
  hyperId = (p + MCMC_setting$x_i[1] + MCMC_setting$x_i[2] + 1)

  for(i in 1:d){
    par_probs = MCMC_obj$par_probs
    par_new = MCMC_obj$par
    initial_new = MCMC_obj$par[1:p]
    param_new = MCMC_obj$par[-(1:p)]
    subd = length(parIdlist[[i]])
    parId = parIdlist[[i]]
    isLog = isLoglist[[i]]
    priorId = priorIdlist[[i]]
    proposeId = proposeIdlist[[i]]

    initialId = parId[parId <= p]
    paramId = parId[parId > p & parId <= (p + x_i[1] + x_i[2] + 1)]

    RawTransParam = MCMC_obj$par[parId]
    RawTransParam[isLog == 1] = log(RawTransParam[isLog == 1])
    RawTransParam_new = RawTransParam
    prior_proposal_offset = 0
    for(i in 1:length(parId)){

      newdiff = 0

      if(priorId[i] %in% c(1:4)){
        pr = MCMC_setting$prior[[ priorId[i] ]]
        po = MCMC_setting$proposal[[proposeId[i]]]
        RawTransParam_new[i] = RawTransParam[i] + runif(1,-po, po)
        newdiff = dnorm(RawTransParam_new[i], pr[1], pr[2], log = T) - dnorm(RawTransParam[i], pr[1], pr[2], log = T)
        par_probs[priorId[i]] = dlnorm(exp(RawTransParam_new[i]), pr[1], pr[2], log = T)
      }

      if(hyper == T){
        if(priorId[i] == 5){
          pr = MCMC_setting$prior[[priorId[i]]]
          po = MCMC_setting$proposal[[proposeId[i]]]
          u = runif(1,-po, po)
          RawTransParam_new[i] = RawTransParam[i] * exp(u)
          newdiff = dlnorm(RawTransParam_new[i], pr[1], pr[2], log = T) - u - dlnorm(RawTransParam[i], pr[1], pr[2], log = T)
          par_probs[priorId[i]] = dlnorm(RawTransParam_new[i], pr[1], pr[2], log = T)
        }
      }

      prior_proposal_offset = prior_proposal_offset + newdiff
    }

    RawTransParam_new[isLog == 1] = exp(RawTransParam_new[isLog == 1])

    if(hyperId %in% parId){
      par_new[(p + x_i[2] + 1): (hyperId - 1)] = exp(log(par_new[(p + x_i[2] + 1): (hyperId - 1)]) *
                                                       RawTransParam[subd] / RawTransParam_new[subd])
    }

    par_new[parId] = RawTransParam_new
    initial_new[initialId] = par_new[initialId]
    param_new = par_new[-(1:p)]
    update_res = NULL
    if(joint){
      update_res = Update_Param_Pref(param_new, initial_new, MCMC_obj$pref_par, MCMC_setting$times, MCMC_obj$OriginTraj,
                                          MCMC_setting$x_r, MCMC_setting$x_i, MCMC_setting$Init, MCMC_setting$SampleInit, MCMC_setting$gridsize, MCMC_obj$coalLog,MCMC_obj$PrefLog,prior_proposal_offset,
                                          MCMC_setting$t_correct, model = MCMC_setting$model,
                                          volz = MCMC_setting$likelihood == "volz", addCoal = enable[1], addPref = enable[2])
    }else{
      update_res = Update_Param(param_new, initial_new, MCMC_setting$times, MCMC_obj$OriginTraj,
                                MCMC_setting$x_r, MCMC_setting$x_i, MCMC_setting$Init, MCMC_setting$gridsize, MCMC_obj$coalLog,prior_proposal_offset,
                                MCMC_setting$t_correct, model = MCMC_setting$model,
                                volz = MCMC_setting$likelihood == "volz")

    }
    if(update_res$accept){
      MCMC_obj$par = par_new
      MCMC_obj$FT = update_res$FT
      MCMC_obj$Ode_Traj_coarse = update_res$Ode
      MCMC_obj$betaN = update_res$betaN
      MCMC_obj$coalLog = update_res$coalLog
      MCMC_obj$LatentTraj = update_res$LatentTraj
      MCMC_obj$par_probs = par_probs
      if(joint){
        #print(update_res$PrefLog - log_incidence(MCMC_obj$LatentTraj, MCMC_setting$SampleInit,MCMC_obj$pref_par))
        MCMC_obj$PrefLog = update_res$PrefLog
      }
    }

  }
  return(list(MCMC_obj = MCMC_obj, par_probs = par_probs))
}


##########################

update_Par_ESlice_combine3 = function(MCMC_obj, MCMC_setting, priorList, ESS_vec,i,enable=c(T,T)){
  p = MCMC_setting$p
  #print(paste("==",MCMC_obj$PrefLog - log_incidence(MCMC_obj$LatentTraj, MCMC_setting$SampleInit, MCMC_obj$pref_par)))
  ESlice_Result = ESlice_par_General_pref(MCMC_obj$par, MCMC_obj$pref_par, MCMC_setting$times, MCMC_obj$OriginTraj,
                                               priorList, MCMC_setting$x_r, MCMC_setting$x_i, MCMC_setting$Init, MCMC_setting$SampleInit,
                                               MCMC_setting$gridsize, ESS_vec, coal_log = MCMC_obj$coalLog, pref_log = MCMC_obj$PrefLog, MCMC_setting$t_correct,
                                               addCoal = enable[1],addPref = enable[2])
  MCMC_obj$par = ESlice_Result$par
  MCMC_obj$LatentTraj = ESlice_Result$LatentTraj
  MCMC_obj$betaN = ESlice_Result$betaN
  MCMC_obj$FT = ESlice_Result$FT

  if(ESlice_Result$CoalLog + ESlice_Result$PrefLog - MCMC_obj$PrefLog - MCMC_obj$coalLog < -25){
    print(paste("ChangePoint slice sampling problem",i," th"))
  }
  MCMC_obj$coalLog = ESlice_Result$CoalLog
  MCMC_obj$Ode_Traj_coarse = ESlice_Result$OdeTraj
  MCMC_obj$par_probs = prob_Eval(MCMC_obj$par, priorList, sum(MCMC_setting$x_i[1:2]) + p + 1)
  MCMC_obj$PrefLog = ESlice_Result$PrefLog
  return(MCMC_obj)
}



############################


updateTraj_general_NC2 = function(MCMC_obj,MCMC_setting,i, pref = F, enable = c(T,T)){
  new_CoalLog = 0
  if(MCMC_setting$likelihood == "structural"){

    Res = ESlice_general_NC_Structural(MCMC_obj$OriginTraj, MCMC_obj$Ode_Traj_coarse,
                                       MCMC_obj$FT, MCMC_setting$Init_Detail,
                                       MCMC_obj$par[(MCMC_obj$p+1):(MCMC_obj$p + MCMC_setting$x_i[1] + MCMC_setting$x_i[2])],
                                       MCMC_setting$x_r, MCMC_setting$x_i,
                                       coal_log = MCMC_obj$coalLog,
                                       model = MCMC_setting$model)
    MCMC_obj$coalLog = Res$CoalLog
  }else{

    if(pref){
      Res = ESlice_general_NC_pref(MCMC_obj$OriginTraj,MCMC_obj$Ode_Traj_coarse,
                                    MCMC_obj$FT, MCMC_obj$par[1:MCMC_obj$p], MCMC_obj$pref_par, MCMC_setting$Init, MCMC_setting$SampleInit,
                                    betaN = betaTs(MCMC_obj$par[(MCMC_obj$p+1):(MCMC_setting$x_i[1]+MCMC_setting$x_i[2]+MCMC_obj$p)],MCMC_obj$LatentTraj[,1], MCMC_setting$x_r,MCMC_setting$x_i),
                                    MCMC_setting$t_correct,lambda = 1,
                                    coal_log = MCMC_obj$coalLog, PrefLog = MCMC_obj$PrefLog, MCMC_setting$gridsize,
                                    volz = (MCMC_setting$likelihood == "volz"), model = MCMC_setting$model,addCoal = enable[1], addPref = enable[2])

    }else{
      Res = ESlice_general_NC(MCMC_obj$OriginTraj,MCMC_obj$Ode_Traj_coarse,
                              MCMC_obj$FT, MCMC_obj$par[1:MCMC_obj$p], MCMC_setting$Init,
                              betaN = betaTs(MCMC_obj$par[(MCMC_obj$p+1):(MCMC_setting$x_i[1]+MCMC_setting$x_i[2]+MCMC_obj$p)],MCMC_obj$LatentTraj[,1], MCMC_setting$x_r,MCMC_setting$x_i),
                              MCMC_setting$t_correct,lambda = 1,
                              coal_log = MCMC_obj$coalLog, MCMC_setting$gridsize,
                              volz = (MCMC_setting$likelihood == "volz"), model = MCMC_setting$model)

    }
    if(MCMC_setting$likelihood == "volz"){
      new_CoalLog = volz_loglik_nh_adj(MCMC_setting$Init, Res$LatentTraj,
                                       betaTs(MCMC_obj$par[(MCMC_obj$p+1):(MCMC_setting$x_i[1]+MCMC_setting$x_i[2]+MCMC_obj$p)],MCMC_obj$LatentTraj[,1], MCMC_setting$x_r,MCMC_setting$x_i),
                                       MCMC_setting$t_correct,
                                       index = MCMC_setting$x_i[3:4],enable = enable[1])

    }else{
      new_CoalLog = coal_loglik(MCMC_setting$Init,LogTraj(Res$LatentTraj ),MCMC_setting$t_correct,
                                MCMC_obj$par[5],MCMC_setting$gridsize)
    }
  }
  if(!pref && new_CoalLog - MCMC_obj$coalLog < -20){
    print(paste("problem with eslice traj" , i))
    print(paste("compare list res", new_CoalLog - Res$CoalLog))
  }
  if(pref && new_CoalLog - MCMC_obj$coalLog + Res$PrefLog - MCMC_obj$PrefLog < -20){
    print(paste("problem with eslice traj" , i))
    print(paste("compare list res", new_CoalLog - Res$CoalLog))
  }
  MCMC_obj$coalLog = new_CoalLog
  MCMC_obj$LatentTraj = Res$LatentTraj
  MCMC_obj$OriginTraj = Res$OriginTraj
  MCMC_obj$logOrigin = Res$logOrigin
  if(pref){
    MCMC_obj$PrefLog = Res$PrefLog
  }

  return(list(MCMC_obj = MCMC_obj))
}

MCMC_setup_Pref_Inference = function(tree,times,t_correct,N,gridsize=50,niter = 1000,burn = 500,thin = 5, changetime,DEMS = c("E", "I"),
                                      prior=list(pop_pr=c(1,10), R0_pr=c(1,7), mu_pr = c(3,0.2), gamma_pr = c(3,0.2), ch_pr = 1,hyper_pr=c(0.01,0.01)),
                                      proposal = list(pop_prop = 1, R0_prop = c(0.01), mu_prop=0.1, gamma_prop = 0.2, ch_prop=0.05),
                                      control = list(), likelihood = "volz", model = "SIR", Index = c(0,1), nparam = 2,REAL = F, cut = list(start = NULL, end = NULL)){

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
  }

  Init = coal_lik_init(coal_obs$samp_times, coal_obs$n_sampled, coal_obs$coal_times, grid,t_correct)
  SampleInit = preferential_lik_init(grid, coal_obs$samp_times, coal_obs$n_sampled, t_correct,gdiff = 1, g_start = cut$start, g_end = cut$end)

  x_i = c(length(changetime),nparam,Index[1],Index[2])
  print("time grids")
  print("_________________")
  MCMC_setting = list(Init = Init, SampleInit = SampleInit, times = Times,t_correct = t_correct,x_r = c(N,changetime),
                      gridsize = gridsize,gridset = gridset, niter = niter,burn = burn,thin = thin,x_i = x_i,
                      prior = prior, proposal = proposal, control = control, p = length(DEMS),
                      reps = 1, likelihood = likelihood, model = model)
  cat("MCMC set up ready \n")

  return(MCMC_setting)
}


MCMC_initialize_Pref = function(MCMC_setting, enable = c(T,T)){

  logMultiNorm = NaN
  coalLog = NaN
  PrefLog = NaN
  ########
  par_probs = numeric(5)
  while(is.nan(logMultiNorm)||is.nan(coalLog) || is.nan(PrefLog)){
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
    pref_par = NULL
    if(!is.null(MCMC_setting$control$pref_par)){
      pref_par = MCMC_setting$control$pref_par
    }
    par_probs[5] = dgamma(hyper, MCMC_setting$prior$hyper_pr[1],MCMC_setting$prior$hyper_pr[2], log = T)

    if(length(MCMC_setting$x_r) > 1){
      if(is.null(MCMC_setting$control$ch)){
        ch = rlnorm(MCMC_setting$x_i[1],0,1/hyper)
      }else{
        ch = MCMC_setting$control$ch
      }
    }
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
    #PrefLog = log_preverential(LatentTraj, MCMC_setting$SampleInit$SampleInit_obs, pref_par = pref_par)
    PrefLog = log_prev(LatentTraj, MCMC_setting$SampleInit, pref_par, Pois = T,enable = enable[2])
    print(paste("coalescent likelihood after initialization ", coalLog))

    if(!is.nan((coalLog))){
      coalLog = ifelse(coalLog <= - 10000000,NaN, coalLog)
    }
    if(!is.nan((PrefLog))){
      PrefLog = ifelse(PrefLog <= - 10000000,NaN, PrefLog)
    }
    plot(LatentTraj[,1],LatentTraj[,3],type="l",xlab = "time",col="blue", lwd = 2)
    #lines(LatentTraj[,1],LatentTraj[,3],col="red", lwd = 2)
  }

  chprobs = -0.5 * sum((log(param[(MCMC_setting$x_i[2] + 1):(MCMC_setting$x_i[1] + MCMC_setting$x_i[2])]) * hyper)^2)
  #LogGamma = dlnorm(gamma,MCMC_setting$prior$gamma_pr[1],MCMC_setting$prior$gamma_pr[2],log = T)
  #LogE = dlnorm(mu, MCMC_setting$prior$mu_pr[1], MCMC_setting$prior$mu_pr[2],log = T)
  plot(LatentTraj[,1],LatentTraj[,3],type="l",xlab = "time",col="blue", lwd = 2)

  paras =  c(state,param)
  MCMC_obj = list(par = paras,pref_par = pref_par,LatentTraj = LatentTraj,logMultiNorm = logMultiNorm,p = MCMC_setting$p,
                  Ode_Traj_coarse = Ode_Traj_coarse, FT = FT, coalLog = coalLog, PrefLog = PrefLog, OriginTraj = OriginTraj,logOrigin = logOrigin,
                  par_probs = par_probs, chprobs = chprobs, betaN = betaN)
  ##########
  cat("Initialize MCMC \n")
  print(paste("population size = ", MCMC_setting$N))
  print(paste(paras))
  return(MCMC_obj)
}





Update_Pref_Pars = function(MCMC_setting, MCMC_obj,update=c(1,1,0),joint = F){
  npar = length(MCMC_obj$pref_par)

  if(joint){
    pref_par_new = MCMC_obj$pref_par
    t_old = MCMC_obj$pref_par[1]
    t_new = t_old + rnorm(1,0,MCMC_setting$proposal$pref_par_prop[1])
    pref_par_new[1] = t_new

    t_old1.5 = MCMC_obj$pref_par[2]
    t_new1.5 = t_old1.5 + rnorm(1,0,MCMC_setting$proposal$pref_par_prop[2])
    pref_par_new[2] = t_new1.5

    t_old2 = log(MCMC_obj$pref_par[3])
    t_new2 = t_old2 + rnorm(1,0,MCMC_setting$proposal$pref_par_prop[3])
    pref_par_new[3] = exp(t_new2)

    PrefLog_new = log_prev(MCMC_obj$LatentTraj, MCMC_setting$SampleInit, pref_par_new,T,T)

    pr_diff = dnorm(t_new,MCMC_setting$prior$a_pr[1],MCMC_setting$prior$a_pr[2],log = T) -
      dnorm(t_old,MCMC_setting$prior$a_pr[1],MCMC_setting$prior$a_pr[2],log = T) +
      dnorm(t_new1.5,MCMC_setting$prior$b_pr[1],MCMC_setting$prior$b_pr[2],log = T) -
      dnorm(t_old1.5,MCMC_setting$prior$b_pr[1],MCMC_setting$prior$b_pr[2],log = T) +
      dnorm(t_new2, MCMC_setting$prior$phi_pr[1], MCMC_setting$prior$phi_pr[2], log = T) -
      dnorm(t_old2, MCMC_setting$prior$phi_pr[1], MCMC_setting$prior$phi_pr[2], log = T)

    u = runif(1,0,1)
    if(log(u) < PrefLog_new - MCMC_obj$PrefLog + pr_diff){
      MCMC_obj$pref_par = pref_par_new
      MCMC_obj$PrefLog = PrefLog_new
    }
  }else{
    if(update[1] == 1){
      # rho
      pref_par_new = MCMC_obj$pref_par
      t_old = MCMC_obj$pref_par[1]
      t_new = t_old + rnorm(1,0,MCMC_setting$proposal$pref_par_prop[1])
      pref_par_new[1] = t_new

      PrefLog_new = log_prev(MCMC_obj$LatentTraj, MCMC_setting$SampleInit, pref_par_new,T,T)

      pr_diff = dnorm(t_new,MCMC_setting$prior$a_pr[1],MCMC_setting$prior$a_pr[2],log = T) -
        dnorm(t_old,MCMC_setting$prior$a_pr[1],MCMC_setting$prior$a_pr[2],log = T)

      u = runif(1,0,1)

      if(log(u) < PrefLog_new - MCMC_obj$PrefLog + pr_diff){
        MCMC_obj$pref_par = pref_par_new
        MCMC_obj$PrefLog = PrefLog_new
      }
    }

    if(update[2] == 1){

      pref_par_new = MCMC_obj$pref_par
      t_old = MCMC_obj$pref_par[2]
      t_new = t_old + rnorm(1,0,MCMC_setting$proposal$pref_par_prop[2])
      pref_par_new[2] = t_new

      PrefLog_new = log_prev(MCMC_obj$LatentTraj, MCMC_setting$SampleInit, pref_par_new,T,T)

      pr_diff = dnorm(t_new, MCMC_setting$prior$b_pr[1],MCMC_setting$prior$b_pr[2],log = T) -
        dnorm(t_old, MCMC_setting$prior$b_pr[1],MCMC_setting$prior$b_pr[2],log = T)

      u = runif(1,0,1)

      if(log(u) < PrefLog_new - MCMC_obj$PrefLog + pr_diff){
        MCMC_obj$pref_par = pref_par_new
        MCMC_obj$PrefLog = PrefLog_new
      }
    }

    if(update[3] == 1){

      pref_par_new = MCMC_obj$pref_par

      t_old = log(MCMC_obj$pref_par[3])
      t_new = t_old + rnorm(1,0,MCMC_setting$proposal$pref_par_prop[3])
      pref_par_new[3] = exp(t_new)

      PrefLog_new = log_prev(MCMC_obj$LatentTraj, MCMC_setting$SampleInit, pref_par_new,T,T)

      pr_diff = dnorm(t_new, MCMC_setting$prior$phi_pr[1], MCMC_setting$prior$phi_pr[2], log = T) -
        dnorm(t_old, MCMC_setting$prior$phi_pr[1], MCMC_setting$prior$phi_pr[2], log = T)

      u = runif(1,0,1)
      if(log(u) < PrefLog_new - MCMC_obj$PrefLog + pr_diff){
        MCMC_obj$pref_par = pref_par_new
        MCMC_obj$PrefLog = PrefLog_new
      }
    }
  }
  return(MCMC_obj)
}




#################
General_MCMC_Pref_ESlice = function(coal_obs,times,t_correct,N,gridsize=1000, niter = 1000, burn = 0, thin = 5, changetime, DEMS=c("S","I"),
                                         prior=list(pop_pr=c(1,1,1,10), R0_pr=c(0.7,0.2), mu_pr = c(3,0.2), gamma_pr = c(3,0.2), hyper_pr = c(0.001,0.001)),
                                         proposal = list(pop_prop = 0.5, R0_prop = c(0.01), mu_prop=0.1, gamma_prop = 0.2, hyper_prop=0.05),
                                         control = list(), ESS_vec = c(1,1,1,1,1,1,1),likelihood = "volz",model = "SIR",
                                         Index = c(0,2), nparam=2, method = "seq",options = list(joint = F, PCOV = NULL,beta = 0.05, burn1 = 5000,
                                                                                                 parIdlist = NULL, priorIdlist = NULL,up = 2000, tune = 0.01, hyper = F), enable = c(T,T), REAL = F,verbose = T, cut = list(start = NULL, end = NULL)){

  MCMC_setting = MCMC_setup_Pref_Inference(coal_obs, times,t_correct,N,gridsize,niter,burn,
                                            thin,changetime, DEMS,prior,proposal,
                                            control,likelihood,model,Index,nparam, REAL = REAL, cut = cut)


  p = MCMC_setting$p
  nparam = sum(MCMC_setting$x_i[1:2]) + 1

  MCMC_obj = MCMC_initialize_Pref(MCMC_setting, enable = enable)

  prlist = list(a = prior[[1]][1:2], b = prior[[2]][1:2],c = prior[[3]][1:2],d = prior[[5]][1:2])
  params_pref = matrix(nrow = niter, ncol = 3)

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
      MCMC_obj = tryCatch({Update_Param_pref_each_Norm(MCMC_obj, MCMC_setting, options$parIdlist,
                                                  options$isLoglist, options$priorIdlist,options$priorIdlist,hyper = F,joint = T, enable = enable)$MCMC_obj
      }, error = function(cond){
        message(cond)
        # Choose a return value in case of error
        return(MCMC_obj)
      })

      # update R0 chpts
      MCMC_obj = tryCatch({update_Par_ESlice_combine3(MCMC_obj,MCMC_setting, priorList = prlist, ESS_vec = ESS_temp, i, enable = enable)},
                          error = function(cond){
                            message(cond)
                            # Choose a return value in case of error
                            return(MCMC_obj)
                          })

      MCMC_obj = Update_Pref_Pars(MCMC_setting, MCMC_obj, update = c(1,1,0), F)
    }else{
      # Metropolis-Hasting step, update I_0 gamma
      MCMC_obj = tryCatch({Update_Param_pref_each_Norm(MCMC_obj, MCMC_setting, options$parIdlist,
                                                  options$isLoglist, options$priorIdlist,options$priorIdlist,hyper = options$hyper,joint = T, enable = enable)$MCMC_obj
      }, error = function(cond){
        message(cond)
        # Choose a return value in case of error
        return(MCMC_obj)
      })

      # update R_0, changepoints and hyperparameter
      MCMC_obj = tryCatch({update_Par_ESlice_combine3(MCMC_obj,MCMC_setting,prlist,ESS_vec,i, enable = enable)},
                          error = function(cond){
                            message(cond)
                            # Choose a return value in case of error
                            return(MCMC_obj)
                          })
      # update trajectories
      MCMC_obj = tryCatch({updateTraj_general_NC2(MCMC_obj,MCMC_setting,i,T, enable = enable)$MCMC_obj},
                          error = function(cond){
                            message(cond)
                            # Choose a return value in case of error
                            return(MCMC_obj)
                          })
      if(enable[2]){
        MCMC_obj = Update_Pref_Pars(MCMC_setting, MCMC_obj, update = c(1,1,0), F)
      }
    }
    if(i %% thin == 0){
      tjs[,,as.integer(i/thin)] = MCMC_obj$LatentTraj
    }
    params[i,] = MCMC_obj$par
    params_pref[i,] = MCMC_obj$pref_par
    l[i] = MCMC_obj$logOrigin
    l1[i] = MCMC_obj$coalLog
    l3[i] = MCMC_obj$PrefLog
  }
  return(list(par = params, pref_par = params_pref, Trajectory = tjs,l=l,l1=l1,l3 = l3, MX = MCMC_setting$PCOV, MCMC_setting = MCMC_setting, MCMC_obj = MCMC_obj))
}


