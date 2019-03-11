#' Simulate a coalescent time based on volz likelihood
#'
#' @param eff_pop Three column matrix for population on a time grid. The first column is time. The second colunmn is number of susceptible and the third is number of infected, with both in log-scale.
#' @param t_correct The time that we observe the first sampling time
#' @param samp_time The time for the sampling event, reverse time scale
#' @param n_sampled Number of sampled event corresponse to each sampling time
#' @param betaN A vector of infection rate beta at the same time grid as eff_pop
#'
#' @return A list of sampling time, number sampled, interval of time points and coalsecent time inverse to the time scale

volz_thin_sir <- function(eff_pop,t_correct,samp_times, n_sampled, betaN,lower = 1, ...)
{
  coal_times = NULL
  lineages = NULL

  curr = 1
  active_lineages = n_sampled[curr]
  time = samp_times[curr]

  traj = 1 / (betaN*eff_pop[,2] /((eff_pop[,3])/2))
  #traj[1]=0.00000000001
  ts = (t_correct-eff_pop[,1])
  active_lineages = n_sampled[curr]
  # time = 0
  i = which.min(ts>=0) - 1
  ### correction
  if(i == 0){
    i = dim(eff_pop)[1]
  }

  #####
  while (time <= max(samp_times) || active_lineages > 1)
  {
    if (active_lineages == 1)
    {
      curr = curr + 1
      #print(curr)
      active_lineages = active_lineages + n_sampled[curr]
      time = samp_times[curr]
    }

    time = time + rexp(1, 0.5*active_lineages*(active_lineages-1)/lower)

    #i = which.min(ts>=0) - 1
    while(time>ts[i]){
      i = i-1
      if(i==0) {
        print("i=0")
        print(length(coal_times))
        i = 1
        break
      }
    }
    # print(ts[i])
    thed = traj[i]
    if (curr < length(samp_times) && time >= samp_times[curr + 1])
    {
      curr = curr + 1
      #print(curr)
      active_lineages = active_lineages + n_sampled[curr]
      time = samp_times[curr]
    }
    else if (runif(1) <= lower/(thed) )
    {
      time = min(c(time,t_correct))
      coal_times = c(coal_times, time)
      lineages = c(lineages, active_lineages)
      active_lineages = active_lineages - 1
    }
  }

  return(list(coal_times = coal_times, lineages = lineages,
              intercoal_times = c(coal_times[1], diff(coal_times)),
              samp_times = samp_times, n_sampled = n_sampled))
}




volz_thin_sir_AR <- function(eff_pop,t_correct,samp_times, n_sampled, betaN,lower = 1, ...)
{
  success = 0
  while(success == 0){
    print(success)
    coal_times = NULL
    lineages = NULL

    curr = 1
    active_lineages = n_sampled[curr]
    time = samp_times[curr]

    traj = 1 / (betaN*eff_pop[,2] /((eff_pop[,3])/2))
    #traj[1]=0.00000000001
    ts = (t_correct-eff_pop[,1])
    active_lineages = n_sampled[curr]
    # time = 0
    i = which.min(ts>=0) - 1
    counter = 0
    success = tryCatch({
      while (time <= max(samp_times) || active_lineages > 1)
      {
        if(counter > 5){
          stop("need resampling")
        }
        if (active_lineages == 1)
        {
          curr = curr + 1
          #print(curr)
          active_lineages = active_lineages + n_sampled[curr]
          time = samp_times[curr]
        }

        time = time + rexp(1, 0.5*active_lineages*(active_lineages-1)/lower)

        #i = which.min(ts>=0) - 1
        while(time>ts[i]){
          i = i-1
          if(i==0) {
            print("i=0")
            print(length(coal_times))
            i = 1
            break
          }
        }
        # print(ts[i])
        thed = traj[i]
        if (curr < length(samp_times) && time >= samp_times[curr + 1])
        {
          curr = curr + 1
          #print(curr)
          active_lineages = active_lineages + n_sampled[curr]
          time = samp_times[curr]
        }
        else if (runif(1) <= lower/(thed) )
        {
          time = min(c(time,t_correct))
          coal_times = c(coal_times, time)
          lineages = c(lineages, active_lineages)
          active_lineages = active_lineages - 1
        }
      }
      1
    },error = function(cond){
      message(cond)
      # Choose a return value in case of error
      return(0)
    })
    print(success)
  }
  return(list(coal_times = coal_times, lineages = lineages,
              intercoal_times = c(coal_times[1], diff(coal_times)),
              samp_times = samp_times, n_sampled = n_sampled))
}


Volz_simulator = function(state, time, param1, x_r1, x_i1, t_correct,samp_times, n_sampled,lower = 1){
  SIR_Traj = simuSIRt_AR(state, time, param1, x_r1, x_i1)
  betaN = betaTs(param1,times = time, x_r = x_r1, x_i = x_i1)
  coal = volz_thin_sir(SIR_Traj, t_correct, samp_times, n_sampled, betaN)
  return(list(traj = SIR_Traj, coal = coal))
}





coalsim_thin_sir <- function(eff_pop,t_correct,samp_times, n_sampled, lambda=1,lower = 0.1, ...)
{
  coal_times = NULL
  lineages = NULL

  curr = 1
  active_lineages = n_sampled[curr]
  time = samp_times[curr]

  traj = eff_pop[,2]/lambda
  #traj[1]=0.00000000001
  ts = (t_correct-eff_pop[,1])
  active_lineages = n_sampled[curr]
  # time = 0
  i = which.min(ts>=0) - 1

  while (time <= max(samp_times) || active_lineages > 1)
  {
    if (active_lineages == 1)
    {
      curr = curr + 1
      active_lineages = active_lineages + n_sampled[curr]
      time = samp_times[curr]
    }

    time = time + rexp(1, 0.5*active_lineages*(active_lineages-1)/lower)

    i = which.min(ts>=0) - 1
    while(time>ts[i]){
      i = i-1
      if(i==0) {
        print("i=0")
        i = 1
        break
      }
    }
    thed = traj[i]
    if (curr < length(samp_times) && time >= samp_times[curr + 1])
    {
      curr = curr + 1
      active_lineages = active_lineages + n_sampled[curr]
      time = samp_times[curr]
    }
    else if (runif(1) <= lower/( thed) )
    {
      time = min(c(time,t_correct))
      coal_times = c(coal_times, time)
      lineages = c(lineages, active_lineages)
      active_lineages = active_lineages - 1
    }
  }

  return(list(coal_times = coal_times, lineages = lineages,
              intercoal_times = c(coal_times[1], diff(coal_times)),
              samp_times = samp_times, n_sampled = n_sampled))
}


coalsim_thin_sir_Hom = function(eff_pop,t_correct,n_sampled,lambda){
  coal_times = NULL
  lineages = NULL

  traj = eff_pop[,2]/lambda
  ts = (t_correct-eff_pop[,1])
  active_lineages = n_sampled
  time = 0
  i = which.min(ts>=0) - 1
  while (active_lineages > 1)
  {
    if(i==0) break
    time = time + rexp(1, 0.5*active_lineages*(active_lineages-1))

    while(time>ts[i]){
      i = i-1
      if(i==0) {
        i = 1
        break
      }
    }
    thred = (traj[i] + traj[i+1]) / 2
    if (runif(1) <= (1/ ( thred)) )
    {
      coal_times = c(coal_times, time)
      lineages = c(lineages, active_lineages)
      active_lineages = active_lineages - 1
    }
  }
  return(list(coal_times = coal_times, lineages = lineages,
              intercoal_times = c(coal_times[1], diff(coal_times)),
              samp_times = 0, n_sampled = n_sampled))

}




coal_loglik_hom = function(eff_pop,t_correct,gene,lambda){
  dt = eff_pop[2,1] - eff_pop[1,1]
  n = gene$lineages[1]
  traj = eff_pop[,2]
  ts = t_correct - eff_pop[,1]
  time = 0
  i = which.min(ts>=0) - 1
  j = 1
  loglik = 0
  while(n>1){
    # determine the time on the grid
    Ck = 0.5 * n * (n-1)
    time = gene$coal_times[j]
    integr = dt / (traj[i] * lambda)
    while(time > ts[i] ){
      i = i - 1
      integr = integr + dt/ (traj[i] * lambda)
    }
    loglik = loglik + log(Ck) - log(lambda) - log(traj[i]) - Ck * integr
    j = j+1
    n=n-1
  }
  return(loglik)
}

#' Generate the suffcient statistics for coalescent likelihood calculation
#'
#' @param Grid A vector of time points that we use to infer the trajectory. grid is in reverse time scale and should match the number the grid of trajectory. The first grid point corresponds to the first sampling event
#'
#' @return A list sufficient statistics for coalescent likelihood calculation:
#'
#'

coal_lik_init_tree = function(tree,grid,grid_enddate = NULL, tree_correct = NULL){
  trunc_tree = tree
  if(!is.null(tree_correct)){
    trunc_tree = TruncTree_Fun(tree,tree_correct)
  }
  tree_dates = sapply(trunc_tree$tip.label, function(x){
    return(strsplit(x,"[|\']")[[1]][7])
  })
  t_correct = as.numeric(as.Date(grid_enddate)-max(as.Date(tree_dates)))/365
  return(list(trunc_tree = trunc_tree,t_correct = t_correct))
}


coal_lik_init = function(samp_times, n_sampled, coal_times, grid, t_correct = NULL )
{
  # t_correct = gridend - last_sample times
  # t_correct gap between your last sampling time and the end of your grid

  samp_times = samp_times + t_correct
  coal_times = coal_times + t_correct

  gridlength = length(grid)
  gridOrigin = grid
  grid = max(grid) - rev(grid) # reverser the grid timeline

  ns = length(samp_times)
  nc = length(coal_times)
  #samp_times = samp_times[ns:1]
  #coal_times = coal_times[nc:1]

  #n0 = which.min(grid < t_correct)
  Tcoal = max(coal_times)
  if(Tcoal > max(grid)){
    stop("coalescent time out of grid")
  }
  n0 = which.min(grid < Tcoal)
  grid = grid[1:n0]

  #ng = length(grid)-1
  ng = n0

  if (length(samp_times) != length(n_sampled))
    stop("samp_times vector of differing length than n_sampled vector.")

  #  if (length(coal_times) != sum(n_sampled) - 1)
  #    stop("Incorrect length of coal_times: should be sum(n_sampled) - 1.")

  #  if (max(samp_times, coal_times) > max(grid))
  #   stop("Grid does not envelop all sampling and/or coalescent times.")

  t = sort(unique(c(samp_times, coal_times, grid)))
  l = rep(0, length(t))

  ######

  #######
  for (i in 1:ns)
    l[t >= samp_times[i]] = l[t >= samp_times[i]] + n_sampled[i]
  for (i in 1:nc)
    l[t >= coal_times[i]] = l[t >= coal_times[i]] - 1

  #print(l)

  if (sum((l < 1) & (t >= min(samp_times))) > 0)
    stop("Number of active lineages falls below 1 after the first sampling point.")

  mask = l > 0
  t = t[mask]
  l = head(l[mask], -1)
  ###########################
  # new code enabling complicated grid setup
  gridstart =  which.max(mask)

  ##########################
  gridrep = rep(0, ng-gridstart+1)
  for (i in gridstart : ng)
    gridrep[i - gridstart + 1] = sum(t > grid[i] & t <= grid[i+1])

  if(sum(gridrep) < length(l)){
    gridrep = c(length(l) - sum(gridrep), gridrep[-length(gridrep)])
    gridstart = gridstart - 1
    ng = ng - 1
  }

  C = 0.5 * l * (l-1)
  D = diff(t)

  y = rep(0, length(D))
  y[t[-1] %in% coal_times] = 1

  buckets = cut(x = samp_times, breaks = t,
                include.lowest = TRUE)
  tab <- aggregate(n_sampled ~ buckets, FUN = sum, labels = FALSE)
  count <- rep(0, length(D))
  count[as.numeric(tab$buckets)] <- tab$n_sampled
  count[head(t, -1) >= max(samp_times)] <- NA

  rep_idx = cumsum(gridrep)
  rep_idx = cbind(rep_idx-gridrep+1,rep_idx)
  grid_idx = c(0,cumsum(gridrep))
  print(ng)
  #if(gridrep[ng] == 0){
  #  grid_idx = grid_idx[1:(length(grid_idx)-1)]
  #}
  temp = gridlength - gridstart + 1
  gridstart = gridlength - ng + 1
  ng = temp
  if(gridOrigin[ng] - max(gridOrigin) + min(samp_times) < 0)
    stop("wrong with the grid setup")
  if(gridrep[length(gridrep)] == 0){
    gridrep = head(gridrep,-1)
    gridstart = gridstart + 1
  }

  if(gridstart > 1){
    gridstart = gridstart - 1
    ng = ng - 1
  }
  return(list(t=t, l=l, C=C, D=D, y=y, count=count, gridrep=gridrep,
              ns=sum(n_sampled), nc=nc, ng=ng, rep_idx=rep_idx,
              args=list(samp_times=samp_times, n_sampled=n_sampled, coal_times=coal_times,
                        grid=grid),
              grid_idx = grid_idx,gridstart = gridstart))
}



preferential_lik_init = function(grid, samp_times, n_sampled, t_correct, gdiff = 1, g_start = NULL, g_end = NULL,incidPref){
  samp_times = samp_times+t_correct
  l = length(grid)
  grid2 = max(grid) - rev(grid)
  start = which.max(grid2 >= min(samp_times))
  if(grid2[start] > min(samp_times)){
    start = start - 1
  }
  j = 0
  while(j <= length(grid)){
    if(grid2[start + j * gdiff] <= max(samp_times)){
      j = j + 1
    }else{
      break;
    }
  }
  end = start + j * gdiff
  n0 = j + 1
  grid2 = grid2[seq(start, end, by = gdiff)]
  sample_cout = numeric(n0)
  for(i in 1:n0){
    sample_cout[i] = sum(rep(samp_times,n_sampled) >= grid2[i])
  }
  sample_cout = c(sample_cout,0)
  sample_cout = -diff(sample_cout)
  sample_cout = rev(sample_cout)
  start = l - start + 1
  end = l - end + 1
  g_idx = seq(end, start, by = gdiff)
  if(sample_cout[1] == 0){
    sample_cout = sample_cout[-1]
  }
  ####y
  times = grid[g_idx]
  if(is.null(g_start)){
    g_start = min(times)
  }
  if(is.null(g_end)){
    g_end = max(times)
  }
  n_idx = which(times >= g_start & times <= g_end)
  times = times[n_idx]
  g_idx = g_idx[n_idx]
  sample_cout = sample_cout[n_idx[-length(n_idx)]]
  ####
  #return(list(sample_cout = sample_cout, g_idx = g_idx - 1, time = grid[g_idx],bound = max(g_idx)))
  return(list(sample_cout = sample_cout, g_idx = g_idx - 1, time = times, bound = max(g_idx)))
}





SEIR_incidence_Traj = function(Traj){
  n = dim(Traj)[1]
  incidenceTraj = matrix(ncol = 2,nrow = n - 1)
  incidenceTraj[,1] = Traj[2:n,1]
  for(i in 2:n){
    NSE = Traj[i-1,2] - Traj[i,2]
    NEI = Traj[i-1,3] - Traj[i,3]+ NSE
    incidenceTraj[i-1,2] = NEI
  }
  return(incidenceTraj)
}


TruncTree_Fun = function(bigTree, date){
  cutID = sapply(bigTree$tip.label, function(x){
    a = strsplit(x,"[|\']")[[1]][7]
    return(as.Date(a) <= as.Date(date))
  })
  return(drop.tip(bigTree, bigTree$tip.label[!cutID],subtree = F))
}


forwardSimulateODE_RST_core = function(param, LatentTraj, x_r, x_i, p, dt, Ngrid, beginning = F, tol=5){
  R_tstart = param[p + x_i[3] + 1] * prod(param[(p + x_i[2] + 1):(p+x_i[2]+x_i[1])])
  if(R_tstart > 0.999){
    return(Inf)
  }

  if(is.matrix(LatentTraj)){
    X1 = tail.matrix(LatentTraj,n=1)
  }else{
    X1 = LatentTraj
  }

  if(beginning == T){
    t = seq(LatentTraj[1,1], X1[1] + dt * Ngrid, by=dt)
    return(ODE_rk45_stop(X1[-1], t, param[-(1:p)], x_r, x_i, tol = tol));
  }else{
    t = seq(X1[1],X1[1] + dt * Ngrid, by=dt)
    return(ODE_rk45_stop(X1[-1], t, param[-(1:p)], x_r, x_i, tol = tol));
  }
}

