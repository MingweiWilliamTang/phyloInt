  #Rcpp::sourceCpp('src/SIRS_period.cpp')
#Rcpp::sourceCpp("src/SIR_BD.cpp")
#source('~/Dropbox/Research/GP/code/LNAPhyloDyn/R/coal_simulation.R', echo=TRUE)
#source('~/Dropbox/Research/GP/code/LNAPhyloDyn/R/SIR_LNA_standard.R', echo=TRUE)

plot3 = function(traj){
  dev.off()
  par(mfrow = c(2,2),mar=c(4,4,1,1))
  plot(traj[,1],traj[,2],type = "l", xlab = "time", ylab = "suscepitible1")
  plot(traj[,1],traj[,3],type = "l", xlab = "time", ylab = "Infected")
  plot(traj[,1],traj[,4],type = "l", xlab = "time", ylab = "Recovered")
  #plot(traj[,1],traj[,5],type = "l", xlab = "time", ylab = "reconvered")
}



effpopfun = function(Traj,beta=0,lambda=1, volz = FALSE){
  if(volz){
    return(1 /(2 * Traj[,3] * beta / Traj[,2]))
  }else{
    return(Traj[,3] / lambda)
  }

}

#' Plot the trajectory of effective population size
#'
#' @param MCMC_res returned MCMC list
#' @param idxs The index choosed to make plot
#' @param volz If use volz likelihood the coalescent modeling
#'
#'
#'
random_effpoptraj_line = function(MCMC_res,idxs,volz = FALSE){
  N = MCMC_res$Trajectory[1,2,1] + MCMC_res$Trajectory[1,3,1]
  #plot(N)
  if(volz){
    for(i in idxs){
      lines(MCMC_res$Trajectory[,1,i],
            effpopfun(MCMC_res$Trajectory[,,i],beta = MCMC_res$par[i,3] * N,volz = T),
            lwd = 0.5, col = "grey")
    }
  }else{
    for(i in idxs){
      lines(MCMC_res$Trajectory[,1,i],
            MCMC_res$Trajectory[,3,i] / MCMC_res$par[i,5],
            lwd = 0.5, col = "yellow")
    }
  }
}



medianCur = function(MCMC_obj,ids,scale=1,col="black",row=3,med = T,volz = F,lwd = 2,lid=7,lty = 2){
  if(volz == F){
    scale = MCMC_obj$par[ids,lid]
  }
  if(med){
  lines(MCMC_obj$Trajectory[,1,1],apply(MCMC_obj$Trajectory[,row,ids]/scale,1,median),
          col=col,lwd=lwd,lty=lty)
  }else{
  lines(MCMC_obj$Trajectory[,1,1],apply(MCMC_obj$Trajectory[,row,ids]/scale,1,mean),
        col=col,lwd=lwd,lty = lty)
  }
}


effpopfun = function(Traj,beta=0,lambda=1, volz = FALSE){
  if(volz){
    return(1/(2 * Traj[,2] * beta / Traj[,3]))
  }else{
    return(Traj[,3] / lambda)
  }
}




CI_Curve = function(MCMC_obj,ids,scale = 1, col = "black", fill_col = "grey", row = 3,method = "qtile",alpha=0.05,fill = T){

  if(method == "NormApp"){
    midcur = apply(MCMC_obj$Trajectory[,row,ids]/scale,1,mean)
    midSd = apply(MCMC_obj$Trajectory[,row,ids]/scale,1,sd)
    qt = qnorm(1-alpha/2)
    upper = midcur + qt * midSd
    lower = midcur - qt * midSd
    #lines(MCMC_obj$Trajectory[,1,1],upper,lty=2,
    #      col=col,lwd=2)

    #lines(MCMC_obj$Trajectory[,1,1],lower,lty=2,
    #      col=col,lwd=2)
  }else if(method == "qtile"){
    qt1 = 1 - alpha/2
    qt2 = alpha/2
    upper = apply(MCMC_obj$Trajectory[,row,ids]/scale,1,function(x){
      return(quantile(x,qt1))
    }
      )
    lower = apply(MCMC_obj$Trajectory[,row,ids]/scale,1,function(x){
      return(quantile(x,qt2))
    })
   # lines(MCMC_obj$Trajectory[,1,1],midcur + qt * midSd,lty=2,
    #      col=col,lwd=2)

    #lines(MCMC_obj$Trajectory[,1,1],midcur - qt * midSd,lty=2,
     #     col=col,lwd=2)
  }
  if(fill == T){
    polygon(x = c(MCMC_obj$Trajectory[,1,1],rev(MCMC_obj$Trajectory[,1,1])),
            y = c(upper,rev(lower)),col = fill_col,border = NA)
  }else{
    polygon(x = c(MCMC_obj$Trajectory[,1,1],rev(MCMC_obj$Trajectory[,1,1])),
            y = c(upper,rev(lower)), col = "black",border = F,angle = c(45,-45),density = c(40))
  }
}

CI_Curve3 = function(times,Trajectory,ids,scale = 1, col = "black", fill_col = "grey", row = 3,method = "qtile",alpha=0.05,fill = T){

  if(method == "NormApp"){
    midcur = apply(Trajectory[,row,ids]/scale,1,mean)
    midSd = apply(Trajectory[,row,ids]/scale,1,sd)
    qt = qnorm(1-alpha/2)
    upper = midcur + qt * midSd
    lower = midcur - qt * midSd
    #lines(MCMC_obj$Trajectory[,1,1],upper,lty=2,
    #      col=col,lwd=2)

    #lines(MCMC_obj$Trajectory[,1,1],lower,lty=2,
    #      col=col,lwd=2)
  }else if(method == "qtile"){
    qt1 = 1 - alpha/2
    qt2 = alpha/2
    upper = apply(Trajectory[,row,ids]/scale,1,function(x){
      return(quantile(x,qt1))
    })
    mid = apply(Trajectory[,row,ids]/scale,1,function(x){
      return(median(x,qt1))
    })
    lower = apply(Trajectory[,row,ids]/scale,1,function(x){
      return(quantile(x,qt2))
    })
    # lines(MCMC_obj$Trajectory[,1,1],midcur + qt * midSd,lty=2,
    #      col=col,lwd=2)

    #lines(MCMC_obj$Trajectory[,1,1],midcur - qt * midSd,lty=2,
    #     col=col,lwd=2)
  }
  if(fill == T){
    polygon(x = c(times,rev(times)),
            y = c(upper,rev(lower)),col = fill_col,border = NA)
    lines(times,mid,col="red",lwd=2,lty=2)
  }
}
CI_Curve_eff = function(MCMC_obj,ids,scale = 1, col = "black", fill_col = "grey", row = 3,method = "qtile",alpha=0.05,fill = T,likelihood ="volz"){
  if(likelihood == "volz"){
    Mx = 1/(2 * MCMC_obj$Trajectory[,2,ids]/MCMC_obj$Trajectory[,3,ids])
  }else{
    Mx = MCMC_obj$Trajectory[,row,ids]
  }
  midcur = apply(Mx/scale,1,mean)
  if(method == "NormApp"){
    midSd = apply(Mx/scale,1,sd)
    qt = qnorm(1-alpha/2)
    upper = midcur + qt * midSd
    lower = midcur - qt * midSd
    #lines(MCMC_obj$Trajectory[,1,1],upper,lty=2,
    #      col=col,lwd=2)

    #lines(MCMC_obj$Trajectory[,1,1],lower,lty=2,
    #      col=col,lwd=2)
  }else if(method == "qtile"){
    qt1 = 1 - alpha/2
    qt2 = alpha/2
    upper = apply(Mx/scale,1,function(x){
      return(quantile(x,qt1))
    }
    )
    lower = apply(Mx/scale,1,function(x){
      return(quantile(x,qt2))
    })
    # lines(MCMC_obj$Trajectory[,1,1],midcur + qt * midSd,lty=2,
    #      col=col,lwd=2)

    #lines(MCMC_obj$Trajectory[,1,1],midcur - qt * midSd,lty=2,
    #     col=col,lwd=2)
  }
  if(fill == T){
    polygon(x = c(MCMC_obj$Trajectory[,1,1],rev(MCMC_obj$Trajectory[,1,1])),
            y = c(upper,rev(lower)),col = fill_col,border = NA)
  }
  lines(MCMC_obj$Trajectory[,1,1],midcur,col = col,lwd=2,lty=2)
}




vlineCI = function(data,cred = T){
  if(cred){
    m = median(data)
    s1 = quantile(data,0.025)
    s2 = quantile(data,1 - 0.025)
    abline(v = s1,col="blue",lwd=2,lty=2)
    abline(v = s2,col="blue",lwd=2,lty=2)
  }
  else{
  s = sd(data)
  m = meam(data)
  abline(v = m + 1.96*s,col="blue",lwd=2,lty=2)
  abline(v = m - 1.96*s,col="blue",lwd=2,lty=2)
  }
  abline(v = m ,col="blue",lwd=2)
}

randomR0_traj = function(times,MCMC_obj,R0_id,col_id,idx,ylim=c(0,2),ylab = "",col = "black",lty=1,fill = rgb(1,0,0,0.3), main = "",xlab = "time",cex.lab=1.3,cex.axis = 1.2,xaxt = "s",add = F,add2 = F){
  R0 = MCMC_obj$par[idx,R0_id]
  R0_traj = matrix(ncol= length(col_id)+2, nrow = length(idx))
  R0_traj[,1] = R0
  for(i in 1:length(col_id)){
    R0_traj[,i+1] = R0_traj[,i] * MCMC_obj$par[idx,col_id[i]]
  }
  i = length(col_id)
  R0_traj[,length(col_id)+2] = R0_traj[,i] * MCMC_obj$par[idx,col_id[i]]
  CIup = apply(R0_traj,2,function(x){
    return(quantile(x,0.975))
  })
  CIlow = apply(R0_traj,2,function(x){
    return(quantile(x,0.025))
  })
  m = apply(R0_traj,2,median)
  if(add){
    lines(times,m,type="l",ylab = ylab,col = col,lty = lty,lwd = 2,ylim=ylim,main=main,xlab = xlab,
          cex.lab= cex.lab,cex.axis = cex.axis, xaxt = xaxt)
    polygon(x = c(times,rev(times)),
            y = c(CIup,rev(CIlow)),col = "black",border = F,angle = c(45,-45),density = c(40))
  }else{
    if(add2){
    lines(times,m,type="l",ylab = ylab,col = col,lty = lty,lwd = 2,ylim=ylim,main=main,xlab = xlab,
        cex.lab= cex.lab,cex.axis = cex.axis, xaxt = xaxt)
    }else{
      plot(times,m,type="l",ylab = ylab,col = col,lty = lty,lwd = 2,ylim=ylim,main=main,xlab = xlab,
            cex.lab= cex.lab,cex.axis = cex.axis, xaxt = xaxt)
    }
    polygon(x = c(times,rev(times)),
            y = c(CIup,rev(CIlow)),col = fill,border = NA)
  }

  #lines(times,m,type="l",col = "red",lwd = 2)
}


randomR0_traj_V = function(times,MCMC_obj,R0_id,col_id,idx,xlim,ylim=c(0,2),main = "",fill_col = rgb(1,0,0,0.3)){
  R0 = MCMC_obj$par[idx,R0_id]
  R0_traj = matrix(ncol= length(col_id)+1, nrow = length(idx))
  R0_traj[,1] = R0
  for(i in 1:length(col_id)){
    R0_traj[,i+1] = R0_traj[,i] * MCMC_obj$par[idx,col_id[i]]
  }
  i = length(col_id)
  #R0_traj[,length(col_id)+2] = R0_traj[,i] * MCMC_obj$par[idx,col_id[i]]
  CIup = apply(R0_traj,2,function(x){
    return(quantile(x,0.975))
    })
  CIlow = apply(R0_traj,2,function(x){
      return(quantile(x,0.025))
    })
  m = apply(R0_traj,2,median)
  times = times[2:(length(m))]
  step1 = stepfun(times, m, f = 0)
  step2 = stepfun(times, CIup, f = 0)
  step3 = stepfun(times, CIlow, f = 0)
  plot(step1,ylab = expression(R0(t)),col = "red",lwd = 2.5,ylim=ylim,
       main=main,verticals = F,xlim=xlim,xlab = "time",lty=2,xaxt = "n",cex.lab=1.3,cex.axis = 1.2)
  #polygon(x = c(times,rev(times)),
   #       y = c(CIup,rev(CIlow)),col = "grey",border = NA)
  polygon(x = c(c(0,rep(times,each=2),1.4), rev(c(0,rep(times,each=2),1.4))),
          c(rep(CIlow,each=2), rev(rep(CIup,each=2))),col = fill_col,border = NA)
  #lines(step2, lty=2,lwd = 1,verticals = F, col = "blue",xlim=xlim)
  #lines(step3, lty=2,lwd = 1,verticals = F, col = "blue",xlim=xlim)
  #lines(times,m,type="l",col = "red",lwd = 2)
}




randomR0s = function(times,MCMC_obj,R0_id,col_id,idx){

  R0 = MCMC_obj$par[idx,R0_id]
  R0_traj = matrix(ncol= length(col_id)+1, nrow = length(idx))
  R0_traj[,1] = R0
  for(i in 1:length(col_id)){
    R0_traj[,i+1] = R0_traj[,i] * MCMC_obj$par[idx,col_id[i]]
  }

  return(R0_traj)
}



coal_like_fast = function(traj, lambda = 1, coal_obj,t_correct,col=3){
  init = coal_lik_init(coal_obj$samp_times,coal_obj$n_sampled,coal_obj$coal_times,grid = traj[,1])
  return(coal_loglik(init,LogTraj(traj),t_correct = t_correct, lambda))
}

histwithCI = function(data,cred = T){
  hist(data);
  vlineCI(data,cred)
}


effpopPlot = function(NP_res,LNA_res,t_correct,idx,row=id,lid=5,ylab="",xlab="",ylim=c(0,10)){
  plot(t_correct - NP_res$x,NP_res$effpopmean,type="l",ylim = ylim,xlab = xlab,ylab=ylab)
  polygon(c(t_correct - NP_res$x, rev(t_correct - NP_res$x)), c(NP_res$effpop975,rev(NP_res$effpop025)),
          col=rgb(0,0,1,0.3),border = F)
  CI_Curve(LNA_res,idx,scale = LNA_res$par[idx,5],fill_col = rgb(1,0,0,0.3))
  medianCur(LNA_res,idx,col="red",lid = 5)
}

volz_eff = function(Traj, param, x_r, x_i){
  x = 1 / (2 * Traj[,2]/Traj[,3])
  betas = betaTs(param, Traj[,1], x_r, x_i)
  return(cbind(Traj[,1],x/betas))
}

#' @title plot effective population size using output from MCMC
#'
#' @param MCMC_obj A list returned from MCMC functions
#' @param ids Indices from the MCMC iteration output
#' @thin Thinning applies on the SIR trajectory
#'
#'
CI_Curve_eff2 = function(MCMC_obj,ids, thin = 1,col = "black", fill_col = "grey", Irow = 3,method = "qtile",alpha=0.05,fill = T,likelihood ="volz",
                         x_r,x_i,p=3){
  if(sum(ids %% thin) == 1){
    stop("Some indices not stored")
  }
  if(likelihood == "volz"){
    Mx = 1/(2 * MCMC_obj$Trajectory[,2,ids / thin]/MCMC_obj$Trajectory[,Irow,ids / thin])
    scale = sapply(ids,function(x){
      return(betaTs(MCMC_obj$par[x,(p+1):(p + x_i[1] + x_i[2])],
                   MCMC_obj$Trajectory[,1,1],x_r, x_i))})
    Mx = Mx/scale
    }else{
    Mx = MCMC_obj$Trajectory[,Irow,ids / thin]
    scale = 1
  }
  midcur = apply(Mx,1,mean)
  if(method == "NormApp"){
    midSd = apply(Mx,1,sd)
    qt = qnorm(1-alpha/2)
    upper = midcur + qt * midSd
    lower = midcur - qt * midSd
    #lines(MCMC_obj$Trajectory[,1,1],upper,lty=2,
    #      col=col,lwd=2)

    #lines(MCMC_obj$Trajectory[,1,1],lower,lty=2,
    #      col=col,lwd=2)
  }else if(method == "qtile"){
    qt1 = 1 - alpha/2
    qt2 = alpha/2
    upper = apply(Mx,1,function(x){
      return(quantile(x,qt1))
    }
    )
    lower = apply(Mx,1,function(x){
      return(quantile(x,qt2))
    })
    # lines(MCMC_obj$Trajectory[,1,1],midcur + qt * midSd,lty=2,
    #      col=col,lwd=2)

    #lines(MCMC_obj$Trajectory[,1,1],midcur - qt * midSd,lty=2,
    #     col=col,lwd=2)
  }
  if(fill == T){
    polygon(x = c(MCMC_obj$Trajectory[,1,1],rev(MCMC_obj$Trajectory[,1,1])),
            y = c(upper,rev(lower)),col = fill_col,border = NA)
  }
  lines(MCMC_obj$Trajectory[,1,1],midcur,col = col,lwd=2,lty=2)
}



NP_plot = function(coal_obj, ngrids, t_correct, main = "", xlim = c(0,100), ylim = c(0,10)){
  #par(mgp = c(2.5,1,0),mar=c(4,4,1,1))
  plot(1, type="n", xlab="time", ylab="effective population size",
       xlim=xlim, ylim=ylim,cex.lab=1.3,cex.axis = 1.2,xaxt = "n")
  NP_res = BNPR(coal_obj, ngrids)
  lines(t_correct - NP_res$x, NP_res$effpopmean, type="l" , lwd=2, col = "black", lty = 3)
  polygon(c(t_correct - NP_res$x,rev(t_correct - NP_res$x)),
          c(NP_res$effpop025,rev(NP_res$effpop975)),
          col = "black",border = F,angle = c(45,-45),density = c(10, 20))
}


random_Incidence = function(MCMC_obj, idx,col = "black", fill_col = rgb(0.15,0.15,0.15,0.3), row = 3,lty = 2,ylim = c(0,1000),ylab = "incidence"){
  n = length(MCMC_obj$Trajectory[,1,1])
  t = SIR_incidence_Traj(MCMC_obj$Trajectory[,,1],0:(n-1))[,1]
  incid_MX = apply(MCMC_obj$Trajectory[,,idx], 3, function(x){
    return(SIR_incidence_Traj(x,0:(n-1))[,2])
  })
  mid = apply(incid_MX,1,median)
  upper = apply(incid_MX,1,function(x){
    return(quantile(x,0.975))
  })
  lower = apply(incid_MX,1,function(x){
    return(quantile(x,0.025))
  })
  plot(t, mid, type = "l", lty = lty, ylim = ylim,ylab = ylab)
  polygon(c(t,rev(t)),
          c(upper,rev(lower)),border = NA, col = fill_col)
}

random_Incidence_pred = function(MCMC_obj,idx,thin, col = "black", fill_col = rgb(0.15,0.15,0.15,0.3), row = 3,lty = 2,ylim = c(0,1000),ylab = "incid.predict",add = F,xaxt = "s"){
  t = SIR_incidence_Traj(MCMC_obj$Trajectory[,,1],1:dim(MCMC_obj$Trajectory)[1]-1)[,1]

  incid_MX = apply(MCMC_obj$Trajectory[,,idx/thin], 3, function(x){
    return(SIR_incidence_Traj(x,1:dim(x)[1]-1)[,2])
  })
  rhosMX = matrix(rep(MCMC_obj$incid_par[idx,1], each = length(t)), ncol = length(idx))

  incid_MX = incid_MX * rhosMX

  mid = apply(incid_MX,1,median)
  upper = apply(incid_MX,1,function(x){
    return(quantile(x,0.975))
  })
  lower = apply(incid_MX,1,function(x){
    return(quantile(x,0.025))
  })
  if(add){
    lines(t, mid, type = "l", lty = lty, ylim = ylim, ylab = ylab, lwd = 2)
  }else{
    plot(t, mid, type = "l", lty = lty, ylim = ylim,ylab = ylab,lwd=2,xaxt = xaxt)
  }
  polygon(c(t,rev(t)),
          c(upper,rev(lower)),border = NA, col = fill_col)
}

random_Incidence_predictive = function(MCMC_obj,idx,thin, col = "black", fill_col = rgb(0.15,0.15,0.15,0.3), row = 3,lty = 2,ylim = c(0,1000),ylab = "incid.predict",add = F,xaxt = "s"){
  t = SIR_incidence_Traj(MCMC_obj$Trajectory[,,1],1:dim(MCMC_obj$Trajectory)[1]-1)[,1]

  incid_MX = apply(MCMC_obj$Trajectory[,,idx/thin], 3, function(x){
    return(SIR_incidence_Traj(x,1:dim(x)[1]-1)[,2])
  })
  rhosMX = MCMC_obj$incid_par[idx,1]
  phiMX = MCMC_obj$incid_par[idx,2]

  incid_MX = incid_MX * rhosMX
  sample_MX = matrix(ncol = length(idx), nrow = length(t))
  for(i in 1:length(rhosMX)){
      sample_MX[,i] = sapply(incid_MX[,i], function(x){
      return(rnbinom(1,size = phiMX[i],mu = rhosMX[i] * x))
    })
  }


  mid = apply(sample_MX,1,function(x){
    return(quantile(x,0.5,na.rm = T))
  })
  upper = apply(sample_MX,1,function(x){
    return(quantile(x,0.975,na.rm = T))
  })
  lower = apply(sample_MX,1,function(x){
    return(quantile(x,0.025,na.rm = T))
  })
  if(add){
    lines(t, mid, type = "l", lty = lty, ylim = ylim, ylab = ylab, lwd = 2)
  }else{
    plot(t, mid, type = "l", lty = lty, ylim = ylim,ylab = ylab,lwd=2,xaxt = xaxt)
  }
  polygon(c(t,rev(t)),
          c(upper,rev(lower)),border = NA, col = fill_col)
}



CI_Curve_irate = function(MCMC_obj,idx,ids2,col_id, dt,col = "black", fill_col = rgb(0.15,0.15,0.15,0.3), row = 3,method = "qtile",lty = 2,alpha=0.05,fill = T){

  if(method == "NormApp"){
    midcur = apply(MCMC_obj$Trajectory[,row,idx]/scale,1,mean)
    midSd = apply(MCMC_obj$Trajectory[,row,idx]/scale,1,sd)
    qt = qnorm(1-alpha/2)
    upper = midcur + qt * midSd
    lower = midcur - qt * midSd
    #lines(MCMC_obj$Trajectory[,1,1],upper,lty=2,
    #      col=col,lwd=2)

    #lines(MCMC_obj$Trajectory[,1,1],lower,lty=2,
    #      col=col,lwd=2)
  }else if(method == "qtile"){
    qt1 = 1 - alpha/2
    qt2 = alpha/2
    R0 = MCMC_obj$par[idx,3] * MCMC_obj$par[idx,4] * dt
    R0_traj = matrix(ncol= length(col_id)+2, nrow = length(idx))
    R0_traj[,1] = R0
    for(i in 1:length(col_id)){
      R0_traj[,i+1] = R0_traj[,i] * MCMC_obj$par[idx,col_id[i]]
    }
    i = length(col_id)
    R0_traj[,length(col_id)+2] = R0_traj[,i] * MCMC_obj$par[idx,col_id[i]]

    upper = apply(MCMC_obj$Trajectory[,row,ids2] * t(R0_traj),1,function(x){
      return(quantile(x,qt1))
    }
    )
    lower = apply(MCMC_obj$Trajectory[,row,ids2] * t(R0_traj),1,function(x){
      return(quantile(x,qt2))
    })
    mid = apply(MCMC_obj$Trajectory[,row,ids2] * t(R0_traj),1,function(x){
      return(median(x))
    })
    # lines(MCMC_obj$Trajectory[,1,1],midcur + qt * midSd,lty=2,
    #      col=col,lwd=2)

    #lines(MCMC_obj$Trajectory[,1,1],midcur - qt * midSd,lty=2,
    #     col=col,lwd=2)
  }
  if(fill == T){
    polygon(x = c(MCMC_obj$Trajectory[,1,1],rev(MCMC_obj$Trajectory[,1,1])),
            y = c(upper,rev(lower)),col = fill_col,border = NA)
    lines(MCMC_obj$Trajectory[,1,1], mid, col = col, lwd = 2,lty = lty)
  }
}

CI_Curve_stat = function(MCMC_obj,ids,scale = 1, row = 3,alpha = 0.05){

    qt1 = 1 - alpha/2
    qt2 = alpha/2
    upper = apply(MCMC_obj$Trajectory[,row,ids]/scale,1,function(x){
      return(quantile(x,qt1))
    }
    )
    lower = apply(MCMC_obj$Trajectory[,row,ids]/scale,1,function(x){
      return(quantile(x,qt2))
    })

    mid = apply(MCMC_obj$Trajectory[,row,ids]/scale,1,median)
    return(list(upper = upper, lower = lower, mid = mid))
}

R0_Curve_stat = function(times, MCMC_obj, R0_id, col_id, idx){
  R0 = MCMC_obj$par[idx,R0_id]
  R0_traj = matrix(ncol= length(col_id)+2, nrow = length(idx))
  R0_traj[,1] = R0
  for(i in 1:length(col_id)){
    R0_traj[,i+1] = R0_traj[,i] * MCMC_obj$par[idx,col_id[i]]
  }
  i = length(col_id)
  R0_traj[,length(col_id)+2] = R0_traj[,i] * MCMC_obj$par[idx,col_id[i]]
  CIup = apply(R0_traj,2,function(x){
    return(quantile(x,0.975))
  })
  CIlow = apply(R0_traj,2,function(x){
    return(quantile(x,0.025))
  })
  m = apply(R0_traj,2,median)
  return(list(times = times, up = CIup, low = CIlow, mid = m))
}



R0_compare = function(RCI, param, x_r, x_i){
  RR = RTs(param, RCI$times, x_r, x_i)
  MAD = mean(abs(RCI$mid - RR))
  MCIW = mean(RCI$up - RCI$low)
  return(list(MAD = MAD, MCIW = MCIW))
}

random_Incidence_pred_stat = function(MCMC_obj, idx, thin, incid = T){
  t = SIR_incidence_Traj(MCMC_obj$Trajectory[,,1],1:dim(MCMC_obj$Trajectory)[1]-1)[,1]

  incid_MX = apply(MCMC_obj$Trajectory[,,idx/thin], 3, function(x){
    return(SIR_incidence_Traj(x,1:dim(x)[1]-1)[,2])
  })
  if(incid == T){
    rhosMX = MCMC_obj$incid_par[idx,1]
    phiMX = MCMC_obj$incid_par[idx,2]
  }else{
    rho_pr = MCMC_obj$MCMC_setting$prior$rho_pr
    phi_pr = MCMC_obj$MCMC_setting$prior$phi_pr
    rhosMX = sigmoid(rnorm(length(idx),rho_pr[1], rho_pr[2]))
    phiMX =  rlnorm(length(idx),phi_pr[1], phi_pr[2])
  }
  incid_MX = incid_MX * rhosMX
  incid_MX = pmax(incid_MX, 0)
  sample_MX = matrix(ncol = length(idx), nrow = length(t))
  for(i in 1:length(rhosMX)){
    sample_MX[,i] = sapply(incid_MX[,i], function(x){
      return(rnbinom(1,size = phiMX[i],mu = rhosMX[i] * x))
    })
  }

  mid = apply(sample_MX,1,function(x){
    return(quantile(x,0.5,na.rm = T))
  })
  upper = apply(sample_MX,1,function(x){
    return(quantile(x,0.975,na.rm = T))
  })
  lower = apply(sample_MX,1,function(x){
    return(quantile(x,0.025,na.rm = T))
  })

  return(list(times = t, up = upper, low = lower, mid = mid))
}

MCMC_summarize = function(times, MCMC_obj, idx, thin, alpha = 0.05){

  I0 = MCMC_obj$par[idx, 2]
  R0 = MCMC_obj$par[idx, 3]
  gamma = MCMC_obj$par[idx, 4]

  p = dim(MCMC_obj$par)[2]
  pb = c(0.025, 0.5, 0.975)
  parDensity = list(I0 = density(I0), R0 = density(R0), gamma = density(gamma),I0q = quantile(I0, pb), R0q = quantile(R0, pb), gammaq = quantile(gamma,pb))
  R0_traj = R0_Curve_stat(times, MCMC_obj, 3, 5:(p - 1), idx)
  I_traj = CI_Curve_stat(MCMC_obj, idx/thin, 1, 3, alpha)
  S_traj = CI_Curve_stat(MCMC_obj, idx/thin, 1, 2, alpha)
  #incid_pred = random_Incidence_pred_stat(MCMC_obj, idx, thin, F)
  return(list(parDensity = parDensity, R0_traj = R0_traj, S_traj = S_traj, I_traj = I_traj, niter = max(idx)))
}




MCMC_summarize_incid = function(times, MCMC_obj, idx, thin, alpha = 0.05, incid = T){

  I0 = MCMC_obj$par[idx, 2]
  R0 = MCMC_obj$par[idx, 3]
  gamma = MCMC_obj$par[idx, 4]

  p = dim(MCMC_obj$par)[2]
  pb = c(0.025, 0.5, 0.975)
  parDensity = list(I0 = density(I0), R0 = density(R0), gamma = density(gamma),I0q = quantile(I0, pb), R0q = quantile(R0, pb), gammaq = quantile(gamma,pb))
  R0_traj = R0_Curve_stat(times, MCMC_obj, 3, 5:(p - 1), idx)
  S_traj = CI_Curve_stat(MCMC_obj, idx/thin, 1, 2, alpha)
  I_traj = CI_Curve_stat(MCMC_obj, idx/thin, 1, 3, alpha)
  incid_par = list(rho = density(MCMC_obj$incid_par[idx,1]), rhoq = quantile(MCMC_obj$incid_par[idx,1], pb), phi = density(MCMC_obj$incid_par[idx,2]), phiq = quantile(MCMC_obj$incid_par[idx,2], pb))
  incid_pred = random_Incidence_pred_stat(MCMC_obj, idx, thin, incid)

  return(list(parDensity = parDensity, R0_traj = R0_traj, I_traj = I_traj,S_traj = S_traj, incid_par = incid_par, incid_pred = incid_pred))
}

MCMC_summarize_pref = function(times, MCMC_obj, idx, thin, alpha, pref = T){

  I0 = MCMC_obj$par[idx, 2]
  R0 = MCMC_obj$par[idx, 3]
  gamma = MCMC_obj$par[idx, 4]

  p = dim(MCMC_obj$par)[2]
  pb = c(0.025, 0.5, 0.975)
  parDensity = list(I0 = density(I0), R0 = density(R0), gamma = density(gamma),I0q = quantile(I0, pb), R0q = quantile(R0, pb), gammaq = quantile(gamma,pb))
  R0_traj = R0_Curve_stat(times, MCMC_obj, 3, 5:(p - 1), idx)
  I_traj = CI_Curve_stat(MCMC_obj, idx/thin, 1, 3, alpha)
  S_traj = CI_Curve_stat(MCMC_obj, idx/thin, 1, 2, alpha)
  Pref_par = list(a = density(MCMC_obj$pref_par[idx,1]), aq = quantile(MCMC_obj$pref_par[idx,1], pb), b = density(MCMC_obj$pref_par[idx,2]), bq = quantile(MCMC_obj$pref_par[idx,2], pb))
  return(list(parDensity = parDensity, R0_traj = R0_traj, I_traj = I_traj, S_traj = S_traj,Pref_par = Pref_par))
}



