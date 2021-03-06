
covtocor = function(MX){
  return(MX / (sqrt(diag(MX)) %*% t(sqrt(diag(MX)))) )
}

partocov = function(MCMC_obj, ids, indices, logvec, COV = 1){
  MX = NULL
  for(i in 1:length(indices)){
    if( logvec[i] == 1){
      v = log( MCMC_obj$par[ids,indices[i]])
    }else{
      v =  MCMC_obj$par[ids,indices[i]]
    }
    MX = cbind(MX,v)
  }
  if(COV == 1){
    return(cov(MX))
  }else{
    return(cor(MX))
  }
}


simuSIRS = function(theta1,theta2,theta3,S,I,R,times){
  SIRS_exact = list()
  SIRS_exact$Pre = matrix(c(1,1,0,0,1,0,0,0,1),ncol=3,byrow = T)
  SIRS_exact$Post = matrix(c(0,2,0,0,0,1,1,0,0),ncol=3,byrow = T)
  SIRS_exact$h = function(x,t,th=c(theta1=theta1,theta2=theta2,theta3 = theta3)){
    with(as.list(c(x,th)),{
      return(c(theta1*X*Y, theta2*Y,theta3 * Z))
    })
  }
  SIRS_FRM = StepFRM(SIRS_exact)
  simu_Traj = simTs(c(X = S, Y = I, Z = R), min(times), max(times), times[2] - times[1], SIRS_FRM)
  plot(times,simu_Traj[,2],type="l")
  return(matrix(cbind(times,simu_Traj),ncol = 4))
}

#SIRS_F = simTs(c(X = 9500,Y = 500,Z = 1500),0,200,0.1,SIRS_FRM)
#plot(seq(0,200,by=0.1),SIRS_F[,2])

simuSIR = function(theta1,theta2,S,I,time){
  R = 0
  SIR = list()
  SIR$Pre=matrix(c(1,0,1,1,0,0),ncol=3)
  SIR$Post = matrix(c(0,0,2,0,0,1),ncol=3)
  SIR$h = function(x,t,th=c(theta1= theta1, theta2 = theta2)){
    with(as.list(c(x,th)),{
      return(c(theta1*X*Y, theta2*Y))
    })
  }
  SIR_F = StepFRM(SIR)
  simu_Traj = simTs(c(X = S, Y = I, Z = R), min(time), max(time), time[2] - time[1], SIR_F)
  plot(time,simu_Traj[,2],type="l")
  return(matrix(cbind(time,simu_Traj),ncol = 4))
}

simuSIRt = function(state, time, param1, x_r1, x_i1){
  R = 0
  SIR = list()
  SIR$Pre=matrix(c(1,0,1,1,0,0),ncol=3)
  SIR$Post = matrix(c(0,0,2,0,0,1),ncol=3)
  SIR$h = function(x, t, param = param1, x_r = x_r1, x_i = x_i1){
    beta1 = betaTs(param, t, x_r, x_i)
    return(c(beta1 * x[1] * x[2], param[2] * x[2]))
  }
  SIR_F = StepFRM(SIR)
  simu_Traj = simTs(c(state,0), min(time), max(time), time[2] - time[1], SIR_F)
  #plot(time,simu_Traj[,2],type="l")
  return(matrix(cbind(time, simu_Traj), ncol = 4))
}


simuSIRt_AR = function(state, time, param1, x_r1, x_i1,upper_bd = NULL, lower_bd = NULL){
  success = 0

  SIR = list()
  SIR$Pre=matrix(c(1,0,1,1,0,0),ncol=3)
  SIR$Post = matrix(c(0,0,2,0,0,1),ncol=3)
  SIR$h = function(x, t, param = param1, x_r = x_r1, x_i = x_i1){
    beta1 = betaTs(param, t, x_r, x_i)
    return(c(beta1 * x[1] * x[2], param[2] * x[2]))
  }

  SIR_F = StepFRM(SIR)

  while(success == 0){
    R = 0
    success = tryCatch({
    simu_Traj = simTs(c(state,0), min(time), max(time), time[2] - time[1], SIR_F)
      1
    }, error = function(cond){
#      message(cond)
      # Choose a return value in case of error
      return(0)
    })
    if(!is.null(upper_bd) && success == 1){
      success = (max(simu_Traj[,2]) <= upper_bd)
    }
    if(!is.null(lower_bd) && success == 1){
      success = (max(simu_Traj[,2]) >= lower_bd)
    }
  }
  return(matrix(cbind(time, simu_Traj), ncol = 4))
}



simuSEIR = function(theta1,theta2,theta3,S,E,I,time){
  R = 0
  SIR = list()
  SIR$Pre=matrix(c(1,0,1,0,1,0,0,0,1),byrow = T,ncol=3)
  SIR$Post = matrix(c(0,1,1,0,0,1,0,0,0),byrow = T,ncol=3)
  SIR$h = function(x,t,th=c(theta1= theta1, theta2 = theta2, theta3 = theta3)){
    with(as.list(c(x,th)),{
      return(c(theta1*X*Z, theta2*Y,theta3 * Z))
    })
  }
  SIR_F = StepFRM(SIR)
  simu_Traj = simTs(c(X = S, Y = E, Z = I), min(time), max(time), time[2] - time[1], SIR_F)
  #plot(time,simu_Traj[,2],type="l")
  return(matrix(cbind(time,simu_Traj),ncol = 4))
}



simuSEIRt = function(state, time, param1, x_r1, x_i1){
  R = 0
  SEIR = list()
  SEIR$Pre=matrix(c(1,0,1,0,1,0,0,0,1),byrow = T,ncol=3)
  SEIR$Post = matrix(c(0,1,1,0,0,1,0,0,0),byrow = T,ncol=3)
  SEIR$h = function(x, t, param = param1, x_r = x_r1, x_i = x_i1){
    beta1 = betaTs(param, t, x_r, x_i)
    return(c(beta1 * x[1] * x[3], param[2] * x[2], param[3] * x[3]))
  }
  SEIR_F = StepFRM(SEIR)
  simu_Traj = simTs(state, min(time), max(time), time[2] - time[1], SEIR_F)
  #plot(time,simu_Traj[,2],type="l")
  return(matrix(cbind(time,simu_Traj),ncol = 4))
}


#' Plot the posterior
#'
#' @param timegrid the time grid for the trajectory
#' @param trMatrix A matrix of trjactory for one population. The column corresponds to each posterior sample, the row corresponds to each time point
#' @param color The
#'
#'
random_trajectory = function(timegrid,trMatrix,color="grey",ylab="infected population",main="",lwd=0.5,ylim = NULL){
  k = dim(trMatrix)[2]
  if(is.null(ylim)){
  plot(timegrid,trMatrix[,1],lwd=lwd,col=color,ylab=ylab,
       ylim=c(min(trMatrix[,1],na.rm = T)-500,max(trMatrix[,1],na.rm = T)+500),
       main=main,type="l")
  }else{
    plot(timegrid,trMatrix[,1],lwd=lwd,col=color,ylab=ylab,
         ylim = ylim,
         main=main,type="l")
  }
  for(i in 2:k){
    lines(timegrid,trMatrix[,i],lwd=lwd,col=color)
  }
}



random_trajectory_line = function(timegrid,trMatrix,color="yellow",lwd=0.5, ylim = NULL){
 # timegrid = seq(t0,t,by=gridsize)
  k = dim(trMatrix)[2]
  if(is.null(ylim)){

    lines(timegrid,trMatrix[,1],lwd=lwd,col=color,
          ylim=c(min(trMatrix[,1])-500,max(trMatrix[,1])+500),
          main="population vs time",type="l")
  }else{
    lines(timegrid,trMatrix[,1],lwd=lwd,col=color,
          ylim=c(min(trMatrix[,1])-500,max(trMatrix[,1])+500),
          main="population vs time",type="l")

  }
  for(i in 2:k){
    lines(timegrid,trMatrix[,i],lwd=lwd,col=color)
  }
}


random_trajectory3D = function(t,gridsize,trCube){
  dev.off()
  par(mfrom(2,2))
  timegrid = seq(0,t,by=gridsize)
  k = dim(trCube)[3]
  picName = c("Susceptible","infectious","removed")
  for(j in 1:3){
    random_trajectory(t,gridsize,trCube[,i,])
    lines(timegrid,apply(trCube[,i,],1,mean),col="blue",lty=2,lwd=2)
    legend("topright",legend=c("samplepath",paste("mean ",picName[i])),col=c("grey","blue"))
  }
}

Convert_R0 = function(MCMC_obj,beta_row_id = 3, gamma_row_id=4,sample_id,S0){
  return(MCMC_obj$par[sample_id,beta_row_id]/MCMC_obj$par[sample_id,gamma_row_id]*S0)
}


Covert_to_incidence_SIR = function(Traj){
  n = dim(Traj)[1]
  Initial = Traj[1,2:3]
  result = matrix(nrow = n, ncol = 3)
  result[,1] = Traj[,1]
  result[,2:3] = c(0,0)
  a = apply(Traj[2:n,2:3], 1, function(x){
    r1 = Initial[1] - x[1]
    r2 = Initial[2] - x[2] + r1
    return(c(r1,r2))
  })
  result[2:n,2:3] = t(a)
  inter_res = cbind(result[1:(n-1),1], diff(result[,2]),diff(result[,3]))
  return(list(cum = result, count = inter_res))
}
