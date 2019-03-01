/*
#include "basic_util.h"
#include "SIR_phylodyn.h"
#define pi 3.14159265358979323846264338327950288

typedef arma::vec (*parat)(double, arma::vec, arma::vec, arma::ivec);
typedef arma::vec (*ODE_fun)(arma::vec, arma::vec, double, arma::vec, arma::ivec, std::string, std::string);
typedef arma::mat (*F_fun)(arma::vec, arma::vec, std::string);
typedef arma::vec (*h_fun)(arma::vec, arma::vec, std::string);
*/

#include<RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include<Rcpp.h>
#include<R.h>
#include<math.h>
#include<stdexcept>
#include<assert.h>
using namespace Rcpp;
using namespace arma;

#define pi 3.14159265358979323846264338327950288

typedef arma::vec (*parat)(double, arma::vec, arma::vec, arma::ivec);
typedef arma::vec (*ODE_fun)(arma::vec, arma::vec, double, arma::vec, arma::ivec, std::string, std::string);
typedef arma::mat (*F_fun)(arma::vec, arma::vec, std::string);
typedef arma::vec (*h_fun)(arma::vec, arma::vec, std::string);



arma::mat chols(arma::mat S){
  arma::mat res(2,2);
  res(0,0) = sqrt(S(0,0));
  res(0,1) = S(1,0) / res(0,0);
  res(1,0) = 0;
  res(1,1) = sqrt(S(1,1) - res(0,1) * res(0,1));
  return res.t();
}

//[[Rcpp::export()]]
arma::mat mvrnormArma(int n, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = randn(ncols,1);
  // arma::mat res = arma::chol(sigma+0.00000000001 * arma::diagmat(ones(3))) * Y;
  arma::mat res;
  if(n == 2){
    res = chols(sigma) * Y;
  }else{
    res = arma::chol(sigma + 0.0000000000001 * arma::diagmat(ones(n)) ).t() * Y;
  }
  if(res.has_nan()){
    Rcout<<"666"<<sigma<<endl;
  }
  return res;
}


//[[Rcpp::export()]]
arma::mat LogTraj(arma::mat Traj){
  int p = Traj.n_cols - 1;
  arma::mat logTraj(Traj.n_rows, Traj.n_cols);
  logTraj.col(0) = Traj.col(0);
  logTraj.cols(1,p) = log(Traj.cols(1,p));
  return logTraj;
}

/*
 *
 * Generate a vector of beta corresponds to each time grid
 *
 * param: parameters for the model
 * index: index vector of length two:
 * first element is the index for infection rate,
 * 2nd element is the index for recovery rate
 *
 * time grid
 * x_r: parameters for change points
 * x_i: index for change points
 */
//[[Rcpp::export()]]
arma::vec betaTs(arma::vec param, arma::vec times, arma::vec x_r, arma::ivec x_i){
  /*
   * x_r = (N,cht1, cht2, ...)
   * x_i = (nch,nparam, index0, index1)
   *
   * return a vector of infection rate beta for different time
   *
   * example:
   * betaTs(c(1.2,40,33,0.6), seq(0,1.5,0.0001), c(1000000,0.5),c(1,3,0,2))
   */
  double R0 = param[x_i(2)];
  int i = 0;
  int nch = x_i[0];
  int m = times.n_elem;
  arma::vec betas(m);

  for(int j = 0; j < m; j++){
    if(i < nch){
      while(x_r[i + 1] <= times[j]){
        R0 *= param[x_i[1] + i];
        i ++;
        if(i == nch) break;
      }
    }
    betas(j) = R0 * param[x_i(3)] / x_r[0];
  }
  return betas;
}


//[[Rcpp::export()]]
arma::vec RTs(arma::vec param, arma::vec times, arma::vec x_r, arma::ivec x_i){
  /*
  * x_r = (N,cht1, cht2, ...)
  * x_i = (nch,nparam, index0, index1)
  *
  * return a vector of infection rate beta for different time
  *
  * example:
  * betaTs(c(1.2,40,33,0.6), seq(0,1.5,0.0001), c(1000000,0.5),c(1,3,0,2))
  */
  double R0 = param[x_i(2)];
  int i = 0;
  int nch = x_i[0];
  int m = times.n_elem;
  arma::vec betas(m);

  for(int j = 0; j < m; j++){
    if(i < nch){
      while(x_r[i + 1] <= times[j]){
        R0 *= param[x_i[1] + i];
        i ++;
        if(i == nch) break;
      }
    }
    betas(j) = R0;
  }
  return betas;
}



//[[Rcpp::export()]]
arma::vec param_transform(double t, arma::vec param, arma::vec x_r, arma::ivec x_i){
  /**
  * x_r = (N,cht1, cht2, ...)
  *
  * x_i = (nch,nparam, index0, index1)
  * theta = (beta,gamma,..., ch1,ch2,...) at time t
  *
  * example:
  * param_transform(0.7,c(1.2,40,33,0.6),c(1000000,0.5),c(1,3,0,2))
  */
  arma::vec param2 = param;
  double R0 = param[x_i(2)];
  int i = 0;
  int nch = x_i[0];
  if(nch > 0){
    while(x_r[i + 1] <= t){
      R0 *= param[x_i[1] + i];
      if(i == nch - 1) break;
      i ++;
    }
  }
  param2[0] = R0 * param[x_i(3)] / x_r[0];
  return param2;
}


/*
 * pointer for transform parameter
 * input trans: transform stype
 * {changepoint, periodic, ... }
 * return a function pointer of transform paramters into ODE parameters
 */

XPtr<parat> transformPtr(std::string trans = "changepoint"){
  if(trans == "changepoint"){
    return(XPtr<parat>(new parat(&param_transform)));
  }else{
    return XPtr<parat>(R_NilValue); // runtime error as NULL no XPtr
  }
}



//[[Rcpp::export()]]
arma::vec ODE_SIR_one(arma::vec states, arma::vec param, double t, arma::vec x_r, arma::ivec x_i,
                      std::string transP = "changepoint", std::string transX = "standard"){
  double dx, dy;
  arma::vec res(2);
  XPtr<parat> param_trans = transformPtr(transP);
  //double th1 = exp(param[0] + param[3] * sin(2 * pi * t / 40.0));
  arma::vec thetas = (*param_trans)(t,param,x_r,x_i);
  if(transX == "standard"){
    dx = - thetas[0] * states[0] * states[1];
    dy = thetas[0] * states[0] * states[1] - thetas[1] * states[1];
  }else if(transX == "log"){
    dx = -thetas[0] * exp(states[1]) - thetas[0] * exp(states[1]-states[0])/2;
    dy = thetas[0] * exp(states[0]) - thetas[1] - thetas[0] * exp(states[0]-states[1])/2 - thetas[1] *exp(-states[1])/2;
  }else{
    dx = 0;
    dy = 0;
  }
  res(0) = dx;
  res(1) = dy;
  return res;
}

//[[]Rcpp::export()]]
arma::vec ODE_SEIR2_one(arma::vec states, arma::vec param, double t,
                        arma::vec x_r, arma::ivec x_i,
                        std::string transP = "changepoint",
                        std::string transX = "standard"){
  /*
   * Approximate the susceptible population with true population
   *
   * state = {E, I}
   * parameters = {beta, mu, gamma}
   *
   * return one step result for SEIR model
   *
   */
  double dy, dz;
  arma::vec res(2);
  XPtr<parat> param_trans = transformPtr(transP);
  arma::vec thetas = (*param_trans)(t,param,x_r,x_i);
  double N = x_r[0];
  dy = thetas[0] * N * states[1] - thetas[1] * states[0];
  dz = thetas[1] * states[0] - thetas[2] * states[1];
  res(0) = dy;
  res(1) = dz;
  return res;
}

//[[Rcpp::export()]]
arma::vec ODE_SEIR_one(arma::vec states, arma::vec param, double t,
                       arma::vec x_r, arma::ivec x_i,
                      std::string transP = "changepoint",
                      std::string transX = "standard"){
  /*
   * states = {S,E,I}
   * params = {beta,mu,gamma}
   *
   */

  double dx, dy, dz; // Susceptible, exposed, Infected
  arma::vec res(3); // stores the result vector
  // transform parameter
  XPtr<parat> param_trans = transformPtr(transP);
  arma::vec thetas = (*param_trans)(t,param,x_r,x_i);
  if(transX == "standard"){
    // theta0: infection rate, theta1: rate from exposed to infected
    // theta2: recover rate

    dx = - thetas[0] * states[0] * states[2];
    dy = thetas[0] * states[0] * states[2] - thetas[1] * states[1];
    dz = thetas[1] * states[1] - thetas[2] * states[2];

  }else if(transX == "log"){

    dx = -thetas[0] * exp(states[2]) - thetas[0] * exp(states[2]-states[0])/2;

    dy = thetas[0] * exp(states[0] + states[2] - states[1]) - thetas[1] -
      thetas[0] * exp(states[0]+states[2] - 2 * states[1])/2 - thetas[1] *exp(-states[1])/2;

    dz = thetas[1] * exp(states[1] - states[2]) - thetas[2] -
      thetas[2] * exp(-states[2])/2 - thetas[1] * exp(states[1] - 2 * states[2])/2;

  }else{
    dx = 0;
    dy = 0;
    dz = 0;
  }
  res(0) = dx;
  res(1) = dy;
  res(2) = dz;
  return res;
}

/*
double log_Preferential(List init, arma::mat Trajectory, arma::vec prefer_param){
  arma::vec lambda = prefer_param(0) + prefer_param(1) * Trajectory.col(2);
  return sum(log(lambda) % as<vec>(init[10]) - lambda);
}
 */

//[[Rcpp::export()]]
arma::vec ODE_SIRS_one(arma::vec states,arma::vec param, double t, arma::vec x_r, arma::ivec x_i,
                       std::string transP = "changepoint", std::string transX = "standard"){
  double dx, dy;
  XPtr<parat> param_trans = transformPtr();
  //double th1 = exp(param[0] + param[3] * sin(2 * pi * t / 40.0));
  arma::vec thetas = (*param_trans)(t,param,x_r,x_i);
  dx = - thetas[0] * states[0] * states[1];
  dy = - thetas[0] * states[0] * states[1] - thetas[1] * states[1];

  arma::vec res(2);
  res(0) = dx;
  res(1) = dy;
  return res;
}






//[[Rcpp::export()]]
arma::mat SIRS_F(arma::vec states,arma::vec thetas,std::string transX){
  arma::mat M;
  //XPtr<parat> param_trans = transformPtr();
  //double th1 = exp(param[0] + param[3] * sin(2 * pi * t / 40.0));
  //arma::vec thetas = (*param_trans)(t,param,x_r,x_i);
  if(transX == "standard"){
    M << - thetas[0] * states[2] << 0 << - thetas[0] * states[0] << arma::endr
      << thetas[0] * states[2] << -thetas[1] << thetas[0] * states[0] << arma::endr
      << 0 << thetas[1] << - thetas[2] << endr;
  }else if(transX == "log"){
    M << thetas[0] * exp(states(2)-states(0))/2 << 0 << -thetas[0]*exp(states(2))-thetas[0]/2*exp(states(2)-states(0)) <<arma::endr;
     // << thetas[0] * exp(states[0] + states[2] -)
    //  <<thetas[0]*exp(states(0))-thetas[0]*exp(states(0)-states(1))/2 << thetas[0]*exp(states(0)-states(1))/2+thetas[1]/2*exp(-states(1))<<arma::endr;
  }
  return M;
}






//[[Rcpp::export()]]
arma::mat SIR_F(arma::vec states,arma::vec thetas,std::string transX){
  arma::mat M;
  //XPtr<parat> param_trans = transformPtr();
  //double th1 = exp(param[0] + param[3] * sin(2 * pi * t / 40.0));
  //arma::vec thetas = (*param_trans)(t,param,x_r,x_i);
  if(transX == "standard"){
    M << - thetas[0] * states[1] << - thetas[0] * states[0] << arma::endr
      << thetas[0] * states[1] << thetas[0] * states[0] - thetas[1] << arma::endr;
  }else if(transX == "log"){
    M << thetas[0] * exp(states(1)-states(0))/2 << -thetas[0]*exp(states(1))-thetas[0]/2*exp(states(1)-states(0)) <<arma::endr
      <<thetas[0]*exp(states(0))-thetas[0]*exp(states(0)-states(1))/2 << thetas[0]*exp(states(0)-states(1))/2+thetas[1]/2*exp(-states(1))<<arma::endr;
  }
  return M;
}


//[[Rcpp::export()]]
arma::mat SEIR2_F(arma::vec states, arma::vec thetas, std::string transX){
  arma::mat M;
  double N = 1000000;
  M << -thetas[1] << thetas[0] * N << arma::endr
    << thetas[1] << -thetas[2] << arma::endr;
  return M;
}



//[[Rcpp::export()]]
arma::mat SEIR_F(arma::vec states,arma::vec thetas,std::string transX){
  arma::mat M;

  //A << -1 << 1<<arma::endr<<0<< -1 << arma::endr;
  //XPtr<parat> param_trans = transformPtr();
  //double th1 = exp(param[0] + param[3] * sin(2 * pi * t / 40.0));
  //arma::vec thetas = (*param_trans)(t,param,x_r,x_i);
  if(transX == "standard"){
    M << - thetas[0] * states[2] << 0 << - thetas[0] * states[0] << arma::endr
      << thetas[0] * states[2] << - thetas[1] << thetas[0] * states[0] << arma::endr
      << 0 << thetas[1] << -thetas[2] << arma::endr;
  }else if(transX == "log"){
    M << thetas[0] * exp(states(2)-states(0))/2 << 0 << -thetas[0]*exp(states(2))-thetas[0]/2*exp(states(2)-states(0)) <<arma::endr
      << thetas[0] * exp(states(0) + states(2) - states(1)) - thetas[0] * exp(states(0) + states(2) - 2*states(1)) / 2
      << -thetas[0] * exp(states(0) + states(2) - states(1)) + thetas[1] * exp(-states[1])/2 + thetas[0] * exp(states(0) + states(2) - 2*states(1))
      <<  thetas[0] * exp(states(0) + states(2) - states(1)) - thetas[0] * exp(-states(0) + states(2) - 2* states(1)) << arma::endr
      << 0 << thetas[1]*exp(states(1)-states(2)) - thetas[1] * exp(states[1] - 2*states[2])/2
      <<  -thetas[1]*exp(states(1)-states(2)) + thetas[2] * exp(-states[2])/2 + thetas[1] * exp(states[1] - 2 * states[2])  << arma::endr;
  }
  return M;
}



//[[Rcpp::export()]]
arma::vec SIR_h(arma::vec states,arma::vec thetas,std::string transX = "standard"){
  arma::vec h(2);
  if(transX == "standard"){
    //XPtr<parat> param_trans = transformPtr();
    //double th1 = exp(param[0] + param[3] * sin(2 * pi * t / 40.0));
    //arma::vec thetas = (*param_trans)(tstates(0),param,x_r,x_i);
    h(0)= thetas[0] * states(0) * states(1);
    h(1) = thetas[1] * states(1);
  }else if(transX == "log"){
    h(0)= thetas[0] * exp(states(0)) * exp(states(1));
    h(1) = thetas[1] * exp(states(1));
  }
  return h;
}

//[[Rcpp::export()]]
arma::vec SEIR2_h(arma::vec states, arma::vec thetas, std::string transX = "standard"){
  arma::vec h(3);
  double N = 1000000;
  h(0) = thetas[0] * N * states(1);
  h(1) = thetas[1] * states(0);
  h(2) = thetas[2] * states(1);
  return h;
}


//[[Rcpp::export()]]
arma::vec SEIR_h(arma::vec states,arma::vec thetas,std::string transX = "standard"){
  arma::vec h(3);
  if(transX == "standard"){
    //XPtr<parat> param_trans = transformPtr();
    //double th1 = exp(param[0] + param[3] * sin(2 * pi * t / 40.0));
    //arma::vec thetas = (*param_trans)(tstates(0),param,x_r,x_i);

    h(0)= thetas[0] * states(0) * states(2);
    h(1) = thetas[1] * states(1);
    h(2) = thetas[2] * states(2);
  }else if(transX == "log"){
    h(0) = thetas[0] * exp(states(0)) * exp(states(1));
    h(1) = thetas[1] * exp(states(1));
    h(2) = thetas[2] * exp(states(2));
  }
  return h;
}



/*
//[[Rcpp::export()]]
arma::vec ODE_SIR_log_one(arma::vec states,arma::vec param, double t, arma::vec x_r, arma::ivec x_i,
                           std::string transP = "changepoint"){
  double dx, dy;
  XPtr<parat> param_trans = transformPtr(transP);
  //double th1 = exp(param[0] + param[3] * sin(2 * pi * t / 40.0));
  arma::vec thetas = (*param_trans)(t,param,x_r,x_i);

  dx = - thetas[0] * states[0] * states[1];
  dy = - thetas[0] * states[0] * states[1] - thetas[1] * states[1];

  arma::vec res(2);
  res(0) = dx;
  res(1) = dy;
  return res;
}
arma::mat F_SIR_log(arma::vec states, arma::vec transParam){
  arma::mat M;
  double X = states[0],Y = states[1];
  M << transParam[0] * exp(Y-X)/2 << -transParam[0]*exp(Y)-transParam[0]/2*exp(Y-X) <<arma::endr
    <<transParam[0]*exp(X)-transParam[0]*exp(X-Y)/2 << transParam[0]*exp(X-Y)/2+transParam[1]/2*exp(-Y)<<arma::endr;
  return M;
}
*/

/*
 * function pointers
 *
 *
 */


/*
 * Pionter of Jacobian matrix in LNA for different EPI model
 *  input model: name of the EPI model
 *  SIR SEIS
 *
 *  return a function pointer of Jacobian matrix F (standard or with log transform)
 *
 */

XPtr<F_fun> F_funPtr(std::string model = "SIR"){
  if(model == "SIR" ){
    return XPtr<F_fun>(new F_fun(&SIR_F));
  }else if(model == "SEIR"){
    return XPtr<F_fun>(new F_fun(&SEIR_F));
  }else if(model == "SEIR2"){
    return XPtr<F_fun>(new F_fun(&SEIR2_F));
  }else{
    return XPtr<F_fun>(R_NilValue); // runtime error as NULL no XPtr
  }
}




/*
 * Pionter of rate vector h for different EPI model
 *  input model: name of the EPI model
 *  SIR SEIS
 *
 *  return a function pointer of rate vector h (standard or with log transform)
 *
 */

XPtr<h_fun> h_fun_Ptr(std::string model = "SIR"){
  if(model == "SIR"){
    return XPtr<h_fun>(new h_fun(&SIR_h));
  }else if(model == "SEIR"){
    return XPtr<h_fun>(new h_fun(&SEIR_h));
  }else if(model == "SEIR2"){
    return XPtr<h_fun>(new h_fun(&SEIR2_h));
  }else{
    return XPtr<h_fun>(R_NilValue); // runtime error as NULL no XPtr
  }
}

/**
 * Function pointer for one-step ODE integrator
 *
 */

XPtr<ODE_fun> ODE_fun_Ptr(std::string model = "SIR"){
  if(model == "SIR"){
    return XPtr<ODE_fun>(new ODE_fun(&ODE_SIR_one));
  }else if(model == "SEIR"){
    return XPtr<ODE_fun>(new ODE_fun(&ODE_SEIR_one));
  }else if(model == "SIRS"){
    return XPtr<ODE_fun>(new ODE_fun(&ODE_SIRS_one));
  }else if(model == "SEIR2"){
    return XPtr<ODE_fun>(new ODE_fun(&ODE_SEIR2_one));
  }else{
    return XPtr<ODE_fun>(R_NilValue); // runtime error as NULL no XPtr
  }
}


/*
 * general functions
 *
 *
 */



//[[Rcpp::export()]]
arma::mat ODE_rk45(arma::vec initial, arma::vec t, arma::vec param,
                           arma::vec x_r, arma::ivec x_i,
                           std::string transP = "changepoint",std::string model = "SIR",
                           std::string transX = "standard"){
  // XPtr<ODEfuncPtr> SIR_ODEfun = putFunPtrInXPtr(funname);
  //double N = initial[0] + initial[1] + initial[2];

  /*
   * example:
   *
    traj0 = ODE_rk45(c(10,10),seq(0,1.5,0.01),c(1.5,40,33,0.6),c(975000,0.6),c(1,3,0,2),model = "SEIR2")
    traj = ODE_rk45(c(10,10),seq(0,1.5,0.01),c(1.5,40,33,0.6),c(1000000,0.6),c(1,3,0,2),model = "SEIR2")
    traj2 = ODE_rk45(c(1000000,10,10),seq(0,1.5,0.01),c(1.5,40,33,0.6),c(1000000,0.6),c(1,3,0,2),model = "SEIR")
    plot(traj[,1],traj[,2],col = "red")
    lines(traj0[,1],traj0[,2],col = "blue")
    lines(traj2[,1],traj2[,3], col = "green")
   *
   */
  int n = t.n_rows;
  int p = initial.size();
  double dt = t[1] - t[0];
  arma::mat OdeTraj(n,p + 1);
  OdeTraj.col(0) = t;
  OdeTraj.submat(0,1,0,p) = initial.t();
  arma::vec X0 = initial, k1=initial, k2=initial, k3=initial, k4=initial,X1=initial;
  XPtr<ODE_fun> ODE_one = ODE_fun_Ptr(model);
  for(int i = 1; i < n; i++){
    X0 = X1;
    k1 = (*ODE_one)(X0,param,t[i-1], x_r,x_i,transP, transX);
    k2 = (*ODE_one)(X0 + k1 * dt / 2,param, t[i-1], x_r,x_i,transP, transX);
    k3 = (*ODE_one)(X0 + k2 * dt / 2,param, t[i-1], x_r,x_i,transP, transX);
    k4 = (*ODE_one)(X0 + k3 * dt / 2,param, t[i-1], x_r,x_i,transP, transX);
    X1 = X0 + (k1/6 + k2/3 + k3/3 + k4/6) * dt;
    OdeTraj.submat(i,1,i,p) = X1.t();
  }
  return OdeTraj;
}


//[[Rcpp::export()]]
double ODE_rk45_stop(arma::vec initial, arma::vec t, arma::vec param,
                   arma::vec x_r, arma::ivec x_i,
                   std::string transP = "changepoint",std::string model = "SIR",
                   std::string transX = "standard", double tol = 5){
  // XPtr<ODEfuncPtr> SIR_ODEfun = putFunPtrInXPtr(funname);
  //double N = initial[0] + initial[1] + initial[2];

  /*
  * example:
  *
  traj0 = ODE_rk45(c(10,10),seq(0,1.5,0.01),c(1.5,40,33,0.6),c(975000,0.6),c(1,3,0,2),model = "SEIR2")
  traj = ODE_rk45(c(10,10),seq(0,1.5,0.01),c(1.5,40,33,0.6),c(1000000,0.6),c(1,3,0,2),model = "SEIR2")
  traj2 = ODE_rk45(c(1000000,10,10),seq(0,1.5,0.01),c(1.5,40,33,0.6),c(1000000,0.6),c(1,3,0,2),model = "SEIR")
  plot(traj[,1],traj[,2],col = "red")
  lines(traj0[,1],traj0[,2],col = "blue")
  lines(traj2[,1],traj2[,3], col = "green")
  *
  */
  int n = t.n_rows;
  double dt = t[1] - t[0];
  arma::vec X0 = initial, k1=initial, k2=initial, k3=initial, k4=initial,X1=initial;
  XPtr<ODE_fun> ODE_one = ODE_fun_Ptr(model);
  int i = 1;
  for(; i < n; i++){
    X0 = X1;
    k1 = (*ODE_one)(X0,param,t[i-1], x_r,x_i,transP, transX);
    k2 = (*ODE_one)(X0 + k1 * dt / 2,param, t[i-1], x_r,x_i,transP, transX);
    k3 = (*ODE_one)(X0 + k2 * dt / 2,param, t[i-1], x_r,x_i,transP, transX);
    k4 = (*ODE_one)(X0 + k3 * dt / 2,param, t[i-1], x_r,x_i,transP, transX);
    X1 = X0 + (k1/6 + k2/3 + k3/3 + k4/6) * dt;
    if(X1(x_i(3)) < tol){
      return t[i];
    }
  }
  return t[i-1]+0.1;
}




//[[Rcpp::export()]]
List SigmaF(arma::mat Traj_par,arma::vec param,
                  arma::vec x_r, arma::ivec x_i,
                  std::string transP = "changepoint", std::string model = "SIR",std::string transX = "standard"){
  /*
   * Input:
   * Traj_par: SIR ode trjactory
   * param: parameters for SIR model
   *
   * Return:
   * A list contains two matrices: Jocobian matrix F and Covariance matrix Sigma
   *
   */

  int p = Traj_par.n_cols - 1;

  arma::mat Q = arma::eye(p,p);
  arma::mat Sigma,F,F0 = arma::eye(p,p); // F0 be the fundemental matrix
  Sigma.zeros(p,p);
  int k = Traj_par.n_rows;
  arma::mat A;
  arma::vec h,state(p);
  if(model == "SIR"){
    A << -1 << 1  <<arma::endr<<0<< -1 << arma::endr;
  }else if(model == "SEIR"){
    A << -1 << 1 << 0 << arma::endr << 0 << -1 << 1 << arma::endr << 0 << 0 << -1 << endr;
  }else if(model == "SEIR2"){
    A << 1 << 0 << arma::endr << -1 << 1 <<arma::endr << 0 << -1 << arma::endr;
  }
  arma::mat H;
  //F0 << 1 <<0 <<arma::endr << 0 << 1 <<endr;
  double t;
  double dt = Traj_par(1,0) - Traj_par(0,0);
  XPtr<parat> param_trans = transformPtr(transP);
  XPtr<h_fun> h_Function = h_fun_Ptr(model);
  XPtr<F_fun> F_function = F_funPtr(model);
  //double th1 = exp(param[0] + param[3] * sin(2 * pi * t / 40.0));
  for(int i = 0; i < k; i++){
    t = Traj_par(i,0);
    arma::vec thetas = (*param_trans)(t,param,x_r,x_i);
    for(int j = 0; j < p; j ++){
      state(j) = Traj_par(i,j+1);
      if(transX == "log"){
        Q(j,j) = exp(-state(j));
      }
    }

    h = (*h_Function)(state, thetas, transX);
    H = diagmat(h);
    F = (*F_function)(state, thetas, transX);
    F0 = F0 + (F * F0) * dt;
    Sigma = Sigma + (Sigma * F.t() + F * Sigma +  Q * A.t() * H * A * Q ) * dt;

  }

  List Res;
  Res["expF"] = F0;
  Res["Sigma"] = Sigma;
  if(Sigma.has_nan()){
    Rcout<<"na in sigma"<<endl;
 //   Rcout<<Traj_par.row(0)<<endl;
  }
  return Res;
}


//[[Rcpp::export()]]
List KF_param(arma::mat OdeTraj, arma::vec param,int gridsize,arma::vec x_r, arma::ivec x_i,
              std::string transP = "changepoint",
              std::string model = "SIR",std::string transX = "standard"){
  /*
   * Parametrize the latent process as a Kalman filter
   * X_{i=1} = A * X_i + epsilon where E(Epsilon) = 0  Var(epsilon) = Sigma
   *
   * Return a list with two cubes
   * 1. cubes of As: A_i = A.slice(i)
   * 2. cubes of Sigmas: Sigma_i = Sigma.slices(i)
   *
   */
  int n = OdeTraj.n_rows;
  //double dt = (OdeTraj(1,0) - OdeTraj(0,0));
  int k = (n-1) / gridsize;
  //int p = OdeTraj.n_cols - 1;
  int p = OdeTraj.n_cols - 1;
  arma::cube Acube(p,p,k);
  arma::cube Scube(p,p,k);
  arma::mat Traj_part;
  for(int i=0;i<k;i++){
    Traj_part = OdeTraj.submat(i*gridsize,0,(i+1)*gridsize-1,p);
    List tempres = SigmaF(Traj_part,param,x_r,x_i,transP, model, transX);
    // Rcout<< tempres<<endl;
    Acube.slice(i) = as<arma::mat>(tempres[0]);
    Scube.slice(i) = as<arma::mat>(tempres[1]);
  }
  List Res;
  Res["A"] = Acube;
  Res["Sigma"] = Scube;
  return Res;
}


//[[Rcpp::export()]]
List KF_param_chol(arma::mat OdeTraj, arma::vec param,int gridsize,arma::vec x_r, arma::ivec x_i,
              std::string transP = "changepoint",
              std::string model = "SIR",std::string transX = "standard"){
  /*
   * Using non-centralized reparametrization:
   * LatentTraj = OdeTraj + COV_matrix * iid gaussian
   * return the components of matrix used in non-centralized reparametrization
   *
   */
  int n = OdeTraj.n_rows;
  //double dt = (OdeTraj(1,0) - OdeTraj(0,0));
  int k = (n-1) / gridsize;
  //int p = OdeTraj.n_cols - 1;
  int p = OdeTraj.n_cols - 1;
  arma::cube Acube(p,p,k);
  arma::cube Lcube(p,p,k);
  arma::mat Traj_part;
  for(int i=0;i<k;i++){
    Traj_part = OdeTraj.submat(i*gridsize,0,(i+1)*gridsize-1,p);
    List tempres = SigmaF(Traj_part,param,x_r,x_i,transP, model, transX);
    // Rcout<< tempres<<endl;
    Acube.slice(i) = as<arma::mat>(tempres[0]);
    try{
      Lcube.slice(i) = arma::chol(as<arma::mat>(tempres[1]) + 0.00000001 * arma::diagmat(ones(p))).t();
    }catch(...){
      throw std::invalid_argument("Invalid input for cholesky decomposition., parametrization fails");
    }
  }
  List Res;
  Res["A"] = Acube;
  Res["Sigma"] = Lcube;
  return Res;
}
/*
List KF_param_chol(arma::mat OdeTraj, arma::vec param,arma::ivec gridsizes,arma::vec x_r, arma::ivec x_i,
                   std::string transP = "changepoint",
                   std::string model = "SIR",std::string transX = "standard"){

  int n = OdeTraj.n_rows;
  //double dt = (OdeTraj(1,0) - OdeTraj(0,0));
  int k = gridsizes.n_elem;
  //int p = OdeTraj.n_cols - 1;
  int p = OdeTraj.n_cols - 1;
  arma::cube Acube(p,p,k);
  arma::cube Lcube(p,p,k);
  arma::mat Traj_part;
  int begin = 0;
  for(int i=0;i<k;i++){
    Traj_part = OdeTraj.submat(begin,0,begin + gridsizes(i) - 1,p);
    List tempres = SigmaF(Traj_part,param,x_r,x_i,transP, model, transX);
    // Rcout<< tempres<<endl;
    Acube.slice(i) = as<arma::mat>(tempres[0]);
    if(i < k - 1) {
      begin += gridsizes(i);
    }
    try{
      Lcube.slice(i) = arma::chol(as<arma::mat>(tempres[1]) + 0.00000001 * arma::diagmat(ones(p))).t();
    }catch(...){
      throw std::invalid_argument("Invalid input for cholesky decomposition., parametrization fails");
    }
  }
  List Res;
  Res["A"] = Acube;
  Res["Sigma"] = Lcube;
  return Res;
}

*/

//[[Rcpp::export()]]
double log_like_traj_general2(arma::mat SdeTraj,arma::mat OdeTraj, List Filter,
                             int gridsize,double t_correct){
  /*
   * SdeTraj: input any random trajectory with proper dimension
   *
   *
   *
   *
   */
  assert(SdeTraj.n_cols == OdeTraj.n_cols);
  arma::cube Acube = as<arma::cube>(Filter[0]);
  arma::cube Scube = as<arma::cube>(Filter[1]);

  int k = SdeTraj.n_rows - 1;
  int p = SdeTraj.n_cols - 1;
  double loglik = 0;
  //  int id1,id0;
  arma::vec Xd1, Xd0;

  Xd1 = (SdeTraj.submat(0,1,0,p)-OdeTraj.submat(0,1,0,p)).t();
  // Xd0 = Xd1 = (SdeTraj.submat(0,1,0,3) - OdeTraj.submat(0,1,0,3)).t();
  for(int i = 0; i < k; i++){
    Xd0 = Xd1;
    //  id0 = i * gridsize;
    //  id1 = (i+1) * gridsize - 1;


    arma::mat A = Acube.slice(i);

    Xd1 = (SdeTraj.submat((i+1),1,(i+1),p)-OdeTraj.submat(i+1,1,i+1,p)).t();

    //arma::vec SigInv = arma::solve(Scube.slice(i),);

    if(SdeTraj(i+1,0) <= t_correct){
      arma::mat INexp = (Xd1 - A * Xd0).t() * solve(Scube.slice(i), Xd1 - A * Xd0);
      // Rcout <<  ((Xd1 - A * Xd0).t() * SigInv * (Xd1 - A * Xd0)) <<endl;
      loglik += -log(arma::det(Scube.slice(i)))/2.0 - 0.5 * INexp(0,0);
    }
    /*
    arma::mat A = Acube.slice(i);
    arma::mat SigInv = inv2(Scube.slice(i)+0.00001 * eye(3,3));
    arma::mat INexp = ((Xd1 - A * Xd0).t() * SigInv * (Xd1 - A * Xd0));
    loglik += log(arma::det(SigInv)) * 1.5 - 0.5 * INexp(0,0);
    */
  }

  return loglik;
}

//[[Rcpp::export()]]
double log_like_traj_general_ez(arma::mat SdeTraj, double t_correct, arma::vec initial, arma::vec t, arma::vec param,
                                arma::vec x_r, arma::ivec x_i,
                                std::string transP = "changepoint",std::string model = "SIR",
                                std::string transX = "standard"){

  arma::mat OdeTraj = ODE_rk45(initial,t,param,
                               x_r,  x_i, transP , model,transX);

  int gridsize = (int) (SdeTraj(1,0) - SdeTraj(0,0) )/ (t[1]-t[0]);
  List Filter = KF_param(OdeTraj, param, gridsize, x_r, x_i,
                         transP, model, transX);
  arma::mat OdeTraj_Coarse(SdeTraj.n_rows, SdeTraj.n_cols);
  for(int i = 0; i < SdeTraj.n_rows; i ++){
    OdeTraj_Coarse.row(i) = OdeTraj_Coarse.row(i * gridsize);
  }

  return log_like_traj_general2(SdeTraj, OdeTraj_Coarse, Filter, gridsize, t_correct);
}


//[[Rcpp::export()]]
double log_like_traj_general_adjust(arma::mat SdeTraj,arma::mat OdeTraj, List Filter_NC,
                             int gridsize,double t_correct){
/*
  arma::cube Acube = as<arma::cube>(Filter_NC[0]);
  arma::cube Lcube = as<arma::cube>(Filter_NC[1]);

  int k = SdeTraj.n_rows - 1;
  double loglik = 0;
  //  int id1,id0;

  // Xd0 = Xd1 = (SdeTraj.submat(0,1,0,3) - OdeTraj.submat(0,1,0,3)).t();
  for(int i = 0; i < k; i++){

    if(SdeTraj(i+1,0) <= t_correct){
      loglik += -log(arma::det(Lcube.slice(i)));
    }
  }

  return loglik;
*/
return 0;
}


//[[Rcpp::export()]]
List Traj_sim_general3(arma::mat OdeTraj, List Filter,double t_correct){
  arma::cube Acube = as<arma::cube>(Filter[0]);
  arma::cube Scube = as<arma::cube>(Filter[1]);
  int k = OdeTraj.n_rows - 1;

  int p = OdeTraj.n_cols - 1;
  arma::vec X0(p),eta0(p),X1(p),eta1(p);
  X1 = OdeTraj.submat(0,1,0,p).t();
  eta1 = X1;
  double loglike = 0;
  arma::mat LNA_traj(k+1,p+1);
  LNA_traj(0,0) = 0;

  LNA_traj.submat(0,1,0,p) = X1.t();
  //Rcout<<"test1"<<endl;

  for(int i = 0; i< k; i++){

    arma::mat Sig = Scube.slice((i));
    //   arma::mat SigInv = inv2(Scube.slice(i).submat(0,0,1,1));
    arma::mat A = Acube.slice(i);

    X0 = X1;

    eta0 = eta1;
    eta1 = OdeTraj.submat(i+1, 1, i+1, p).t();
    //    Rcout<<"test2"<<endl;
    arma::mat noise = mvrnormArma(p,Sig);

    X1 = A * (X0 - eta0) + eta1 + noise;

    if(OdeTraj(i+1,0) <= t_correct){
      arma::mat l1 = (-0.5) * noise.t() * arma::solve(Sig,noise);
      //     arma::mat l1  = (-0.5) * noise.t() * inv2(Sig + eye(3,3)*0.00001) * noise;
      //    loglike += -1.5 * log(det(Sig+eye(3,3)*0.00001)) + l1(0,0);
      loglike += -log(arma::det(Sig))/2.0 + l1(0,0);
    }
    //    Rcout<<"test3"<<endl;
    LNA_traj(i+1,0) = OdeTraj((i+1), 0);
    LNA_traj.submat(i+1,1,i+1,p) = X1.t();
  }

  List Res;
  Res["SimuTraj"] = LNA_traj;
  Res["loglike"] = loglike;
  return Res;
}





//[[Rcpp::export()]]
List Traj_sim_general_noncentral(arma::mat OdeTraj, List Filter_NC,double t_correct){
  /*
   * Simulate Epi trajecotry using non-centralized parametrization
   *
   * Return: A list
   * 1.Simulated trajectory
   * 2. The underlying normal variable
   * loglikelihood of the trajectory
   *
   *
   *
   */
  arma::cube Acube = as<arma::cube>(Filter_NC[0]);
  arma::cube Lcube = as<arma::cube>(Filter_NC[1]);

  int k = OdeTraj.n_rows - 1;
  int p = OdeTraj.n_cols - 1;
  arma::vec X0(p),X1(p);
  for(int i = 0; i < p; i ++){
    X0(i) = 0;
  }
  double loglike = 0;
  double logTraj = 0;
  arma::mat LNA_traj = OdeTraj;
  arma::mat OriginLatent(k,p);
  for(int i = 0; i< k; i++){
    arma::mat epsilons_i = arma::randn(p,1);
    OriginLatent.row(i) = epsilons_i.t();
    arma::mat L_i = Lcube.slice((i));
    //   arma::mat SigInv = inv2(Scube.slice(i).submat(0,0,1,1));
    arma::mat A = Acube.slice(i);
    X1 = A * X0 + L_i * epsilons_i;
    LNA_traj.submat(i+1,1,i+1,p) = X1.t() + LNA_traj.submat(i+1,1,i+1,p);
    X0 = X1;
    //Rcout<<noise.submat(0,0,1,0)<<endl;
    if(OdeTraj(i+1,0) <= t_correct){
      arma::mat l1 = (-0.5) * epsilons_i.t() * epsilons_i;
      loglike += -log(arma::det(L_i));
      logTraj += l1(0,0);
    }
  }

  List Res;
  Res["SimuTraj"] = LNA_traj;
  Res["OriginTraj"] = OriginLatent;
  Res["logMultiNorm"] = loglike;
  Res["logOrigin"] = logTraj;
  return Res;
}



//[[Rcpp::export()]]
arma::mat TransformTraj(arma::mat OdeTraj,arma::mat OriginLatent, List Filter_NC){
  int n = OriginLatent.n_rows;
  int p = OriginLatent.n_cols;
  arma::vec X0(p),X1(p);

   for(int i = 0; i < p; i ++){
    X0(i) = 0;
  }

  arma::cube Acube = as<arma::cube>(Filter_NC[0]);
  arma::cube Lcube = as<arma::cube>(Filter_NC[1]);
  arma::mat LNA_traj = OdeTraj;

  for(int i = 0; i< n; i++){
    arma::mat epsilons_i = OriginLatent.row(i).t();
    arma::mat L_i = Lcube.slice((i));

    arma::mat A = Acube.slice(i);
    X1 = A * X0 + L_i * epsilons_i;
    LNA_traj.submat(i+1,1,i+1,p) = X1.t() + LNA_traj.submat(i+1,1,i+1,p);
    X0 = X1;
  }
  return LNA_traj;
}





//[[Rcpp::export()]]
List Traj_sim_ezG2(arma::vec initial, arma::vec times, arma::vec param,
                  int gridsize,arma::vec x_r,arma::ivec x_i,double t_correct,
                  std::string transP = "changepoint",std::string model = "SIR",
                  std::string transX = "standard"){
  /*
   example:
   Traj_sim_ezG2(c(1000000,10,10),seq(0,1.5,0.0001),c(1.5,40,33,0.6),500,
   c(1000000,0.6),c(1,3,0,2),1.1,model = "SEIR")

   */

  //int p = initial.n_elem;
  int p = initial.n_elem;
  int k = times.n_rows / gridsize;
  arma::mat OdeTraj_thin = ODE_rk45(initial,times, param, x_r, x_i,
                                transP, model, transX);
  arma::mat OdeTraj(k+1,p+1);
  for(int i = 0; i< k + 1; i++){
    OdeTraj.submat(i,0,i,p) = OdeTraj_thin.submat(i*gridsize,0,i*gridsize,p);
  }

  List Filter = KF_param(OdeTraj_thin, param, gridsize, x_r, x_i,
                         transP, model, transX);

  arma::cube Acube = as<arma::cube>(Filter[0]);
  arma::cube Scube = as<arma::cube>(Filter[1]);
  k = OdeTraj.n_rows-1;
  //Rcout << Acube.slice(9) <<endl;
  arma::vec X0,X1 = initial.subvec(0,p-1),eta0(p),eta1 = initial.subvec(0,p-1);

  double loglike = 0;
  arma::mat LNA_traj(k+1,p+1);
  LNA_traj(0,0) = 0;

  LNA_traj.submat(0,1,0,p) = initial.subvec(0,p-1).t();
  //  Rcout<<"test1"<<endl;
  for(int i = 0; i< k; i++){

    arma::mat Sig = Scube.slice((i));
    if(Sig(0,0)<0){
      Rcout<<i<<endl;
    }
    //   arma::mat SigInv = inv2(Scube.slice(i).submat(0,0,1,1));
    arma::mat A = Acube.slice(i);

    X0 = X1;

    eta0 = eta1;
    eta1 = OdeTraj.submat(i+1, 1, i+1, p).t();
    arma::mat noise = mvrnormArma(p,Sig);

    X1 = A * (X0 - eta0) + eta1 + noise;
    //Rcout<<noise.submat(0,0,1,0)<<endl;
    if(OdeTraj((i+1), 0) <= t_correct){
      arma::mat l1 = (-0.5) * noise.t() * solve(Sig,noise);
      //     arma::mat l1  = (-0.5) * noise.t() * inv2(Sig + eye(3,3)*0.00001) * noise;
      //    loglike += -1.5 * log(det(Sig+eye(3,3)*0.00001)) + l1(0,0);
      loglike += -log(arma::det(Sig))/2 + l1(0,0);
    }
    //    Rcout<<"test3"<<endl;
    LNA_traj(i+1,0) = OdeTraj((i+1), 0);
    LNA_traj.submat(i+1,1,i+1,p) = X1.t();
  }

  List Res;
  Res["SimuTraj"] = LNA_traj;
  Res["loglike"] = loglike;
  return Res;
}



//[[Rcpp::export()]]
List Traj_sim_ezG_NC(arma::vec initial, arma::vec times, arma::vec param,
                  int gridsize,arma::vec x_r,arma::ivec x_i,double t_correct,
                  std::string transP = "changepoint",std::string model = "SIR",
                  std::string transX = "standard"){
  /*
   example:
   Traj_sim_ezG2(c(1000000,10,10),seq(0,1.5,0.0001),c(1.5,40,33,0.6),500,
   c(1000000,0.6),c(1,3,0,2),1.1,model = "SEIR")

   */

  //int p = initial.n_elem;
  int p = initial.n_elem;
  int k = times.n_rows / gridsize;
  arma::mat OdeTraj_thin = ODE_rk45(initial,times, param, x_r, x_i,
                                transP, model, transX);
  arma::mat OdeTraj(k+1,p+1);
  for(int i = 0; i< k + 1; i++){
    OdeTraj.submat(i,0,i,p) = OdeTraj_thin.submat(i*gridsize,0,i*gridsize,p);
  }

  List Filter = KF_param_chol(OdeTraj_thin, param, gridsize, x_r, x_i,
                         transP, model, transX);
  return Traj_sim_general_noncentral(OdeTraj, Filter, t_correct);
}









// Kingman's coalescent model
//[[Rcpp::export()]]
double coal_loglik3(List init, arma::mat f1, double t_correct, double lambda, int Index, std::string transX = "standard"){

  int n0 = 0;
  while(f1(n0,0) < t_correct){
    n0 ++;
  }
  //Rcout<<f1(n0,0)<<endl;
  //  Rcout<<n0<<endl;
  arma::vec f2(n0 + 1);
  if(transX == "standard"){
    for(int i = n0; i>0; i--){
      f2(n0-i) = log(f1(i,Index + 1));
    }
  }else if(transX == "log"){
    for(int i = n0; i>0; i--){
      f2(n0 - i) = f1(i, Index + 1);
    }
  }

  if(as<int>(init[9]) != f2.n_rows){
    Rcout<<"Incorrect length for f"<<endl;
  }

  arma::vec gridrep;
  gridrep = as<vec>(init[6]);
  int k = sum(as<arma::vec>(init[6]));

  arma::vec f(k);

  int start = 0;
  for(int i = 0; i< f2.n_rows; i++){
    for(int j = 0; j<gridrep(i);j++){
      f(start ++) = f2(i);
    }
  }
  arma::vec ll = -lambda * (as<vec>(init[2]) % as<vec>(init[3]) % arma::exp(-f)) +\
    as<vec>(init[4]) % (log(lambda) -f);
  return sum(ll);
}


//[[Rcpp::export()]]
double volz_loglik_nh3(List init, arma::mat f1, arma::vec betaN, double t_correct, arma::ivec index,
                       std::string transX = "standard"){

  int p = f1.n_cols - 1;
  int n2 = f1.n_rows;
  int n0 = as<int>(init[9]) - 1;
  int L = 0;
  while(f1(L,0) < t_correct){
      L ++;
      if(L == n2 - 1) break;
    }

  if(f1.submat(0,1,n0,p).min() < 0){
    return -10000000;
  }
  // Rcout<<f1 <<endl;
  arma::mat f2(n0,p);
  arma::vec betanh(n0);
  for(int i = 1; i<= n0; i++){
    f2.submat((n0-i),0,(n0-i),p - 1) = f1.submat(L-i+1,1,L-i+1,p);
    betanh(n0-i) = betaN(L-i+1);
  }
  // Rcout<<f2.n_rows<<"\t"<<as<int>(init[9])<<endl;
  if(n0 != f2.n_rows){
    Rcout<<"Incorrect length for f"<<endl;
  }

  arma::vec gridrep;
  gridrep = as<vec>(init[6]);
  int k = sum(as<arma::vec>(init[6]));

  arma::vec f(k);
  arma::vec s(k);
  arma::vec e(k);
  arma::vec betas(k);
  int start = 0;
  for(int i = 0; i < n0; i++){
    for(int j = 0; j < gridrep(i);j++){
      if(transX == "standard"){
        f(start) = log(f2(n0-i-1,index(1)));
        s(start) = log(f2(n0-i-1,index(0)));
        e(start) = log(f2(n0-i-1,1));
      }else if(transX == "log"){
        f(start) = f2(f2.n_rows-i-1,index(1));
        s(start) = f2(f2.n_rows-i-1,index(0));
        e(start) = f2(f2.n_rows-i-1,1);
      }
      betas(start) = betanh(n0-i-1);
      start ++;
    }
  }
  arma::vec ll = -2 * (as<vec>(init[2]) % as<vec>(init[3]) % betas % (arma::exp(- f)) % (arma::exp(s))) +\
    as<vec>(init[4]) % (log(betas) -f + s);
  //Rcout<< ll <<endl;
  return sum(ll);
}

//[[Rcpp::export()]]
double volz_loglik_nh2(List init, arma::mat f1, arma::vec betaN, double t_correct, arma::ivec index,
                       std::string transX = "standard"){

  // std::string transX = "standard"
  // the original volz_loglik_nh2
  /***
   * unused
   */
  Rcout << "wrong" << endl;
  int p = f1.n_cols - 1;
  int n2 = f1.n_rows; // 38
  int n0 = as<int>(init[9]) - 1; //37
  int L = 0;
  while(f1(L,0) < t_correct){
    L ++;
    if(L == n2 - 1) break;     //37
  }
  if(f1.submat(0,1,n0,p).min() < 0){
    return -10000000;
  }
  // Rcout<<f1 <<endl;
  arma::mat f2(n0,p);
  arma::vec betanh(n0);

  // original
  for(int i = 1; i<= n0; i++){
    f2.submat((n0-i),0,(n0-i),p - 1) = f1.submat(L-i+1,1,L-i+1,p);
    betanh(n0-i) = betaN(L-i+1);
  }

  // modified
  for(int i = 1; i<= n0; i++){
    f2.submat((n0-i),0,(n0-i),p - 1) = 0.5 * (f1.submat(L-i+1,1,L-i+1,p) + f1.submat(L-i,1,L-i ,p));
    betanh(n0-i) = betaN(L-i+1);
  }

  if(as<int>(init[9]) != (f2.n_rows+1)){
    Rcout<<"Incorrect length for f"<<endl;
  }

  arma::vec gridrep;
  gridrep = as<vec>(init[6]);
  int k = sum(as<arma::vec>(init[6]));

  arma::vec f(k);
  arma::vec s(k);
  arma::vec e(k);
  arma::vec betas(k);
  int start = 0;
  for(int i = 0; i < f2.n_rows; i++){
    for(int j = 0; j < gridrep(i);j++){
      if(transX == "standard"){
        f(start) = log(f2(f2.n_rows-i-1,index(1)));
        s(start) = log(f2(f2.n_rows-i-1,index(0)));
        e(start) = log(f2(f2.n_rows-i-1,1));
      }else if(transX == "log"){
        f(start) = f2(f2.n_rows-i-1,index(1));
        s(start) = f2(f2.n_rows-i-1,index(0));
        e(start) = f2(f2.n_rows-i-1,1);
      }
      betas(start) = betanh(f2.n_rows-i-1);
      start ++;
    }
  }

  arma::vec ll = -2 * (as<vec>(init[2]) % as<vec>(init[3]) % betas % (arma::exp(- f)) % (arma::exp(s))) +\
    as<vec>(init[4]) % (log(betas) -f + s);
  //Rcout<< ll <<endl;
  if(ll.has_nan()){
    return -10000000;
  }
  return sum(ll);
}


//[[Rcpp::export()]]
double volz_loglik_nh_adj(List init, arma::mat f1, arma::vec betaN, double t_correct, arma::ivec index,
                       bool enable){
  double res = 0;
  if(enable){
    int p = f1.n_cols - 1; // 3
    int n0 = as<int>(init[9]) - 1; // 37
    int nstart = as<int>(init[13]) - 1; //
    //while(f1(L,0) < t_correct){
    //  L ++;
    //  if(L == n2 - 1) break;
    //}
    if(f1.submat(nstart,1,n0,p).min() < 0){
      return -10000000;
    }
    arma::mat f2(n0-nstart+1,p);
    arma::vec betanh(n0-nstart+1);
    for(int i = 1; i<= n0-nstart+1; i++){
      f2.submat((n0-nstart-i+1),0,(n0-nstart-i+1),p - 1) = f1.submat(n0-i+1,1,n0-i+1,p);
      betanh(n0-nstart-i+1) = betaN(n0-i+1);
    }
    // Rcout<<f2.n_rows<<"\t"<<as<int>(init[9])<<endl;
    /*
    if(as<int>(init[9]) != (f2.n_rows+1)){
      Rcout<<"Incorrect length for f"<<endl;
    }
    */
    arma::vec gridrep;
    gridrep = as<vec>(init[6]);
    int k = sum(as<arma::vec>(init[6]));

    arma::vec f(k);
    arma::vec s(k);
    arma::vec e(k);
    arma::vec betas(k);
    int start = 0;
    for(int i = 0; i < f2.n_rows; i++){
      for(int j = 0; j < gridrep(i);j++){
        f(start) = log(f2(f2.n_rows-i-1,index(1)));
        s(start) = log(f2(f2.n_rows-i-1,index(0)));
        e(start) = log(f2(f2.n_rows-i-1,1));
        betas(start) = betanh(f2.n_rows-i-1);
        start ++;
      }
    }
    arma::vec ll = -2 * (as<vec>(init[2]) % as<vec>(init[3]) % betas % (arma::exp(- f)) % (arma::exp(s))) +\
      as<vec>(init[4]) % (log(betas) -f + s);
    //Rcout<< ll <<endl;
    if(ll.has_nan()){
      return -10000000;
    }
    res = sum(ll);
  }
  return res;
}


//[[Rcpp::export()]]
arma::mat Ode_Coarse_Slicer(arma::mat Ode_thin, int gridsize){
  int n = Ode_thin.n_rows;
  int p = Ode_thin.n_cols;
  int d = (n / gridsize) + 1;
  arma::mat Ode_Coarse(d,p);
  for(int i = 0; i < d; i ++){
    Ode_Coarse.row(i) = Ode_thin.row(i * gridsize);
  }
  return Ode_Coarse;
}
/*
arma::mat Ode_Coarse_Slicer(arma::mat Ode_thin, arma::ivec gridsizes){
  int n = Ode_thin.n_rows;
  int p = Ode_thin.n_cols;
  int d = gridsizes.n_elem + 1;
  arma::mat Ode_Coarse(d,p);
  int begin = 0;
  for(int i = 0; i < d; i ++){
    Ode_Coarse.row(i) = Ode_thin.row(begin);
    if(i < d - 1) {
      begin += gridsizes(i);
    }
  }
  return Ode_Coarse;
}

*/


arma::mat SIR_incidence_Traj(arma::mat Traj){
  // arma::vec traj = Traj.submat(0,1,0,2);
  int n = Traj.n_rows;
  arma::mat incidenceTraj(n - 1,3);
  incidenceTraj.col(0) = Traj.submat(1, 0 ,n - 1, 0);
  for(int i = 1; i < n; i ++){
    double NSI = Traj(i-1,1) - Traj(i,1);
    double NIR = NSI + Traj(i - 1,2) - Traj(i, 2);
    incidenceTraj(i - 1, 1) = NSI;
    incidenceTraj(i - 1, 2) = NIR;
  }
  return incidenceTraj;
}

//[[Rcpp::export()]]
arma::mat SIR_incidence_Traj(arma::mat Traj, arma::ivec gridrep){
  // arma::vec traj = Traj.submat(0,1,0,2);
  int n = gridrep.n_rows;
  arma::mat incidenceTraj(n - 1,3);
  for(int i = 1; i < n; i ++){
    double NSI = Traj(gridrep(i - 1),1) - Traj(gridrep(i),1);
    double NIR = NSI + Traj(gridrep(i - 1),2) - Traj(gridrep(i), 2);
    incidenceTraj(i - 1, 1) = NSI;
    incidenceTraj(i - 1, 2) = NIR;
    incidenceTraj(i - 1,0) = Traj(gridrep(i),0);
  }
  return incidenceTraj;
}

double incidence_core(arma::mat Incidence, arma::mat Incidence_obs, arma::vec Incid_Par, bool Pois = false){
  int n = Incidence_obs.n_rows;
  double logprob = 0;
  for(int i = 0; i < n ; i ++){
    if(Incidence(i,1) >= 0){
      double NSI = (Incidence(i,1) > 0) ? Incidence(i,1) : 0;
      double mu = Incid_Par[0] * NSI + 0.00000001;
      if(Pois){
        logprob += R::dpois(Incidence_obs(i,1), mu, 1);
      }else{
        logprob += R::dnbinom(Incidence_obs(i,1), Incid_Par(1), Incid_Par(1) / (Incid_Par(1) + mu), 1);
      }
    }else{
      logprob = -10000000;
      break;
    }
  }
  return logprob;
}


double pref_core(arma::mat Incidence, arma::vec Samp_count, arma::vec Pref_par, bool Pois = false){
  int n = Samp_count.n_elem;
  double logprob = 0;
  double beta0 = Pref_par[1];
  double beta1 = Pref_par[0];
  for(int i = 0; i < n; i ++){
    if(Incidence(i,1) > 0){
      double NSI = (Incidence(i,1) > 0) ? Incidence(i,1) : 0;
      double mu = exp(beta0 + beta1 * log(NSI+1));
      if(Pois){
        logprob += R::dpois(Samp_count(i), mu, 1);
      }else{
        logprob += R::dnbinom(Samp_count(i), Pref_par(2), Pref_par(2) / (Pref_par(2) + mu), 1);
      }
    }else{
      logprob = -10000000;
      break;
    }
  }
  return logprob;
}

//[[Rcpp::export()]]
double log_pref(arma::mat Traj, List & SampleInit, arma::vec Pref_par, bool Pois = false, bool enable = true){
  double res = 0;
  if(enable){
    arma::vec Samp_count = as<arma::vec>(SampleInit[0]);
    arma::ivec grid_idx = as<arma::ivec>(SampleInit[1]);
    arma::mat Incid_Traj = SIR_incidence_Traj(Traj, grid_idx);
    res = pref_core(Incid_Traj, Samp_count, Pref_par, Pois);
  }
  return res;
}
//[[Rcpp::export()]]
arma::mat SIR_Traj(arma::mat Traj, arma::ivec gridrep){
  // arma::vec traj = Traj.submat(0,1,0,2);
  int n = gridrep.n_rows - 1;
  arma::mat incidenceTraj(n,2);
  for(int i = 0; i < n; i ++){
    incidenceTraj(i,0) = Traj(gridrep(i), 0);
    incidenceTraj(i,1) = Traj(gridrep(i), 2);
  }
  return incidenceTraj;
}

//[[Rcpp::export()]]
double log_prev(arma::mat Traj, List & SampleInit, arma::vec Pref_par, bool Pois = false, bool enable = true){
  double res = 0;
  if(enable){
    arma::vec Samp_count = as<arma::vec>(SampleInit[0]);
    arma::ivec grid_idx = as<arma::ivec>(SampleInit[1]);
    arma::mat Prev_Traj = SIR_Traj(Traj, grid_idx);
    res = pref_core(Prev_Traj, Samp_count, Pref_par, Pois);
  }
  return res;
}

double log_incidence(arma::mat Traj, arma::mat Incidence_obs, arma::vec Incid_Par, bool Pois = false){
  /*
   * Traj: Trajectory first column is t, other columns are S I trajectories
   *  Incidence_obs: Incidence data that matches the time grid
   * rho: report rate
   *
   */
  // transform the SIR trajectory to incidence data
  arma::mat Incidence = SIR_incidence_Traj(Traj);
  int n = Incidence_obs.n_rows;
  double logprob = 0;
  for(int i = 0; i < n ; i ++){
    if(Incidence(i,1) >= 0){
      double NSI = (Incidence(i,1) > 0) ? Incidence(i,1) : 0;
      double mu = Incid_Par[0] * NSI + 0.00000001;
      if(Pois){
        logprob += R::dpois(Incidence_obs(i,1), mu, 1);
      }else{
        logprob += R::dnbinom(Incidence_obs(i,1), Incid_Par(1), Incid_Par(1) / (Incid_Par(1) + mu), 1);
      }
    }else{
      logprob = -10000000;
      break;
    }
  }
  return logprob;
}


//[[Rcpp::export()]]
double log_incidence_approx(arma::vec betat, arma::mat Traj, arma::mat Incidence_obs, arma::vec incid_par, arma::vec dts, bool Pois = false){
  /*
  * Traj: Trajectory first column is t, other columns are S I trajectories
  *  Incidence_obs: Incidence data that matches the time grid
  * rho: report rate
  * dts: times to generate the incidence data
  *
  */

  // transform the SIR trajectory to incidence data
  int i = 0; // index for traj
  int j = 0; // index for incidence
  int n = Incidence_obs.n_rows;
  int m = Traj.n_rows;
  double logProb = 0;
  double Dt = Traj(1,0) - Traj(0,0);
  while(Traj(i,0) > Incidence_obs(j,0)){
    j ++;
  }
  while(i < m){
    if(j == n){
      break;
    }
    if(Traj(i,0) < Incidence_obs(j,0)){
      if(Traj(i,0) + Dt >= Incidence_obs(j,0)){
        double rate = incid_par(0) * betat(i) * Traj(i,1) * Traj(i,2) * dts(j);
        if(Pois){
          logProb += R::dpois(Incidence_obs(j,1), rate, 1);
        }else{
          logProb += R::dnbinom(Incidence_obs(j,1), incid_par(1), incid_par(1) / (rate + incid_par(1)), 1);
        }
        j ++;
      }else{
        i ++;
      }
    }
  }
  return logProb;
}



//[[Rcpp::export()]]
List New_Param_List(arma::vec param, arma::vec initial, int gridsize, arma::vec t, arma::vec x_r, arma::ivec x_i,
                    std::string transP = "changepoint",
                    std::string model = "SIR", std::string transX = "standard"){

  arma::mat OdeTraj_thin = ODE_rk45(initial,t, param,
                                   x_r, x_i, transP, model, transX);

  List FT_new = KF_param_chol(OdeTraj_thin, param, gridsize, x_r, x_i, transP, model, transX);

  arma::mat Ode_Coarse = Ode_Coarse_Slicer(OdeTraj_thin, gridsize);
  arma::vec betaNs = betaTs(param,Ode_Coarse.col(0),x_r, x_i);
  List result;
  result["FT"] = FT_new;
  result["Ode"] = Ode_Coarse;
  result["betaN"] = betaNs;

  return result;
}


//[[Rcpp::export()]]
List Update_Param(arma::vec param, arma::vec initial, arma::vec t, arma::mat OriginTraj,
                  arma::vec x_r, arma::ivec x_i, List init, int gridsize,
                  double coal_log=0, double prior_proposal_offset = 0, double t_correct = 0, std::string transP = "changepoint",
                  std::string model = "SIR", std::string transX = "standard", bool volz = true){

  arma::mat OdeTraj_thin = ODE_rk45(initial,t, param,
                                    x_r, x_i, transP, model, transX);

  List FT_new = KF_param_chol(OdeTraj_thin, param, gridsize, x_r, x_i, transP, model, transX);

  arma::mat Ode_Coarse = Ode_Coarse_Slicer(OdeTraj_thin, gridsize);
  arma::vec betaNs = betaTs(param,Ode_Coarse.col(0),x_r, x_i);

  arma::mat NewTraj = TransformTraj(Ode_Coarse, OriginTraj, FT_new);
  double coal_log_new = volz_loglik_nh2(init, NewTraj,betaNs, t_correct, x_i.subvec(2,3) ,transX);

  double a = coal_log_new - coal_log + prior_proposal_offset;
  List Result;
  if(log(R::runif(0,1)) < a){

    Result["accept"] = true;
    Result["FT"] = FT_new;
    Result["Ode"] = Ode_Coarse;
    Result["betaN"] = betaNs;
    Result["coalLog"] = coal_log_new;
    Result["LatentTraj"] = NewTraj;

  }else{
    Result["accept"] = false;
  }
  return Result;
}


List Update_Param_JointData(arma::vec param, arma::vec initial, arma::vec incid_par, arma::vec t, arma::mat OriginTraj,
                  arma::vec x_r, arma::ivec x_i, List init, arma::mat Incidence_data, int gridsize,
                  double coal_log=0, double incid_log = 0, double prior_proposal_offset = 0, double t_correct = 0, std::string transP = "changepoint",
                  std::string model = "SIR", std::string transX = "standard", bool volz = true){

  arma::mat OdeTraj_thin = ODE_rk45(initial,t, param,
                                    x_r, x_i, transP, model, transX);

  List FT_new = KF_param_chol(OdeTraj_thin, param, gridsize, x_r, x_i, transP, model, transX);

  arma::mat Ode_Coarse = Ode_Coarse_Slicer(OdeTraj_thin, gridsize);

  arma::vec betaNs = betaTs(param,Ode_Coarse.col(0),x_r, x_i);

  arma::mat NewTraj = TransformTraj(Ode_Coarse, OriginTraj, FT_new);

  double coal_log_new = volz_loglik_nh2(init, NewTraj,betaNs, t_correct, x_i.subvec(2,3) ,transX);

  double incid_log_new = log_incidence(NewTraj, Incidence_data, incid_par);

  double a = coal_log_new - coal_log + incid_log_new - incid_log + prior_proposal_offset;

  List Result;

  if(log(R::runif(0,1)) < a){

    Result["accept"] = true;
    Result["FT"] = FT_new;
    Result["Ode"] = Ode_Coarse;
    Result["betaN"] = betaNs;
    Result["coalLog"] = coal_log_new;
    Result["LatentTraj"] = NewTraj;
    Result["incidLog"] = incid_log_new;
  }else{
    Result["accept"] = false;
  }
  return Result;
}










//[[Rcpp::export()]]
arma::vec Param_Slice_update(arma::vec param, arma::vec x_r, arma::ivec x_i, double theta, arma::vec newChs, double rho = 1){

  arma::vec param_new = param;
  arma::vec OdeChs = arma::log(param.subvec(x_i(1),x_i(1) + x_i(0) - 1));
  arma::vec Chs_prime = OdeChs * cos(theta) + newChs * sin(theta);
  param_new.subvec(x_i(1), x_i(0) + x_i(1) - 1) = arma::exp(Chs_prime);
  return param_new;
}


//[[Rcpp::export()]]
arma::vec Param_Slice_update_all(arma::vec par_old, arma::vec x_r, arma::ivec x_i, double theta, arma::vec newChs, double rho = 1){

  arma::vec par_new = par_old;
  arma::vec OdeChs = arma::log(par_old.subvec(1,x_i(1) + x_i(0) + 2));
  arma::vec Chs_prime = OdeChs * cos(theta) + newChs * sin(theta);
  par_new.subvec(1, x_i(0) + x_i(1) + 1) = arma::exp(Chs_prime);
  return par_new;
}


//[[Rcpp::export()]]
arma::vec Param_Slice_update_all2(arma::vec par_old, arma::vec x_r, arma::ivec x_i, double theta, arma::vec newChs, arma::ivec ESS_vec, List priorList){
  /*
   * new_chs: vector of normal ((log(I0_new) - pr11)/pr12, (log(R0_new) - pr21_/ pr22, (log(gamma_new) - pr31)/ pr32,
   * log(changepoint_new) / hyper_new, log(hyper_new))
   *
   * return: vector par_new that can be attached to MCMC_obj
   */

  arma::vec par_new = par_old;


  // extract parameters for prior distribution
  arma::vec pr(8);
  pr.subvec(0,1) = as<arma::vec>(priorList[0]);
  pr.subvec(2,3) = as<arma::vec>(priorList[1]);
  pr.subvec(4,5) = as<arma::vec>(priorList[2]);
  pr.subvec(6,7) = as<arma::vec>(priorList[3]);


  // extract Normal(0,1) variable from par_old
  arma::vec OdeChs = arma::log(par_old.subvec(1,x_i(1) + x_i(0) + 2));
  // I_0 R_0 gamma
  for(int i = 0; i < 3; i ++){
      OdeChs(i) = (OdeChs(i) - pr(2 * i)) / pr(2 * i + 1);
  }
   // hyper
  OdeChs(x_i(0) + x_i(1) + 1) =  (OdeChs(x_i(0) + x_i(1) + 1) - pr(6)) / pr(7);
  // changepoints
  OdeChs.subvec(3,2 + x_i(0)) = OdeChs.subvec(3,2 + x_i(0)) * par_new(x_i(0) + x_i(1) + 2);

  // rotation in Elliptical slice sampler
  arma::vec Chs_prime = OdeChs * cos(theta) + newChs * sin(theta);

  for(int i = 0; i < 3; i ++){
    if(ESS_vec(i) != 0){
      par_new(1 + i) = exp(Chs_prime(i) * pr(2 * i + 1) + pr(2 * i));
    }
  }
  if(ESS_vec(4) != 0){
    par_new(x_i(0) + x_i(1) + 2) = exp(Chs_prime(x_i(0) + x_i(1) + 1) * pr(7) + pr(6));
  }
  if(ESS_vec(5) != 0){
    par_new.subvec(4, x_i(0) + x_i(1) + 1) = arma::exp(Chs_prime.subvec(3,x_i(0) + x_i(1)) / par_new(x_i(0) + x_i(1) + 2));
  }

  return par_new;
}

//[[Rcpp::export()]]
List ESlice_change_points(arma::vec param, arma::vec initial, arma::vec t, arma::mat OriginTraj,
                          arma::vec x_r, arma::ivec x_i, List init, int gridsize,
                          double coal_log=0, double t_correct = 0, std::string transP = "changepoint",
                          std::string model = "SIR", std::string transX = "standard", bool volz = true){

  int nch = x_i[0];
  arma::ivec Index(2);
  if(model == "SIR"){
    Index(0) = 0; Index(1) = 1;
  }else if(model == "SEIR"){
    Index(0) = 0; Index(1) = 2;
  }

  //param.subvec(x_i(1),x_i(1) + x_i(0) - 1)
  double u = R::runif(0,1);

  if(log(u) < -25){
    Rcout << "really small u" << endl;
    Rcout << u << endl;
  }

  double logy = coal_log + log(u);

  double theta = R::runif(0,2*pi);
  double theta_min = theta - 2*pi;
  double theta_max = theta;

  arma::vec newChs = arma::randn(nch,1) / param(x_i(0) + x_i(1));
  arma::vec param_new = Param_Slice_update(param, x_r, x_i, theta, newChs);

  List param_list = New_Param_List(param_new, initial, gridsize, t, x_r, x_i,
                                   transP, model, transX);

  List FT_new = as<Rcpp::List>(param_list[0]);

  arma::vec betaN = as<arma::vec>(param_list[2]);

  arma::mat OdeTraj_new  = as<arma::mat>(param_list[1]);
  arma::mat NewTraj = TransformTraj(OdeTraj_new, OriginTraj, FT_new);

  double loglike = volz_loglik_nh2(init, NewTraj,betaN,t_correct,Index ,transX);
  int i = 0;

  while(loglike <= logy){
    i += 1;

    if(i>20){
      theta = 0;
      param_new = param;
      param_list = New_Param_List(param, initial, gridsize, t, x_r, x_i,
                                  transP, model, transX);
      FT_new = as<Rcpp::List>(param_list[0]);
      betaN = as<arma::vec>(param_list[2]);
      OdeTraj_new  = as<arma::mat>(param_list[1]);
      NewTraj = TransformTraj(OdeTraj_new, OriginTraj, FT_new);
      loglike = volz_loglik_nh2(init, NewTraj,betaN,t_correct,Index ,transX);
      break;
    }
    if(theta < 0){
      theta_min = theta;
    }else{
      theta_max = theta;
    }

    theta = R::runif(theta_min,theta_max);

    param_new = Param_Slice_update(param, x_r, x_i, theta, newChs);

    param_list = New_Param_List(param_new, initial, gridsize, t, x_r, x_i,
                                     transP, model, transX);

    FT_new = as<Rcpp::List>(param_list[0]);
    betaN = as<arma::vec>(param_list[2]);
    OdeTraj_new  = as<arma::mat>(param_list[1]);
    NewTraj = TransformTraj(OdeTraj_new, OriginTraj, FT_new);

    loglike = volz_loglik_nh2(init, NewTraj, betaN, t_correct, Index ,transX);
  }
  arma::vec delta = arma::log(param_new.subvec(x_i(1),x_i(0) + x_i(1) - 1)) * param_new(x_i(0) + x_i(1));
  double LogChProb = 0;
  for(int j = 0; j < delta.n_elem; j ++){
    LogChProb += -0.5 * delta(j) * delta(j);
  }
  List result;
  result["betaN"] = betaN;
  result["FT"] = FT_new;
  result["OdeTraj"] = OdeTraj_new;
  result["param"] = param_new;
  result["LatentTraj"] = NewTraj;
  result["CoalLog"] = loglike;
  result["LogChProb"] = LogChProb;
  return result;
}


//[[Rcpp::export()]]
List ESlice_par(arma::vec param, arma::vec initial, arma::vec t, arma::mat OriginTraj,List priorList,
                          arma::vec x_r, arma::ivec x_i, List init, int gridsize,
                          double coal_log=0, double t_correct = 0, std::string transP = "changepoint",
                          std::string model = "SIR", std::string transX = "standard", bool volz = true){

  int nch = x_i[0];

  arma::ivec Index(2);

  arma::vec I_pr = as<arma::vec>(priorList[0]);
  arma::vec R0_pr = as<arma::vec>(priorList[1]);
  arma::vec gamma_pr = as<arma::vec>(priorList[2]);

  if(model == "SIR"){
    Index(0) = 0; Index(1) = 1;
  }else if(model == "SEIR"){
    Index(0) = 0; Index(1) = 2;
  }

  //param.subvec(x_i(1),x_i(1) + x_i(0) - 1)
  double u = R::runif(0,1);

  if(log(u) < -25){
    Rcout << "really small u" << endl;
    Rcout << u << endl;
  }

  double logy = coal_log + log(u);

    //+ R::dnorm(initial(1), I_pr(0), I_pr(1),1) + \
    //R::dnorm(param(0), R0_pr(0), R0_pr(1),1) + \
    //R::dnorm(param(1), gamma_pr(0), gamma_pr(1),1);

  double theta = R::runif(0,2*pi);
  double theta_min = theta - 2*pi;
  double theta_max = theta;
  arma::vec new_par_raw = arma::randn(3 + nch, 1); // 47 * 1

  new_par_raw[0] = (new_par_raw[0]) * I_pr[1] + I_pr[0];
  new_par_raw[1] = (new_par_raw[1]) * R0_pr[1] + R0_pr[0];
  new_par_raw[2] = (new_par_raw[2]) * gamma_pr[1] + gamma_pr[0];

  new_par_raw.subvec(3, 2 + nch) = new_par_raw.subvec(3, 2 + nch) / param(x_i(0) + x_i(1));
  arma::vec par_old(initial.n_elem + x_i(0) + x_i(1) + 1); // 49

  par_old.subvec(0,1) = initial;
  par_old.subvec(2, x_i(0) + x_i(1) + 2) = param;
  arma::vec par_new = Param_Slice_update_all(par_old, x_r, x_i, theta, new_par_raw); // 49

  List param_list = New_Param_List(par_new.subvec(2, x_i(0) + x_i(1) + 2), par_new.subvec(0,1), gridsize, t, x_r, x_i,
                                   transP, model, transX);
  List FT_new = as<Rcpp::List>(param_list[0]);

  arma::vec betaN = as<arma::vec>(param_list[2]);

  arma::mat OdeTraj_new  = as<arma::mat>(param_list[1]);
  arma::mat NewTraj = TransformTraj(OdeTraj_new, OriginTraj, FT_new);

  double loglike = volz_loglik_nh2(init, NewTraj,betaN,t_correct,Index ,transX);
  int i = 0;

  while(loglike <= logy){
    i += 1;

    if(i>20){
      theta = 0;
      par_new.subvec(2, x_i(0) + x_i(1) + 2) = param;
      par_new.subvec(0,1) = initial;
      param_list = New_Param_List(param, initial, gridsize, t, x_r, x_i,
                                  transP, model, transX);
      FT_new = as<Rcpp::List>(param_list[0]);
      betaN = as<arma::vec>(param_list[2]);
      OdeTraj_new  = as<arma::mat>(param_list[1]);
      NewTraj = TransformTraj(OdeTraj_new, OriginTraj, FT_new);
      loglike = volz_loglik_nh2(init, NewTraj,betaN,t_correct,Index ,transX);
      break;
    }
    if(theta < 0){
      theta_min = theta;
    }else{
      theta_max = theta;
    }

    theta = R::runif(theta_min,theta_max);

    par_new = Param_Slice_update_all(par_old, x_r, x_i, theta, new_par_raw);

    param_list = New_Param_List(par_new.subvec(2, x_i(0) + x_i(1) + 2), par_new.subvec(0,1), gridsize, t, x_r, x_i,
                                transP, model, transX);

    FT_new = as<Rcpp::List>(param_list[0]);
    betaN = as<arma::vec>(param_list[2]);
    OdeTraj_new  = as<arma::mat>(param_list[1]);
    NewTraj = TransformTraj(OdeTraj_new, OriginTraj, FT_new);

    loglike = volz_loglik_nh2(init, NewTraj, betaN, t_correct, Index ,transX);
  }
  List result;
  arma::vec param_new = par_new.subvec(2, x_i(0) + x_i(1) + 2);
  arma::vec initial_new = par_new.subvec(0,1);
  arma::vec delta = arma::log(param_new.subvec(x_i(1),x_i(0) + x_i(1) - 1)) * param_new(x_i(0) + x_i(1));
  double LogChProb = 0;
  for(int j = 0; j < delta.n_elem; j ++){
    LogChProb += -0.5 * delta(j) * delta(j);
  }
  result["betaN"] = betaN;
  result["FT"] = FT_new;
  result["OdeTraj"] = OdeTraj_new;
  result["param"] = param_new;
  result["LatentTraj"] = NewTraj;
  result["CoalLog"] = loglike;
  result["initial"] = initial_new;
  result["LogChProb"] = LogChProb;
  return result;
}



//[[Rcpp::export()]]
List ESlice_par_General(arma::vec par_old, arma::vec t, arma::mat OriginTraj,List priorList,
                arma::vec x_r, arma::ivec x_i, List init, int gridsize, arma::ivec ESS_vec,
                double coal_log = 0, double t_correct = 0, std::string transP = "changepoint",
                std::string model = "SIR", std::string transX = "standard", bool volz = true){

  /*
   *
   *
   *  priorList: a list of priors for parameters for 1. I_0, 2. R_0, 3. gamma, 4, mu (no exist in SIR model) 5 hyper
   * all of them follows lognormal distribution
   *
   * ESS vec: vector that stores the index for each
   *
   *  test R code:
   *
   *  par_old = resres5MJP$par[125000,]
      A = ESlice_par_General(par_old, seq(0,100,length.out = 2001), resres5MJP$MCMC_obj$OriginTraj, priorList = list(a=c(1,1),b=c(0.7,0.2), c=c(2,1), d=c(3,0.2)),
                          resres5MJP$MCMC_setting$x_r, resres5MJP$MCMC_setting$x_i, resres5MJP$MCMC_setting$Init,50,ESS_vec = c(0,1,0,0,1,1,0),
        resres5MJP$MCMC_obj$coalLog, 90)

      A$par[c(1,2,4)] - par_old[c(1,2,4)]
   *
   */

  arma::ivec Index(2);

  int u_param = par_old.n_elem - 1; // length(par) - 1

  // extract prior distribution for each parameter

//  arma::vec R0_pr = as<arma::vec>(priorList[1]);
//  arma::vec gamma_pr = as<arma::vec>(priorList[2]);
//  arma::vec hyper_pr = as<arma::vec>(priorList[3]);


  if(model == "SIR"){
    Index(0) = 0; Index(1) = 1;
  }else if(model == "SEIR"){
    Index(0) = 0; Index(1) = 2;
  }

  arma::vec new_par_raw = arma::randn(u_param, 1); // 47 * 1

  int N = t.n_elem / gridsize + 1;
  if(ESS_vec(6) == 1){
    arma::mat OriginTraj_new_raw = arma::randn(OriginTraj.n_rows, OriginTraj.n_cols);
  }

  arma::mat OriginTraj_new = OriginTraj;


  double u = R::runif(0,1);

  if(log(u) < -25){
    Rcout << "really small u" << endl;
    Rcout << u << endl;
  }

  double logy = coal_log + log(u);
  double theta = 2 * pi;
  double theta_min = 0;
  double theta_max = 2 * pi;

  arma::vec par_new = par_old;
  List param_list;
  List FT_new;
  arma::vec betaN(N);
  arma::mat OdeTraj_new(N,3);
  arma::mat NewTraj(N,3);
  double loglike;
  int i = 0;

  do{

    i += 1;

    // use the original values if no upates in 20 iterations
    if(i>20){
      Rcout<<"theta = "<< theta <<endl;
      theta = 0;
      par_new = par_old;
      param_list = New_Param_List(par_new.subvec(2, x_i(0) + x_i(1) + 2), par_new.subvec(0,1), gridsize, t, x_r, x_i,
                                   transP, model, transX);
      FT_new = as<Rcpp::List>(param_list[0]);
      betaN = as<arma::vec>(param_list[2]);
      OdeTraj_new  = as<arma::mat>(param_list[1]);
      OriginTraj_new = OriginTraj;
      NewTraj = TransformTraj(OdeTraj_new, OriginTraj_new, FT_new);
      loglike = volz_loglik_nh2(init, NewTraj,betaN,t_correct,Index ,transX);
      break;
    }

    if(i == 1){ // first iteration
      theta = R::runif(0,2*pi);
      theta_min = theta - 2*pi;
      theta_max = theta;
    }else{

      if(theta < 0){
        theta_min = theta;
      }else{
        theta_max = theta;
      }
      theta = R::runif(theta_min,theta_max);
    }

    par_new = Param_Slice_update_all2(par_old, x_r, x_i, theta, new_par_raw, ESS_vec, priorList); // 49
    int error = 0;

    try{
    param_list = New_Param_List(par_new.subvec(2, x_i(0) + x_i(1) + 2), par_new.subvec(0,1), gridsize, t, x_r, x_i,
                                   transP, model, transX);
    }catch(...){
      error = 1;
    }

    if(error == 1) { // if error happens in parameters, go back
      loglike = -10000000;
      continue;
    }

    FT_new = as<Rcpp::List>(param_list[0]);

     betaN = as<arma::vec>(param_list[2]);

     OdeTraj_new  = as<arma::mat>(param_list[1]);

  // Update OriginTraj if ESS_vec(6) == 1
    if(ESS_vec(6) == 1){
      OriginTraj_new = OriginTraj * cos(theta) + OriginTraj_new * sin(theta);
    }

    NewTraj = TransformTraj(OdeTraj_new, OriginTraj_new, FT_new);

    loglike = volz_loglik_nh2(init, NewTraj, betaN, t_correct, Index ,transX);

  }
  while(loglike <= logy);


  List result;
  double logOrigin = 0;
  for(int j = 0; j < OriginTraj_new.n_rows; j ++){
    for(int k = 0; k < OriginTraj_new.n_cols; k ++){
      logOrigin -= 0.5 * OriginTraj_new(j,k) * OriginTraj_new(j,k);
    }
  }
  result["betaN"] = betaN;
  result["FT"] = FT_new;
  result["OdeTraj"] = OdeTraj_new;
  result["par"] = par_new;
  result["LatentTraj"] = NewTraj;
  result["CoalLog"] = loglike;
  result["OriginTraj"] = OriginTraj_new;
  result["logOrigin"] = logOrigin;
  return result;
}



List ESlice_par_General_JointData(arma::vec par_old, arma::vec incid_par, arma::vec t, arma::mat OriginTraj,List priorList,
                        arma::vec x_r, arma::ivec x_i, List init, arma::mat Incidence_obs,int gridsize, arma::ivec ESS_vec,
                        double coal_log = 0, double incid_log = 0, double t_correct = 0, std::string transP = "changepoint",
                        std::string model = "SIR", std::string transX = "standard", bool volz = true){

  /*
  *
  *
  *  priorList: a list of priors for parameters for 1. I_0, 2. R_0, 3. gamma, 4, mu (no exist in SIR model) 5 hyper
  * all of them follows lognormal distribution
  *
  * ESS vec: vector that stores the index for each
  *
  *  test R code:
  *
  *  par_old = resres5MJP$par[125000,]
  A = ESlice_par_General(par_old, seq(0,100,length.out = 2001), resres5MJP$MCMC_obj$OriginTraj, priorList = list(a=c(1,1),b=c(0.7,0.2), c=c(2,1), d=c(3,0.2)),
                         resres5MJP$MCMC_setting$x_r, resres5MJP$MCMC_setting$x_i, resres5MJP$MCMC_setting$Init,50,ESS_vec = c(0,1,0,0,1,1,0),
  resres5MJP$MCMC_obj$coalLog, 90)

  A$par[c(1,2,4)] - par_old[c(1,2,4)]
  *
  */

  arma::ivec Index(2);

  int u_param = par_old.n_elem - 1; // length(par) - 1

  // extract prior distribution for each parameter

  //  arma::vec R0_pr = as<arma::vec>(priorList[1]);
  //  arma::vec gamma_pr = as<arma::vec>(priorList[2]);
  //  arma::vec hyper_pr = as<arma::vec>(priorList[3]);


  if(model == "SIR"){
    Index(0) = 0; Index(1) = 1;
  }else if(model == "SEIR"){
    Index(0) = 0; Index(1) = 2;
  }

  arma::vec new_par_raw = arma::randn(u_param, 1); // 47 * 1

  int N = t.n_elem / gridsize + 1;
  if(ESS_vec(6) == 1){
    arma::mat OriginTraj_new_raw = arma::randn(OriginTraj.n_rows, OriginTraj.n_cols);
  }

  arma::mat OriginTraj_new = OriginTraj;


  double u = R::runif(0,1);

  if(log(u) < -25){
    Rcout << "really small u" << endl;
    Rcout << u << endl;
  }

  double logy = coal_log + incid_log + log(u);
  double theta = 2 * pi;
  double theta_min = 0;
  double theta_max = 2 * pi;

  arma::vec par_new = par_old;
  List param_list;
  List FT_new;
  arma::vec betaN(N);
  arma::mat OdeTraj_new(N,3);
  arma::mat NewTraj(N,3);
  double coal_log_new;
  double incid_log_new;
  int i = 0;

  do{

    i += 1;

    // use the original values if no upates in 20 iterations
    if(i>20){
      //Rcout<<"theta = "<< theta <<endl;
      theta = 0;
      par_new = par_old;
      param_list = New_Param_List(par_new.subvec(2, x_i(0) + x_i(1) + 2), par_new.subvec(0,1), gridsize, t, x_r, x_i,
                                  transP, model, transX);
      FT_new = as<Rcpp::List>(param_list[0]);
      betaN = as<arma::vec>(param_list[2]);
      OdeTraj_new  = as<arma::mat>(param_list[1]);
      OriginTraj_new = OriginTraj;
      NewTraj = TransformTraj(OdeTraj_new, OriginTraj_new, FT_new);
      coal_log_new = volz_loglik_nh2(init, NewTraj,betaN,t_correct,Index ,transX);
      incid_log_new = log_incidence(NewTraj, Incidence_obs, incid_par);
      //Rcout << "coal diff" <<coal_log_new - coal_log << "incid diff" << incid_log_new - incid_log << endl;
      break;
    }

    if(i == 1){ // first iteration
      theta = R::runif(0,2*pi);
      theta_min = theta - 2*pi;
      theta_max = theta;
    }else{

      if(theta < 0){
        theta_min = theta;
      }else{
        theta_max = theta;
      }
      theta = R::runif(theta_min,theta_max);
    }

    par_new = Param_Slice_update_all2(par_old, x_r, x_i, theta, new_par_raw, ESS_vec, priorList); // 49
    int error = 0;

    try{
      param_list = New_Param_List(par_new.subvec(2, x_i(0) + x_i(1) + 2), par_new.subvec(0,1), gridsize, t, x_r, x_i,
                                  transP, model, transX);
    }catch(...){
      error = 1;
    }

    if(error == 1) { // if error happens in parameters, go back
      coal_log_new = -10000000;
      incid_log_new = -10000000;
      continue;
    }

    FT_new = as<Rcpp::List>(param_list[0]);

    betaN = as<arma::vec>(param_list[2]);

    OdeTraj_new  = as<arma::mat>(param_list[1]);

    // Update OriginTraj if ESS_vec(6) == 1
    if(ESS_vec(6) == 1){
      OriginTraj_new = OriginTraj * cos(theta) + OriginTraj_new * sin(theta);
    }

    NewTraj = TransformTraj(OdeTraj_new, OriginTraj_new, FT_new);

    coal_log_new = volz_loglik_nh2(init, NewTraj, betaN, t_correct, Index ,transX);
    incid_log_new = log_incidence(NewTraj, Incidence_obs, incid_par);
  }
  while(coal_log_new + incid_log_new <= logy);


  List result;
  double logOrigin = 0;
  for(int j = 0; j < OriginTraj_new.n_rows; j ++){
    for(int k = 0; k < OriginTraj_new.n_cols; k ++){
      logOrigin -= 0.5 * OriginTraj_new(j,k) * OriginTraj_new(j,k);
    }
  }
  result["betaN"] = betaN;
  result["FT"] = FT_new;
  result["OdeTraj"] = OdeTraj_new;
  result["par"] = par_new;
  result["LatentTraj"] = NewTraj;
  result["CoalLog"] = coal_log_new;
  result["IncidLog"] = incid_log_new;
  result["OriginTraj"] = OriginTraj_new;
  result["logOrigin"] = logOrigin;
  return result;
}
//[[Rcpp::export()]]
List ESlice_general_NC(arma::mat f_cur, arma::mat OdeTraj, List FTs, arma::vec state,
                       List init, arma::vec betaN, double t_correct, double lambda=10,
                       double coal_log=-99999999, int gridsize = 100, bool volz = false, std::string model = "SIR",
                       std::string transX = "standard"){
  // OdeTraj is the one with low resolution
  arma::ivec Index(2);
  if(model == "SIR"){
    Index(0) = 0; Index(1) = 1;
  }else if(model == "SEIR"){
    Index(0) = 0; Index(1) = 2;
  }

  int p = f_cur.n_cols;
  arma::mat newTraj(f_cur.n_rows, f_cur.n_cols);
  double logy;

  arma::mat v_traj = arma::randn(f_cur.n_rows, f_cur.n_cols);
  double u = R::runif(0,1);
  if(log(u) < -25){
    Rcout << "really small u" << endl;
    Rcout << u << endl;
  }
  //  if(funname == "standard"){
  // logy = coal_loglik(init,LogTraj(f_cur),t_correct,lambda,gridsize) + log(u);
  if(coal_log != -99999999){
    logy = coal_log + log(u);
  }else{
    if(volz){
      logy = volz_loglik_nh2(init, TransformTraj(OdeTraj, f_cur, FTs), betaN, t_correct,Index,transX) + log(u);
    }else{
      logy = coal_loglik3(init,f_cur,t_correct,lambda,Index(1), transX) + log(u);
    }
  }

  double theta = R::runif(0,2 * pi);

  double theta_min = theta - 2*pi;
  double theta_max = theta;

  arma::mat f_prime = f_cur * cos(theta) + v_traj * sin(theta);
  newTraj = TransformTraj(OdeTraj,f_prime, FTs);
  int i = 0;
  double loglike;
  if(volz){
    loglike = volz_loglik_nh2(init, newTraj,betaN,t_correct,Index ,transX);
  }else{
    loglike = coal_loglik3(init,LogTraj(newTraj),t_correct,lambda,Index(1), transX);
  }
  while(newTraj.cols(1,p).min() <0 || loglike <= logy){
    // shrink the bracket
    i += 1;
    if(i>20){
      newTraj = TransformTraj(OdeTraj,f_cur, FTs);
      loglike = volz_loglik_nh2(init, newTraj,betaN,t_correct,Index ,transX);
      f_prime = f_cur;
      Rcout << "theta = "<< theta << endl;
      break;
    }
    if(theta < 0){
      theta_min = theta;
    }else{
      theta_max = theta;
    }
    theta = R::runif(theta_min,theta_max);
    f_prime = f_cur * cos(theta) + v_traj * sin(theta);
    newTraj = TransformTraj(OdeTraj,f_prime, FTs);

    if(volz){
      loglike = volz_loglik_nh2(init, newTraj, betaN, t_correct, Index ,transX);
    }else{
      loglike = coal_loglik3(init,LogTraj(newTraj),t_correct,lambda,Index(1), transX);
    }
    if(isnan(loglike)){
      Rcout << "NaN in Slice Sampler " << endl;
    }
  }
  List Res;
  double logOrigin = 0;
  for(int j = 0; j < f_cur.n_rows - 1; j ++){
    for(int k = 0; k < p; k ++){
      logOrigin -= 0.5 * f_prime(j,k) * f_prime(j,k);
    }
  }
  Res["LatentTraj"] = newTraj;
  Res["OriginTraj"] = f_prime;
  Res["logOrigin"] = logOrigin;
  Res["CoalLog"] = loglike;
  return Res;
}




List ESlice_general_NC_joint(arma::mat f_cur, arma::mat OdeTraj, List FTs, arma::vec state, arma::vec incid_par,
                         List init, arma::mat Incidence, arma::vec betaN, double t_correct, double lambda=10,
                         double coal_log=-99999999, double IncidLog = -9999999, int gridsize = 100, bool volz = false, std::string model = "SIR",
                         std::string transX = "standard"){
  // OdeTraj is the one with low resolution
  arma::ivec Index(2);
  if(model == "SIR"){
    Index(0) = 0; Index(1) = 1;
  }else if(model == "SEIR"){
    Index(0) = 0; Index(1) = 2;
  }

  int p = f_cur.n_cols;
  arma::mat newTraj(f_cur.n_rows, f_cur.n_cols);
  double logy;

    arma::mat v_traj = arma::randn(f_cur.n_rows, f_cur.n_cols);
    double u = R::runif(0,1);
    if(log(u) < -25){
      Rcout << "really small u" << endl;
      Rcout << u << endl;
    }
    //  if(funname == "standard"){
    // logy = coal_loglik(init,LogTraj(f_cur),t_correct,lambda,gridsize) + log(u);
    logy = coal_log + IncidLog + log(u);

    double theta = R::runif(0,2 * pi);

    double theta_min = theta - 2*pi;
    double theta_max = theta;

    arma::mat f_prime = f_cur * cos(theta) + v_traj * sin(theta);
    newTraj = TransformTraj(OdeTraj,f_prime, FTs);
    int i = 0;
    double loglike;
    double IncidLog_new;
    if(volz){
      loglike = volz_loglik_nh2(init, newTraj,betaN,t_correct,Index ,transX);
    }else{
      loglike = coal_loglik3(init,LogTraj(newTraj),t_correct,lambda,Index(1), transX);
    }
    IncidLog_new = log_incidence(newTraj, Incidence, incid_par);

    while(newTraj.cols(1,p).min() <0 || (loglike + IncidLog_new) <= logy){
      // shrink the bracket
      i += 1;
      if(i>20){
        newTraj = TransformTraj(OdeTraj,f_cur, FTs);
        loglike = volz_loglik_nh2(init, newTraj,betaN,t_correct,Index ,transX);
        IncidLog_new = log_incidence(newTraj, Incidence, incid_par);
        f_prime = f_cur;
        //Rcout << "theta = "<< theta << endl;
        //Rcout << "coal diff0" <<loglike - coal_log << "incid diff0" << IncidLog_new - IncidLog << endl;
        break;
      }
      if(theta < 0){
        theta_min = theta;
      }else{
        theta_max = theta;
      }
      theta = R::runif(theta_min,theta_max);
      f_prime = f_cur * cos(theta) + v_traj * sin(theta);
      newTraj = TransformTraj(OdeTraj,f_prime, FTs);

      if(volz){
        loglike = volz_loglik_nh2(init, newTraj, betaN, t_correct, Index ,transX);
      }else{
        loglike = coal_loglik3(init,LogTraj(newTraj),t_correct,lambda,Index(1), transX);
      }
      IncidLog_new = log_incidence(newTraj, Incidence, incid_par);
      if(isnan(loglike)){
        Rcout << "NaN in Slice Sampler " << endl;
      }
      if(isnan(IncidLog_new)){
        Rcout << "NaN in Slice Sampler " << endl;
      }
  }
  List Res;
  double logOrigin = 0;
  for(int j = 0; j < f_cur.n_rows - 1; j ++){
    for(int k = 0; k < p; k ++){
      logOrigin -= 0.5 * f_prime(j,k) * f_prime(j,k);
    }
  }
  Res["LatentTraj"] = newTraj;
  Res["OriginTraj"] = f_prime;
  Res["logOrigin"] = logOrigin;
  Res["CoalLog"] = loglike;
  Res["IncidLog"] = IncidLog_new;
  return Res;
}



//[[Rcpp::export()]]
arma::mat ESlice_general2(arma::mat f_cur, arma::mat OdeTraj, List FTs, arma::vec state,
                         List init, arma::vec betaN, double t_correct, double lambda=10,
                         int reps=1, int gridsize = 100, bool volz = false, std::string model = "SIR",
                         std::string transX = "standard"){
  // OdeTraj is the one with low resolution
  arma::ivec Index(2);
  if(model == "SIR"){
    Index(0) = 0; Index(1) = 1;
  }else if(model == "SEIR"){
    Index(0) = 0; Index(1) = 2;
  }

  int p = f_cur.n_cols - 1;
  arma::mat newTraj(f_cur.n_rows, f_cur.n_cols);
  double logy;

    // centered the old trajectory without time grid
    arma::mat f_cur_centered = f_cur.cols(1,p) - OdeTraj.cols(1,p);
    List v = Traj_sim_general3(OdeTraj,FTs,t_correct);
    arma::mat v_traj = as<mat>(v[0]).cols(1,p) -  OdeTraj.cols(1,p);
    if(v_traj.has_nan()){
      Rcout<<"NA in proposed trajectory"<<endl;
    }
    double u = R::runif(0,1);
    //  if(funname == "standard"){
    // logy = coal_loglik(init,LogTraj(f_cur),t_correct,lambda,gridsize) + log(u);
    if(volz){
      logy = volz_loglik_nh2(init, f_cur, betaN, t_correct,Index,transX) + log(u);
    }else{
      logy = coal_loglik3(init,f_cur,t_correct,lambda,Index(1), transX) + log(u);
    }
    //   }else{
    //     logy = coal_loglik(init,f_cur,t_correct,lambda,gridsize) + log(u);
    //    }
    double theta = R::runif(0,2 * pi);

    double theta_min = theta - 2*pi;
    double theta_max = theta;

    arma::mat f_prime = f_cur_centered * cos(theta) + v_traj * sin(theta);
    newTraj.col(0) = f_cur.col(0);
    newTraj.cols(1,p) = f_prime + OdeTraj.cols(1,p);
    int i = 0;
    double loglike;
    if(volz){
      loglike = volz_loglik_nh2(init, newTraj,betaN,t_correct,Index ,transX);
    }else{
      loglike = coal_loglik3(init,LogTraj(newTraj),t_correct,lambda,Index(1), transX);
    }
    while(newTraj.cols(1,p).min() <0 || loglike <= logy){
      // shrink the bracket
      i += 1;
      if(i>20){
        newTraj = f_cur;
        Rcout<<"theta = "<< theta <<endl;
        break;
      }
      if(theta < 0){
        theta_min = theta;
      }else{
        theta_max = theta;
      }
      theta = R::runif(theta_min,theta_max);
      f_prime = f_cur_centered * cos(theta) + v_traj * sin(theta);
      // newTraj.col(0) = f_cur.col(0);
      newTraj.cols(1,p) = f_prime + OdeTraj.cols(1,p);
      if(volz){
        loglike = volz_loglik_nh2(init, newTraj, betaN, t_correct, Index ,transX);
      }else{
        loglike = coal_loglik3(init,LogTraj(newTraj),t_correct,lambda,Index(1), transX);
      }
    }
    f_cur = newTraj;
  return newTraj;
}


/*
 *
 * overloading functions
 *
 */

//[[Rcpp::export()]]
double log_incidence(arma::mat Traj, List IncidenceList, arma::vec Incid_Par, bool Pois, bool enable){
  /*
   * Traj: Trajectory first column is t, other columns are S I trajectories
   *  Incidence_obs: Incidence data that matches the time grid
   * rho: report rate
   *
   */
  double res = 0;
  if(enable){
    arma::mat Incidence_obs = as<arma::mat>(IncidenceList[0]);
    arma::ivec gridrep = as<arma::ivec>(IncidenceList[1]);
    res = incidence_core(SIR_incidence_Traj(Traj, gridrep), Incidence_obs, Incid_Par, Pois);
  }
  return res;
}


//[[Rcpp::export()]]
List Update_Param_JointData(arma::vec param, arma::vec initial, arma::vec incid_par, arma::vec t, arma::mat OriginTraj,
                            arma::vec x_r, arma::ivec x_i, List init, List Incidence_data, int gridsize,
                            double coal_log=0, double incid_log = 0, double prior_proposal_offset = 0, double t_correct = 0, std::string transP = "changepoint",
                            std::string model = "SIR", std::string transX = "standard", bool volz = true, bool addCoal = true, bool addIncid = true){

  arma::mat OdeTraj_thin = ODE_rk45(initial,t, param,
                                    x_r, x_i, transP, model, transX);

  List FT_new = KF_param_chol(OdeTraj_thin, param, gridsize, x_r, x_i, transP, model, transX);

  arma::mat Ode_Coarse = Ode_Coarse_Slicer(OdeTraj_thin, gridsize);

  arma::vec betaNs = betaTs(param,Ode_Coarse.col(0),x_r, x_i);

  arma::mat NewTraj = TransformTraj(Ode_Coarse, OriginTraj, FT_new);

  double coal_log_new = volz_loglik_nh_adj(init, NewTraj,betaNs, t_correct, x_i.subvec(2,3) ,addCoal);

  double incid_log_new = log_incidence(NewTraj, Incidence_data, incid_par, false,addIncid);

  double a = coal_log_new - coal_log + incid_log_new - incid_log + prior_proposal_offset;

  List Result;

  if(log(R::runif(0,1)) < a){

    Result["accept"] = true;
    Result["FT"] = FT_new;
    Result["Ode"] = Ode_Coarse;
    Result["betaN"] = betaNs;
    Result["coalLog"] = coal_log_new;
    Result["LatentTraj"] = NewTraj;
    Result["incidLog"] = incid_log_new;
  }else{
    Result["accept"] = false;
  }
  return Result;
}





//[[Rcpp::export()]]
List ESlice_par_General_JointData(arma::vec par_old, arma::vec incid_par, arma::vec t, arma::mat OriginTraj,List priorList,
                                  arma::vec x_r, arma::ivec x_i, List init, List Incidence_obs,int gridsize, arma::ivec ESS_vec,
                                  double coal_log = 0, double incid_log = 0, double t_correct = 0, std::string transP = "changepoint",
                                  std::string model = "SIR", std::string transX = "standard", bool volz = true,bool addCoal = true, bool addIncid = true){

  /*
  *
  *
  *  priorList: a list of priors for parameters for 1. I_0, 2. R_0, 3. gamma, 4, mu (no exist in SIR model) 5 hyper
  * all of them follows lognormal distribution
  *
  * ESS vec: vector that stores the index for each
  *
  *  test R code:
  *
  *  par_old = resres5MJP$par[125000,]
  A = ESlice_par_General(par_old, seq(0,100,length.out = 2001), resres5MJP$MCMC_obj$OriginTraj, priorList = list(a=c(1,1),b=c(0.7,0.2), c=c(2,1), d=c(3,0.2)),
                         resres5MJP$MCMC_setting$x_r, resres5MJP$MCMC_setting$x_i, resres5MJP$MCMC_setting$Init,50,ESS_vec = c(0,1,0,0,1,1,0),
  resres5MJP$MCMC_obj$coalLog, 90)

  A$par[c(1,2,4)] - par_old[c(1,2,4)]
  *
  */

  arma::ivec Index(2);

  int u_param = par_old.n_elem - 1; // length(par) - 1

  // extract prior distribution for each parameter

  //  arma::vec R0_pr = as<arma::vec>(priorList[1]);
  //  arma::vec gamma_pr = as<arma::vec>(priorList[2]);
  //  arma::vec hyper_pr = as<arma::vec>(priorList[3]);


  if(model == "SIR"){
    Index(0) = 0; Index(1) = 1;
  }else if(model == "SEIR"){
    Index(0) = 0; Index(1) = 2;
  }

  arma::vec new_par_raw = arma::randn(u_param, 1); // 47 * 1

  int N = t.n_elem / gridsize + 1;
  if(ESS_vec(6) == 1){
    arma::mat OriginTraj_new_raw = arma::randn(OriginTraj.n_rows, OriginTraj.n_cols);
  }

  arma::mat OriginTraj_new = OriginTraj;


  double u = R::runif(0,1);

  if(log(u) < -25){
    Rcout << "really small u" << endl;
    Rcout << u << endl;
  }

  double logy = coal_log + incid_log + log(u);
  double theta = 2 * pi;
  double theta_min = 0;
  double theta_max = 2 * pi;

  arma::vec par_new = par_old;
  List param_list;
  List FT_new;
  arma::vec betaN(N);
  arma::mat OdeTraj_new(N,3);
  arma::mat NewTraj(N,3);
  double coal_log_new;
  double incid_log_new;
  int i = 0;

  do{

    i += 1;

    // use the original values if no upates in 20 iterations
    if(i>20){
      //Rcout<<"theta = "<< theta <<endl;
      theta = 0;
      par_new = par_old;
      param_list = New_Param_List(par_new.subvec(2, x_i(0) + x_i(1) + 2), par_new.subvec(0,1), gridsize, t, x_r, x_i,
                                  transP, model, transX);
      FT_new = as<Rcpp::List>(param_list[0]);
      betaN = as<arma::vec>(param_list[2]);
      OdeTraj_new  = as<arma::mat>(param_list[1]);
      OriginTraj_new = OriginTraj;
      NewTraj = TransformTraj(OdeTraj_new, OriginTraj_new, FT_new);
      coal_log_new = volz_loglik_nh_adj(init, NewTraj,betaN,t_correct,Index ,addCoal);
      incid_log_new = log_incidence(NewTraj, Incidence_obs, incid_par,false,addIncid);
      //Rcout << "coal diff = " <<coal_log_new - coal_log << "incid diff = " << incid_log_new - incid_log << endl;
      break;
    }

    if(i == 1){ // first iteration
      theta = R::runif(0,2*pi);
      theta_min = theta - 2*pi;
      theta_max = theta;
    }else{

      if(theta < 0){
        theta_min = theta;
      }else{
        theta_max = theta;
      }
      theta = R::runif(theta_min,theta_max);
    }

    par_new = Param_Slice_update_all2(par_old, x_r, x_i, theta, new_par_raw, ESS_vec, priorList); // 49
    int error = 0;

    try{
      param_list = New_Param_List(par_new.subvec(2, x_i(0) + x_i(1) + 2), par_new.subvec(0,1), gridsize, t, x_r, x_i,
                                  transP, model, transX);
    }catch(...){
      error = 1;
    }

    if(error == 1) { // if error happens in parameters, go back
      coal_log_new = -10000000;
      incid_log_new = -10000000;
      continue;
    }

    FT_new = as<Rcpp::List>(param_list[0]);

    betaN = as<arma::vec>(param_list[2]);

    OdeTraj_new  = as<arma::mat>(param_list[1]);

    // Update OriginTraj if ESS_vec(6) == 1
    if(ESS_vec(6) == 1){
      OriginTraj_new = OriginTraj * cos(theta) + OriginTraj_new * sin(theta);
    }

    NewTraj = TransformTraj(OdeTraj_new, OriginTraj_new, FT_new);

    coal_log_new = volz_loglik_nh_adj(init, NewTraj, betaN, t_correct, Index ,addCoal);
    incid_log_new = log_incidence(NewTraj, Incidence_obs, incid_par,false,addIncid);
  }
  while(coal_log_new + incid_log_new <= logy);


  List result;
  double logOrigin = 0;
  for(int j = 0; j < OriginTraj_new.n_rows; j ++){
    for(int k = 0; k < OriginTraj_new.n_cols; k ++){
      logOrigin -= 0.5 * OriginTraj_new(j,k) * OriginTraj_new(j,k);
    }
  }
  result["betaN"] = betaN;
  result["FT"] = FT_new;
  result["OdeTraj"] = OdeTraj_new;
  result["par"] = par_new;
  result["LatentTraj"] = NewTraj;
  result["CoalLog"] = coal_log_new;
  result["IncidLog"] = incid_log_new;
  result["OriginTraj"] = OriginTraj_new;
  result["logOrigin"] = logOrigin;
  return result;
}

//[[Rcpp::export()]]
List ESlice_general_NC_joint(arma::mat f_cur, arma::mat OdeTraj, List FTs, arma::vec state, arma::vec incid_par,
                             List init, List Incidence, arma::vec betaN, double t_correct, double lambda=10,
                             double coal_log=-99999999, double IncidLog = -9999999, int gridsize = 100, bool volz = false, std::string model = "SIR",
                             std::string transX = "standard",bool addCoal = true, bool addIncid = true){
  // OdeTraj is the one with low resolution
  arma::ivec Index(2);
  if(model == "SIR"){
    Index(0) = 0; Index(1) = 1;
  }else if(model == "SEIR"){
    Index(0) = 0; Index(1) = 2;
  }
  int p = f_cur.n_cols;
  int incid_bound = as<int>(Incidence[2]);
  int init_bound = as<int>(init[9]) - 1;
  int idx_bound = (incid_bound > init_bound ? incid_bound + 1 : init_bound + 1);
  idx_bound = (idx_bound > f_cur.n_rows) ? f_cur.n_rows : idx_bound;
  arma::mat newTraj(f_cur.n_rows, f_cur.n_cols);
  double logy;

  arma::mat v_traj = arma::randn(f_cur.n_rows, f_cur.n_cols);
  double u = R::runif(0,1);
  if(log(u) < -25){
    Rcout << "really small u" << endl;
    Rcout << u << endl;
  }
  //  if(funname == "standard"){
  // logy = coal_loglik(init,LogTraj(f_cur),t_correct,lambda,gridsize) + log(u);
  logy = coal_log + IncidLog + log(u);

  double theta = R::runif(0,2 * pi);

  double theta_min = theta - 2*pi;
  double theta_max = theta;

  arma::mat f_prime = f_cur * cos(theta) + v_traj * sin(theta);
  newTraj = TransformTraj(OdeTraj,f_prime, FTs);
  int i = 0;
  double loglike;
  double IncidLog_new;
  if(volz){
    loglike = volz_loglik_nh_adj(init, newTraj,betaN,t_correct,Index ,addCoal);
  }else{
    loglike = coal_loglik3(init,LogTraj(newTraj),t_correct,lambda,Index(1), transX);
  }
  IncidLog_new = log_incidence(newTraj, Incidence, incid_par, false,addIncid);

  while(newTraj.submat(0, 1, idx_bound, p).min() <0 || (loglike + IncidLog_new) <= logy){
    // shrink the bracket
    i += 1;
    if(i>20){
      newTraj = TransformTraj(OdeTraj,f_cur, FTs);
      loglike = volz_loglik_nh_adj(init, newTraj,betaN,t_correct,Index , addCoal);
      IncidLog_new = log_incidence(newTraj, Incidence, incid_par,false,addIncid);
      f_prime = f_cur;
      //Rcout << "theta = "<< theta << endl;
      //Rcout << "coal diff0" <<loglike - coal_log << "incid diff0" << IncidLog_new - IncidLog << endl;
      break;
    }
    if(theta < 0){
      theta_min = theta;
    }else{
      theta_max = theta;
    }
    theta = R::runif(theta_min,theta_max);
    f_prime = f_cur * cos(theta) + v_traj * sin(theta);
    newTraj = TransformTraj(OdeTraj,f_prime, FTs);

    if(volz){
      loglike = volz_loglik_nh_adj(init, newTraj, betaN, t_correct, Index , addCoal);
    }else{
      loglike = coal_loglik3(init,LogTraj(newTraj),t_correct,lambda,Index(1), transX);
    }
    IncidLog_new = log_incidence(newTraj, Incidence, incid_par, false,addIncid);
    if(isnan(loglike)){
      Rcout << "NaN in Slice Sampler " << endl;
    }
    if(isnan(IncidLog_new)){
      Rcout << "NaN in Slice Sampler " << endl;
    }
  }
  List Res;
  double logOrigin = 0;
  for(int j = 0; j < f_cur.n_rows - 1; j ++){
    for(int k = 0; k < p; k ++){
      logOrigin -= 0.5 * f_prime(j,k) * f_prime(j,k);
    }
  }
  Res["LatentTraj"] = newTraj;
  Res["OriginTraj"] = f_prime;
  Res["logOrigin"] = logOrigin;
  Res["CoalLog"] = loglike;
  Res["IncidLog"] = IncidLog_new;
  return Res;
}




/*
 * functions for preferential sampling
 *
 *
 */



//[[Rcpp::export()]]
List Update_Param_Pref(arma::vec param, arma::vec initial, arma::vec pref_par, arma::vec t, arma::mat OriginTraj,
                       arma::vec x_r, arma::ivec x_i, List init, List SampleInit, int gridsize,
                       double coal_log=0, double pref_log = 0, double prior_proposal_offset = 0, double t_correct = 0, std::string transP = "changepoint",
                       std::string model = "SIR", std::string transX = "standard", bool volz = true, bool addCoal = true, bool addPref = true){

  arma::mat OdeTraj_thin = ODE_rk45(initial,t, param,
                                    x_r, x_i, transP, model, transX);

  List FT_new = KF_param_chol(OdeTraj_thin, param, gridsize, x_r, x_i, transP, model, transX);

  arma::mat Ode_Coarse = Ode_Coarse_Slicer(OdeTraj_thin, gridsize);

  arma::vec betaNs = betaTs(param,Ode_Coarse.col(0),x_r, x_i);

  arma::mat NewTraj = TransformTraj(Ode_Coarse, OriginTraj, FT_new);

  double coal_log_new = volz_loglik_nh_adj(init, NewTraj,betaNs, t_correct, x_i.subvec(2,3) ,addCoal);

  double pref_log_new = log_prev(NewTraj, SampleInit, pref_par, true, addPref);

  double a = coal_log_new - coal_log + pref_log_new - pref_log + prior_proposal_offset;

  List Result;

  if(log(R::runif(0,1)) < a){

    Result["accept"] = true;
    Result["FT"] = FT_new;
    Result["Ode"] = Ode_Coarse;
    Result["betaN"] = betaNs;
    Result["coalLog"] = coal_log_new;
    Result["LatentTraj"] = NewTraj;
    Result["PrefLog"] = pref_log_new;
  }else{
    Result["accept"] = false;
  }
  return Result;
}


//[[Rcpp::export()]]
List ESlice_par_General_pref(arma::vec par_old, arma::vec pref_par, arma::vec t, arma::mat OriginTraj,List priorList,
                             arma::vec x_r, arma::ivec x_i, List init, List SampleInit,int gridsize, arma::ivec ESS_vec,
                             double coal_log = 0, double pref_log = 0, double t_correct = 0, std::string transP = "changepoint",
                             std::string model = "SIR", std::string transX = "standard", bool volz = true,bool addCoal = true, bool addPref = true){

  /*
  *
  *
  *  priorList: a list of priors for parameters for 1. I_0, 2. R_0, 3. gamma, 4, mu (no exist in SIR model) 5 hyper
  * all of them follows lognormal distribution
  *
  * ESS vec: vector that stores the index for each
  *
  *  test R code:
  *
  *  par_old = resres5MJP$par[125000,]
  A = ESlice_par_General(par_old, seq(0,100,length.out = 2001), resres5MJP$MCMC_obj$OriginTraj, priorList = list(a=c(1,1),b=c(0.7,0.2), c=c(2,1), d=c(3,0.2)),
                         resres5MJP$MCMC_setting$x_r, resres5MJP$MCMC_setting$x_i, resres5MJP$MCMC_setting$Init,50,ESS_vec = c(0,1,0,0,1,1,0),
  resres5MJP$MCMC_obj$coalLog, 90)

  A$par[c(1,2,4)] - par_old[c(1,2,4)]
  *
  */

  arma::ivec Index(2);

  int u_param = par_old.n_elem - 1; // length(par) - 1

  // extract prior distribution for each parameter

  //  arma::vec R0_pr = as<arma::vec>(priorList[1]);
  //  arma::vec gamma_pr = as<arma::vec>(priorList[2]);
  //  arma::vec hyper_pr = as<arma::vec>(priorList[3]);


  if(model == "SIR"){
    Index(0) = 0; Index(1) = 1;
  }else if(model == "SEIR"){
    Index(0) = 0; Index(1) = 2;
  }

  arma::vec new_par_raw = arma::randn(u_param, 1); // 47 * 1

  int N = t.n_elem / gridsize + 1;
  if(ESS_vec(6) == 1){
    arma::mat OriginTraj_new_raw = arma::randn(OriginTraj.n_rows, OriginTraj.n_cols);
  }

  arma::mat OriginTraj_new = OriginTraj;


  double u = R::runif(0,1);

  if(log(u) < -25){
    Rcout << "really small u" << endl;
    Rcout << u << endl;
  }

  double logy = coal_log + pref_log + log(u);
  double theta = 2 * pi;
  double theta_min = 0;
  double theta_max = 2 * pi;

  arma::vec par_new = par_old;
  List param_list;
  List FT_new;
  arma::vec betaN(N);
  arma::mat OdeTraj_new(N,3);
  arma::mat NewTraj(N,3);
  double coal_log_new;
  double pref_log_new;
  int i = 0;

  do{

    i += 1;

    // use the original values if no upates in 20 iterations
    if(i>20){
      //Rcout<<"theta = "<< theta <<endl;
      theta = 0;
      par_new = par_old;
      param_list = New_Param_List(par_new.subvec(2, x_i(0) + x_i(1) + 2), par_new.subvec(0,1), gridsize, t, x_r, x_i,
                                  transP, model, transX);
      FT_new = as<Rcpp::List>(param_list[0]);
      betaN = as<arma::vec>(param_list[2]);
      OdeTraj_new  = as<arma::mat>(param_list[1]);
      OriginTraj_new = OriginTraj;
      NewTraj = TransformTraj(OdeTraj_new, OriginTraj_new, FT_new);
      coal_log_new = volz_loglik_nh_adj(init, NewTraj,betaN,t_correct,Index ,addCoal);
      pref_log_new = log_prev(NewTraj, SampleInit, pref_par,true,addPref);
      //Rcout << "coal diff = " <<coal_log_new - coal_log << "Pref diff = " << pref_log_new - pref_log << endl;
      break;
    }

    if(i == 1){ // first iteration
      theta = R::runif(0,2*pi);
      theta_min = theta - 2*pi;
      theta_max = theta;
    }else{

      if(theta < 0){
        theta_min = theta;
      }else{
        theta_max = theta;
      }
      theta = R::runif(theta_min,theta_max);
    }

    par_new = Param_Slice_update_all2(par_old, x_r, x_i, theta, new_par_raw, ESS_vec, priorList); // 49
    int error = 0;

    try{
      param_list = New_Param_List(par_new.subvec(2, x_i(0) + x_i(1) + 2), par_new.subvec(0,1), gridsize, t, x_r, x_i,
                                  transP, model, transX);
    }catch(...){
      error = 1;
    }

    if(error == 1) { // if error happens in parameters, go back
      coal_log_new = -10000000;
      pref_log_new = -10000000;
      continue;
    }

    FT_new = as<Rcpp::List>(param_list[0]);

    betaN = as<arma::vec>(param_list[2]);

    OdeTraj_new  = as<arma::mat>(param_list[1]);

    // Update OriginTraj if ESS_vec(6) == 1
    if(ESS_vec(6) == 1){
      OriginTraj_new = OriginTraj * cos(theta) + OriginTraj_new * sin(theta);
    }

    NewTraj = TransformTraj(OdeTraj_new, OriginTraj_new, FT_new);

    coal_log_new = volz_loglik_nh_adj(init, NewTraj, betaN, t_correct, Index ,addCoal);
    pref_log_new = log_prev(NewTraj, SampleInit, pref_par,true,addPref);
  }
  while(coal_log_new + pref_log_new <= logy);


  List result;
  double logOrigin = 0;
  for(int j = 0; j < OriginTraj_new.n_rows; j ++){
    for(int k = 0; k < OriginTraj_new.n_cols; k ++){
      logOrigin -= 0.5 * OriginTraj_new(j,k) * OriginTraj_new(j,k);
    }
  }
  result["betaN"] = betaN;
  result["FT"] = FT_new;
  result["OdeTraj"] = OdeTraj_new;
  result["par"] = par_new;
  result["LatentTraj"] = NewTraj;
  result["CoalLog"] = coal_log_new;
  result["PrefLog"] = pref_log_new;
  result["OriginTraj"] = OriginTraj_new;
  result["logOrigin"] = logOrigin;
  return result;
}

//[[Rcpp::export()]]
List ESlice_general_NC_pref(arma::mat f_cur, arma::mat OdeTraj, List FTs, arma::vec state, arma::vec pref_par,
                            List init, List SampleInit, arma::vec betaN, double t_correct, double lambda=10,
                            double coal_log = -99999999, double PrefLog = -9999999, int gridsize = 100, bool volz = false, std::string model = "SIR",
                            std::string transX = "standard",bool addCoal = true, bool addPref = true){

  // OdeTraj is the one with low resolution

  arma::ivec Index(2);
  if(model == "SIR"){
    Index(0) = 0; Index(1) = 1;
  }else if(model == "SEIR"){
    Index(0) = 0; Index(1) = 2;
  }
  int p = f_cur.n_cols;
  int Pref_bound = as<int>(SampleInit[3]);
  int init_bound = as<int>(init[9]) - 1;
  int idx_bound = (Pref_bound > init_bound ? Pref_bound + 1 : init_bound + 1);
  idx_bound = (idx_bound > f_cur.n_rows) ? f_cur.n_rows : idx_bound;
  arma::mat newTraj(f_cur.n_rows, f_cur.n_cols);
  double logy;

  arma::mat v_traj = arma::randn(f_cur.n_rows, f_cur.n_cols);
  double u = R::runif(0,1);
  if(log(u) < -25){
    Rcout << "really small u" << endl;
    Rcout << u << endl;
  }
  //  if(funname == "standard"){
  // logy = coal_loglik(init,LogTraj(f_cur),t_correct,lambda,gridsize) + log(u);
  logy = coal_log + PrefLog + log(u);

  double theta = R::runif(0,2 * pi);

  double theta_min = theta - 2*pi;
  double theta_max = theta;

  arma::mat f_prime = f_cur * cos(theta) + v_traj * sin(theta);
  newTraj = TransformTraj(OdeTraj,f_prime, FTs);
  int i = 0;
  double loglike;
  double PrefLog_new;
  if(volz){
    loglike = volz_loglik_nh_adj(init, newTraj,betaN,t_correct,Index ,addCoal);
  }else{
    loglike = coal_loglik3(init,LogTraj(newTraj),t_correct,lambda,Index(1), transX);
  }
  PrefLog_new = log_prev(newTraj, SampleInit, pref_par, true,addPref);

  while(newTraj.submat(0, 1, idx_bound, p).min() <0 || (loglike + PrefLog_new) <= logy){
    // shrink the bracket
    i += 1;
    if(i>20){
      newTraj = TransformTraj(OdeTraj,f_cur, FTs);
      loglike = volz_loglik_nh_adj(init, newTraj,betaN,t_correct,Index , addCoal);
      PrefLog_new = log_prev(newTraj, SampleInit, pref_par,true,addPref);
      f_prime = f_cur;
      //Rcout << "theta = "<< theta << endl;
      //Rcout << "coal diff0" <<loglike - coal_log << "Pref diff0" << PrefLog_new - PrefLog << endl;
      break;
    }
    if(theta < 0){
      theta_min = theta;
    }else{
      theta_max = theta;
    }
    theta = R::runif(theta_min,theta_max);
    f_prime = f_cur * cos(theta) + v_traj * sin(theta);
    newTraj = TransformTraj(OdeTraj,f_prime, FTs);

    if(volz){
      loglike = volz_loglik_nh_adj(init, newTraj, betaN, t_correct, Index , addCoal);
    }else{
      loglike = coal_loglik3(init,LogTraj(newTraj),t_correct,lambda,Index(1), transX);
    }
    PrefLog_new = log_prev(newTraj, SampleInit, pref_par, true,addPref);
    if(isnan(loglike)){
      Rcout << "NaN in Slice Sampler " << endl;
    }
    if(isnan(PrefLog_new)){
      Rcout << "NaN in Slice Sampler " << endl;
    }
  }
  List Res;
  double logOrigin = 0;
  for(int j = 0; j < f_cur.n_rows - 1; j ++){
    for(int k = 0; k < p; k ++){
      logOrigin -= 0.5 * f_prime(j,k) * f_prime(j,k);
    }
  }
  Res["LatentTraj"] = newTraj;
  Res["OriginTraj"] = f_prime;
  Res["logOrigin"] = logOrigin;
  Res["CoalLog"] = loglike;
  Res["PrefLog"] = PrefLog_new;
  return Res;
}
