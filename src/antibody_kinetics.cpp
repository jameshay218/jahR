#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector kinetics_power_function_cpp(NumericVector times, NumericVector pars) {
  int n_times = times.size();
  NumericVector titers(n_times);
  double y0 = pars["y0"];
  double peak = pars["peak"];
  double k = pars["k"];
  double r = pars["r"];
  double v = pars["v"];
  double t_peak = pars["t_peak"];
  
  peak = (peak - 1.0) * exp(-y0*k) + 1.0;
  double mu = (1.0/t_peak) * log(peak);
    
  double peak_reached = peak*y0;
  
  for(int t = 0; t < n_times; t++){
    if(times[t] <= t_peak){
      titers[t] =  y0*exp(mu*times[t]);
    } else {
      titers[t] = peak_reached * pow((1.0 + v*(r - 1.0)*(pow(peak_reached,(r-1.0)))*(times[t]-t_peak)),(-1/(r-1)));
    }
  }
  return titers;
}

