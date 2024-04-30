## Can use the `with` function and named vectors to make this shorter
SIR_odes_with <- function(t, States, pars){
  with(as.list(c(States,pars)),{
    
    dS <- -beta*S*I
    dI <- beta*S*I - gamma*I
    dR <- gamma*I
    inc <- beta*S*I
    
    return(list(c(dS, dI, dR, inc)))
  })
}

#' Simple SIR model implementation
#' 
#' @param pars named vector of model parameters, give R0 and gamma
#' @param times vector of times to solve model over
#' @param initial_states named vector of initial states, giving values for S, I and R. If the sum is more than 1, uses frequency dependent model
#' @return a data frame with the solved SIR model, including compartment sizes over time, incidence and growth rate of incidence
#' @export
sir_model <- function(pars=c(R0=1.5,gamma=1/5),times=seq(0,365,by=0.1),initial_states=c(S=1-0.00001, I=0.00001, R=0)){
  S0 <- initial_states["S"]
  I0 <- initial_states["I"]
  R0 <- initial_states["R"]
  initial_states <- c(initial_states,"inc"=0)
  basic_repro_number <- pars["R0"]
  gamma <- pars["gamma"]
  if(sum(initial_states) > 1){
    print("Frequency dependent")
    beta <- basic_repro_number*gamma/sum(S0, I0, R0)
  } else {
    print("Density dependent")
    
    beta <- basic_repro_number*gamma ## Get transmission rate from above pars
  }
  pars <- c(beta=beta,gamma=gamma)
  names(pars) <- c("beta","gamma")
  
  results <- as.data.frame(deSolve::ode(y=initial_states, times=times, func=SIR_odes_with, parms= pars))
  results <- results %>% mutate(inc = inc - lag(inc, 1))
  results <- results %>% mutate(gr = log(inc/lag(inc,1)))
  ## Plot the results
  plot(results$I~results$time,type='l',col="red",
       ylim=c(0,sum(initial_states)),xlab="Time (days)", ylab="Per capita",
       main="SIR dynamics")
  lines(results$R~results$time,col="blue")
  lines(results$S~results$time,col="green")
  legend("topright", title="Key",
         legend=c("Susceptible","Infected","Recovered"),
         col=c("green","red","blue"),
         lty=c(1,1,1))
  
  ## What is the final cumulative incidence?
  print(paste0("Final size: ", initial_states["S"]-min(results$S)))
  
  ## Check against theory
  simeq <- function(A, R0){
    y <- A - (1-exp(-R0*A))
  }
  print(paste0("Nleqslv final size: ", nleqslv::nleqslv(x=0.5,fn=simeq,R0=1.5)$x))
  return(results)
}
