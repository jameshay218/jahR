#' Power function boosting and waning model
#' 
#' Solve the model described in Teunis et al. 2016 (Epidemics) and used in Aiemjoy et al. 2022 (Lancet Microbe). There is an addition of a titer-dependent boosting term as used in Ranjeva et al. 2019 (Nature Communications). 
#' @param times vector of times to solve model over
#' @param pars vector of named model parameters
#' @details Model parameters:
#' \itemize{
#'  \item{"t_peak"}{time of peak antibody level post exposure}
#'  \item{"peak"}{peak antibody level on natural scale}
#'  \item{"k"}{antibody level dependent boosting term}
#'  \item{"y0"}{starting antibody level on natural scale -- Note this needs to start at 1}
#'  \item{"v"}{antibody decay rate}
#'  \item{"r"}{shape of antibody decay function}
#'  }
#' @family kinetics_models
#' @return vector of solved antibody levels
#' @export
#' @examples
#' N <- 25
#' r <- rlnorm(N, log(1.3),0.01) ## Shape parameter
#' v <- rexp(N,500) ## Decay rate
#' y0 <- 1## Start titre
#' k <- 0 ## Titre dependency
#' peak <- 10^(rnorm(N,5, 1.5)) ## Peak titre
#' 
#' times <- seq(0,365,by=0.1)
#' ys <- matrix(nrow=N,ncol=length(times))
#' for(i in 1:N){
#'   ys[i,] <- kinetics_power_function(times, c(t_peak=21,peak=peak[i],r=r[i],v=v[i],y0=1,k=0))
#' }
#' ys <- log10(ys)
#' 
#' plot(times,(ys[1,]),type='l',ylim=c(0,10))
#' for(i in 2:N){
#'   lines(times,(ys[i,]))
#' }
kinetics_power_function <- function(times, pars){
  t1 <- times <= pars["t_peak"]
  t2 <- times > pars["t_peak"]
  y <- numeric(length(times))
  peak <- pars["peak"]
  k <- pars["k"]
  y0 <- pars["y0"]
  r <- pars["r"]
  v <- pars["v"]
  t_peak <- pars["t_peak"]
  ## Adjust to titre ceiling effect
  peak <- (peak - 1) * exp(-y0*k) + 1
  mu <- (1/t_peak) * log(peak)
  
  peak_reached <- peak*y0
  y[t1] <- y0*exp(mu*times[t1])
  y[t2] <- peak_reached * (1 + v*(r - 1)*(peak_reached^(r-1))*(times[t2]-t_peak))^(-1/(r-1))
  y
}



#' Internal boost function for Neil's simple boosting/waning model
boost_function <- function(times, pars){
  a <- pars["a"]
  b <- pars["b"]
  t_peak <- pars["t_peak"]
  y <- ((a+b)/b) * exp(a*(times-t_peak))
  y
}
#' Internal wane function for Neil's simple boosting/waning model, with added asymptote
#' @export
wane_function_asymptote <- function(times, pars){
  a <- pars["a"]
  b <- pars["b"]
  c <- pars["c"]
  t_peak <- pars["t_peak"]
  y <- ((a+b)/a - c) * exp(-b*(times-t_peak)) + c
  y
}
#' Internal wane function for Neil's simple boosting/waning model
#' @export
wane_function <- function(times, pars){
  a <- pars["a"]
  b <- pars["b"]
  t_peak <- pars["t_peak"]
  y <- ((a+b)/a) * exp(-b*(times-t_peak))
  y
}

#' Simple boosting and waning kinetics
#' 
#' Implements the biphasic boosting and waning model from Singanayagam et al. 2021 (Lancet ID). A phenomonelogical assuming exponential growth and decline, where the weighting of the two components shifts from growth to decline, swapping over at the peak time. Can have an asymptote or not in the waning rate. NOTE: when using the asymptote waning function, the function *does not* peak at t_peak, but slightly after.
#' @inheritParams kinetics_power_function
#' @param boost pointer to the function used to calculate boosting (defaults to boost_function)
#' @param wane pointer to the function used to calculate waning (defaults to wane_function_asymptote)
#' @details Model parameters:
#' \itemize{
#'  \item{"t_peak"}{time of peak antibody level post exposure}
#'  \item{"ymax"}{peak antibody level on natural scale}
#'  \item{"a"}{boosting rate}
#'  \item{"b"}{waning rate}
#'  \item{"c"}{asymptote of waning (can set to 0, or use wane_function instead of wane_function_asymptote to decline to 0)}
#'  }
#' @family kinetics_models
#' @return vector of antibody levels
#' @export
#' @examples
#' pars<- c(a=0.3,b=0.05,c=0.75,ymax=100,t_peak=18)
#' times <- seq(0,100,by=1)
#' y <- kinetics_simple(times,pars)
#' plot(times, y, type='l')
kinetics_simple <- function(times, pars, boost=boost_function, wane=wane_function_asymptote){
  y_boost <- boost(times, pars)
  y_wane <- wane(times, pars)
  pars["ymax"]*y_boost*y_wane/(y_boost + y_wane)
}


#' Internal boosting function for Ranjeva et al.
#' @export
f_rise <- function(times, pars){
  y0 <- pars["y0"]
  peak <- pars["peak"]
  k <- pars["k"]
  peak <- exp(-y0*pars["k"])*peak
  a <- (1/pars["t_peak"]) * log(peak + 1)
  exp(a*(times)) + y0 - 1
}
#' Internal waning function for Ranjeva et al.
#' @export
f_wane <- function(times, pars){
  y_peak <- f_rise(pars["t_peak"], pars)
  y_baseline <- pars["y0"] + pars["frac"]*pars["peak"]*exp(-pars["y0"]*pars["k"])
  w <- pars["w"]
  (y_peak - y_baseline)*exp(-w*(times-pars["t_peak"]))
}

#' Boosting and waning model with ceiling effect and asymptoting boost
#' 
#' Implements the antibody kinetics model from Ranjeva et al. 2019 (Nature Communications). Key features of this model are a short-term, titer-dependent boost that asymptotes towards a peak value rather than continuing to grow exponentially. The waning phase is biphasic -- an initial proportion of the boost is lost and the asymptotes to a steady state.
#' @inheritParams kinetics_simple
#' @param boost pointer to the antibody boosting function (defaults to f_rise)
#' @param wane pointer to the antibody waning function (defaults to f_wane)
#' #' @details Model parameters:
#' \itemize{
#'  \item{"t_peak"}{time of peak antibody level post exposure}
#'  \item{"peak"}{peak antibody level on natural scale}
#'  \item{"k"}{antibody level dependent boosting term}
#'  \item{"y0"}{starting antibody level on natural scale}
#'  \item{"w"}{antibody decay rate}
#'  \item{"frac"}{proportion of the intial boost lost during short-term waning}
#'  }
#' @family kinetics_models
#' @export
#' @return vector of antibody levels
#' @examples
#' pars <- c(y0=1,peak=80,k=0,t_peak=28,w=0.008,frac=0.5)
#' times <- seq(0,1000,by=0.1)
#' y <- kinetics_ranjeva(times,pars)
#' plot(times, y, type='l')
kinetics_ranjeva <- function(times, pars, boost=f_rise, wane=f_wane){
  y0 <- pars["y0"]
  y_baseline <- y0 + pars["frac"]*pars["peak"]*exp(-y0*pars["k"])
  y <- numeric(length(times))
  t1 <- times < pars["t_peak"]
  t2 <- times >= pars["t_peak"]
  y[t1] <- boost(times[t1], pars)
  y[t2] <- y_baseline + wane(times[t2], pars)
  y
}


#' Boosting and waning model based on within-host ODEs
#' 
#' Implements the analytical solution to the ODE model described in Pelleau et al. 2021 (JID).
#' @inheritParams kinetics_power_function
#' @details Model parameters:
#' \itemize{
#'  \item{"A"}{initial antibody level}
#'  \item{"beta"}{growth rate of antibody secreting cells}
#'  \item{"rho"}{proportion of B cells which are short-lived}
#'  \item{"delta"}{time at which antibody levels start increasing}
#'  \item{"c_l"}{decay rate of long-lived B cells}
#'  \item{"c_s"}{decay rate of short-lived B cells}
#'  \item{"r"}{antibody decay rate}
#'  }
#' @family kinetics_models
#' @export
#' @return vector of antibody levels
#' @examples
#' pars <- c("A"=0,"beta"=100,"rho"=0.25,"c_s"=0.15,"c_l"=0.0075,"delta"=0,"r"=0.001)
#' times <- seq(0,1000,by=0.1)
#' y <- kinetics_ode(times,pars)
#' plot(times, y, type='l')
kinetics_ode <- function(times, pars){
  A <- pars["A"]
  beta <- pars["beta"]
  rho <- pars["rho"]
  c_s <- pars["c_s"]
  c_l <- pars["c_l"]
  delta <- pars["delta"]
  r <- pars["r"]
  
  short_term <- exp(-c_s*(times-delta))- exp(-r*(times-delta))
  long_term <- exp(-c_l*(times-delta))- exp(-r*(times-delta))
  
  y <- A + beta * ( rho * short_term / (r-c_s) + (1-rho)*long_term/(r-c_l))
  y
}


#' Gamma function kinetics
#' 
#' Very simple scaled gamma function to represent kinetics. Used in Zhao et al. 2018, but not particular interpretable otherwise.
#' @inheritParams kinetics_power_function
#' @details Model parameters:
#' \itemize{
#'  \item{"y0"}{initial antibody level}
#'  \item{"peak"}{parameter to scale function in y-axis}
#'  \item{"shape"}{gamma distribution shape}
#'  \item{"scale"}{gamma distribution scale}
#'  }
#' @family kinetics_models
#' @export
#' @return vector of antibody levels
#' @examples
#' pars <- c("y0"=0,"peak"=300,"shape"=1.2,"scale"=50)
#' times <- seq(0,1000,by=0.1)
#' y <- kinetics_gamma(times,pars)
#' plot(times, y, type='l')
kinetics_gamma <- function(times, pars){
  pars["y0"] + pars["peak"]*dgamma(times, shape=pars["shape"],scale=pars["scale"])
}


#' Piecewise linear model of boosting and biphasic waning
#' 
#' Simple piecewise linear model of linear antibody boosting on log scale followed by linear biphasic waning. Model as in Hay et al. 2019 (PLOS Comp Biol).
#' @inheritParams kinetics_power_function
#' @details Model parameters:
#' \itemize{
#'  \item{"t_i"}{Timing of infection}
#'  \item{"tp"}{Time of peak antibody level post infection}
#'  \item{"mu"}{Magnitude of boost}
#'  \item{"dp"}{Proportion of boost lost in short-term waning}
#'  \item{"ts"}{Duration of short term waning phase}
#'  \item{"m"}{Long-term antibody waning rate}
#'  \item{"lower_bound"}{Lower bound of allowable antibody levels}
#'  \item{"y0"}{True starting titer level}
#'  }
#' @family kinetics_models
#' @export
#' @return vector of antibody levels
#' @examples
#' pars <- c("tp"=21,"mu"=8,"dp"=0.2,"ts"=50,"m"=0.0075,"lower_bound"=0, "y0"=0)
#' times <- seq(0,100,by=0.1)
#' y <- kinetics_piecewise(times,pars)
#' plot(times, y, type='l')
kinetics_piecewise <- function(times, pars){
  t_i <- 0 # pars["t_i"]
  tp <- pars["tp"]
  mu <- pars["mu"]
  dp <- pars["dp"]
  ts <- pars["ts"]
  m <- pars["m"]
  lower_bound <- pars["lower_bound"]
  eff_y0 <- pars["y0"]
  
  y <- numeric(length(times))
  i <- 1
  ## Loops through all times and calculate titre based on time relative to time of infection
  for(t in times){
    tmp <- 0
    if(t <= t_i) tmp = 0
    else if(t > t_i & t <= (t_i + tp)) tmp = (mu/tp)*(t-t_i)
    else if(t > (tp+t_i) & t <=(ts + t_i+tp)) tmp = ((-(dp*mu)/ts)*(t) + ((mu*dp)/ts)*(t_i+tp) + mu)
    else tmp = (-m*(t)+m*(t_i+tp+ts)+(1-dp)*mu)
    y[i] <- tmp + eff_y0
    if(y[i] < lower_bound) y[i] = lower_bound
    i <- i + 1
  }
  y
}

#' @export
create_kinetics_posterior_function <- function(parTab, data, PRIOR_FUNC,model_func, use_log=FALSE, ver="posterior",
                                               random_effects=TRUE,titer_range=NULL, solve_likelihood=TRUE){
  par_names <-parTab$names
  times <- data$time
  titers <- data$titer
  model_titers <- numeric(length(titers))
  
  unique_indivs <- unique(data$i)
  model_par_indices <- vector(mode="list",length=length(unique_indivs))  
  data_indices <- vector(mode="list",length=length(unique_indivs))  
  for(indiv in unique_indivs){
    model_par_indices[[indiv]] <- which(parTab$i == indiv)
    data_indices[[indiv]] <- which(data$i == indiv)
  }
  
  par_means <- which(par_names %like% "mean")
  par_vars <- which(par_names %like% "var")

  f <- function(pars){
    names(pars) <- par_names
    if(solve_likelihood | ver == "model"){
      for(indiv in unique_indivs){
        tmp_pars <- pars[model_par_indices[[indiv]]]
        
        if(random_effects){
            tmp_pars <- exp(pars[par_means] + pars[par_vars]*tmp_pars)
        } else {
            tmp_pars <- exp(pars[par_means])
        }
        model_titers[data_indices[[indiv]]] <- model_func(times[data_indices[[indiv]]], tmp_pars)
      }
      if(use_log){
        model_titers <- log10(model_titers)
      }
      if(ver == "model") return(model_titers)
      
   
      if(!is.null(titer_range)){
        lik <- serosolver::likelihood_func_fast_continuous(c("obs_sd"=pars["obs_sd"],"MIN_TITRE"=titer_range[1],"MAX_TITRE"=titer_range[2]), titers, model_titers)
      } else {
        lik <- sum(dnorm(titers, model_titers, sd=pars["obs_sd"],log=TRUE))
      }
    } else {
      lik <- 0
    }

    if(!is.null(PRIOR_FUNC)) lik <- lik + PRIOR_FUNC(pars)
    lik
  }
  return(f)
}

#' @export
create_multi_kinetics_posterior_function <- function(parTab, data, PRIOR_FUNC,model_func, use_log=FALSE, ver="posterior",
                                               random_effects=TRUE,titer_range=NULL,solve_likelihood=TRUE){
  par_names <- parTab$names
  times <- data$time
  titers <- data$titer
  model_titers <- numeric(length(titers))
  
  unique_indivs <- unique(data$i)
  model_par_indices <- vector(mode="list",length=length(unique_indivs))  
  data_indices <- vector(mode="list",length=length(unique_indivs))  
  unique_exposures <- vector(mode="list",length=length(unique_indivs))  
  exposure_table <- unique(parTab[,c("i","x","tinf")])
  
  for(indiv in unique_indivs){
    ## How many exposures does this individual have?
    exposure_times <- exposure_table[exposure_table$i == indiv, "tinf"]
    n_exposures <- length(exposure_times)
    
    ## Vector of exposure times
    unique_exposures[[indiv]] <- numeric(n_exposures)
    ## Need a set of parameters for each exposures
    model_par_indices[[indiv]] <- vector(mode="list",length=n_exposures)  
    ## Need a set of times to solve model over for each exposure
    data_indices[[indiv]] <- vector(mode="list",length=n_exposures)  
    
    exposure_times_tmp <- c(exposure_times, max(data[data$i == indiv,"time"]) + 1)
    
    ## For each exposure, record its time, which parameter indices correspond to it, and 
    ## which time indices in the data data.frame should be solved for this exposure
    for(exposure in seq_along(exposure_times)){
      unique_exposures[[indiv]][exposure] <- exposure_times[exposure]
      model_par_indices[[indiv]][[exposure]] <- which(parTab$i == indiv & parTab$x == exposure)
      data_indices[[indiv]][[exposure]] <- which(data$i == indiv & data$time >= exposure_times_tmp[exposure] & 
                                                   data$time < exposure_times_tmp[exposure + 1])
    }
  }
  
  par_means <- which(par_names %like% "mean")
  par_vars <- which(par_names %like% "var")
  
  f <- function(pars){
    names(pars) <- par_names
    if(solve_likelihood | ver == "model"){
      for(indiv in unique_indivs){
        ## First exposure
        ## Get parameter indices for exposure 1
        tmp_pars <- pars[model_par_indices[[indiv]][[1]]]
        tinf <- unique_exposures[[indiv]][[1]]
        par_names_indiv <- names(tmp_pars)
        ## Transform to natural scale
        if(random_effects){
          tmp_pars <- exp(pars[par_means] + pars[par_vars]*tmp_pars)
        } else {
          tmp_pars <- exp(pars[par_means])
        }
        names(tmp_pars) <- par_names_indiv
        ## Solve model for times governed by the first exposure
        model_titers[data_indices[[indiv]][[1]]] <- model_func(times[data_indices[[indiv]][[1]]]-tinf, tmp_pars)
        
        ## If there are subsequent exposures, solve y0 for the next exposure
        if(length(unique_exposures[[indiv]]) > 1){
          for(j in 2:length(unique_exposures[[indiv]])){
            tinf_prev <- tinf
            tinf <- unique_exposures[[indiv]][[j]]
            
            ## Find titer predicted from end of previous exposure's parameters
            y0 <-  model_func(unique_exposures[[indiv]][j]-tinf_prev, tmp_pars)
            
            ## Extract parameters for this exposure ID
            tmp_pars <- pars[model_par_indices[[indiv]][[j]]]
            par_names_indiv <- names(tmp_pars)
            
            ## Transform parameters
            if(random_effects){
              tmp_pars <- exp(pars[par_means] + pars[par_vars]*tmp_pars)
            } else {
              tmp_pars <- exp(pars[par_means])
            }
            names(tmp_pars) <- par_names_indiv
            
            ## Set y0 for this exposure
            tmp_pars["y0"] <- y0
            
            ## Solve model for this exposure's times
            model_titers[data_indices[[indiv]][[j]]] <- model_func(times[data_indices[[indiv]][[j]]]-tinf, tmp_pars)
          }
        }
      }
      
      if(use_log){
        model_titers <- log10(model_titers)
      }
      if(ver == "model") return(model_titers)
      if(!is.null(titer_range)){
          lik <- sum(serosolver::likelihood_func_fast_continuous(c("error"=unname(pars["obs_sd"]),"MIN_TITRE"=titer_range[1],"MAX_TITRE"=titer_range[2]), titers, model_titers))
        } else {
          lik <- sum(dnorm(titers, model_titers, sd=pars["obs_sd"],log=TRUE))
        }
    } else {
      lik <- 0
    }
    if(!is.null(PRIOR_FUNC)) lik <- lik + PRIOR_FUNC(pars)
    lik
  }
  return(f)
}


## Get prediction intervals
#' @export
get_kinetics_prediction_intervals <- function(model_func, n_samp, indivs=NULL, chain, dat, expand_times=TRUE, times=NULL, use_log=TRUE, add_noise=FALSE,titer_range=NULL){
  unique_samps <- unique(chain$sampno)
  samps <- sample(unique_samps, n_samp)
  
  
  if(is.null(indivs)){
    unique_indivs <- unique(dat$i)
  } else {
    dat <- dat[dat$i %in% indivs,]
    unique_indivs <- indivs
  }
  
  
  if(expand_times){
    solved_model1 <- expand_grid(i=unique_indivs, time=times,titer=0)
    times_use <- times
  } else {
    solved_model1 <- dat[,c("i","time","titer")]
    times_use <- dat$time
  }
  solved_model2 <- matrix(NA, nrow=nrow(solved_model1),ncol=n_samp)
  
  f <- create_multi_kinetics_posterior_function(parTab_use, solved_model1,NULL,model_func,ver = "model",use_log=use_log)
  
  for(x in seq_along(samps)){
    tmp_pars <- get_index_pars(chain, samps[x])
    solved_model2[,x] <- f(tmp_pars)
  }
  colnames(solved_model2) <- seq_along(samps)
  min_titer <- min(solved_model2, na.rm=TRUE)
  
  res <- bind_cols(solved_model1, solved_model2)
  res <- res %>% pivot_longer(-c(i,time,titer)) %>% rename(obs=titer)
  res <- res %>% rename(samp=name,titer=value)
  res <- res %>% mutate(titer = if_else(is.nan(titer), min_titer, titer))
  
  res_summary <- res %>% group_by(i, time) %>%
    summarize(lower95=quantile(titer,0.025),upper95=quantile(titer,0.975),median=median(titer))
  return(list(res,res_summary))
}

#' @export
plot_model_fits <- function(model_func, n_samp, indivs=NULL, chain, dat, expand_times=TRUE, times=NULL, use_log=TRUE,plot_draws=TRUE){
  tmp <- get_kinetics_prediction_intervals(model_func,n_samp, indivs, chain,  dat,expand_times, times, use_log)
  
  if(plot_draws){
    p1 <- ggplot(tmp[[1]]) + 
      geom_line(aes(x=time,y=titer,group=samp),alpha=0.1) + 
      geom_point(data=use_dat%>% filter(i %in% indivs),aes(x=time,y=titer),col="blue") +
      facet_wrap(~i)
  } else {
    p1 <- ggplot(tmp[[2]]) + 
      geom_ribbon(aes(x=time,ymin=lower95,ymax=upper95),alpha=0.25) +
      geom_line(aes(x=time,y=median)) + 
      geom_point(data=dat %>% filter(i %in% indivs),aes(x=time,y=titer),col="blue") +
      facet_wrap(~i)
  }
  p1 <- p1 + ylab("Titer") + xlab("Time") + theme_bw()
  p1
}

#' @export
plot_model_fits_mean <- function(model_func, chain, parTab, nsamp=100,times=seq(0,365,by=1),use_log=TRUE, plot_draws=TRUE){
  chain_means <- extract_population_distributions(chain, parTab, nsamp)
  chain <- chain_means %>% select(sampno, name, mean) %>% pivot_wider(names_from=name,values_from=mean)
  
  preds <- matrix(0, nrow=nsamp, ncol=length(times))
  for(i in 1:nrow(chain)){
    pars <- as.numeric(chain[i,])
    names(pars) <- colnames(chain)
    pars <- exp(pars[names(pars) != "sampno"])
    preds[i,] <- model_func(times, pars)
  }
  if(use_log){
    preds <- log10(preds)
  }
  preds <- reshape2::melt(preds)
  colnames(preds) <- c("samp","time","titer")
  preds_summary <- preds %>% group_by(time) %>%
    summarize(lower95=quantile(titer,0.025),upper95=quantile(titer,0.975),median=median(titer))
  
  if(plot_draws){
    p1 <- ggplot(preds) + 
      geom_line(aes(x=time,y=titer,group=samp),alpha=0.1) 
  } else {
    p1 <- ggplot(preds_summary) + 
      geom_ribbon(aes(x=time,ymin=lower95,ymax=upper95),alpha=0.25) +
      geom_line(aes(x=time,y=median))
  }
  p1 <- p1 + ylab("Titer") + xlab("Time") + theme_bw()
  p1
}

#' @export
plot_posteriors <- function(chain, parTab, nsamp=100, ver="means",indivs=NULL,transform=FALSE){
  if(!(ver %in% c("means","var","indivs"))){
    message("Error - ver must be one of means, var and indivs")
    return(0)
  }
  unique_samps <- unique(chain$sampno)
  samps <- sample(unique_samps, nsamp)
  chain <- chain[chain$sampno %in% samps,]
  
  colnames(chain) <- c("sampno", parTab[parTab$i == 0, "names"],paste0(parTab[parTab$i != 0,"names"],".",parTab[parTab$i != 0,"i"]),"lnlike")
  if(ver=="means"){
    use_cols <- which((parTab$names %like% "mean" | parTab$names == "obs_sd") & parTab$fixed == 0)
    chain <- chain[,c(1, use_cols + 1)]
    chain_melted <- chain %>% pivot_longer(-sampno)
  } else if(ver=="var"){
    use_cols <- which((parTab$names %like% "var" | parTab$names == "obs_sd") & parTab$fixed == 0)
    chain <- chain[,c(1, use_cols + 1)]
    chain_melted <- chain %>% pivot_longer(-sampno)
   
  } else {
    cols_means <- which(parTab$names %like% "mean")
    cols_vars <- which(parTab$names %like% "var")
    use_cols_indiv <- which(parTab$i %in% indivs & parTab$fixed == 0)
    chain_means <- chain[,c(1, cols_means+1)]
    chain_var <- chain[,c(1, cols_vars+1)]
    chain_indiv <- chain[,c(1, use_cols_indiv+1)]

    chain_means <- chain_means %>% pivot_longer(-sampno) %>% 
      mutate(name=str_remove(name, "_mean")) %>%
      rename(mean=value)
    chain_var <- chain_var %>% pivot_longer(-sampno) %>% 
      mutate(name=str_remove(name, "_var")) %>%
      rename(var=value)
    chain_indiv <- chain_indiv%>% pivot_longer(-sampno) %>%
      mutate(var = str_split(name, "[.]",simplify=TRUE)[,1]) %>%
      mutate(i = str_split(name, "[.]",simplify=TRUE)[,2]) %>%
      select(-name)%>% rename(name=var)
    chain_melted <- left_join(chain_indiv, chain_means) %>% left_join(chain_var)
    chain_melted <- chain_melted %>% mutate(value = exp(mean + var*value))
  }
  if(transform){
    chain_melted$value <- exp(chain_melted$value)
  }
  
  if(ver == "indivs"){
    p_trace <- ggplot(chain_melted) + 
      geom_line(aes(x=sampno, y=value)) + 
      facet_wrap(i~name,scales="free") +
      xlab("Sample") + ylab("Value")+
      theme_bw()
    p_density <- ggplot(chain_melted) + 
      geom_density(aes(x=value),fill="blue",alpha=0.25) + 
      facet_wrap(i~name,scales="free") +
      xlab("Value") + ylab("Density") +
      theme_bw()
  } else {
    p_trace <- ggplot(chain_melted) + geom_line(aes(x=sampno, y=value)) + facet_wrap(~name,scales="free") +
      xlab("Sample") + ylab("Value")+
      theme_bw()
    p_density <- ggplot(chain_melted) + geom_density(aes(x=value),fill="blue",alpha=0.25) + facet_wrap(~name,scales="free")+
      xlab("Value") + ylab("Density") +
      theme_bw()
  }

  return(list(p_trace, p_density))
}


#' @export
extract_population_distributions <- function(chain, parTab, nsamp){
  unique_samps <- unique(chain$sampno)
  samps <- sample(unique_samps, nsamp)
  chain <- chain[chain$sampno %in% samps,]
  
  colnames(chain) <- c("sampno", parTab[parTab$i == 0, "names"],paste0(parTab[parTab$i != 0,"names"],".",parTab[parTab$i != 0,"i"]),"lnlike")

    cols_means <- which(parTab$names %like% "mean")
    cols_vars <- which(parTab$names %like% "var")
   
    chain_means <- chain[,c(1, cols_means+1)]
    chain_var <- chain[,c(1, cols_vars+1)]
    chain_means <- chain_means %>% pivot_longer(-sampno) %>% 
      mutate(name=str_remove(name, "_mean")) %>%
      rename(mean=value)
    chain_var <- chain_var %>% pivot_longer(-sampno) %>% 
      mutate(name=str_remove(name, "_var")) %>%
      rename(var=value)
    chain_melted <- left_join(chain_means,chain_var) 
  chain_melted
}
