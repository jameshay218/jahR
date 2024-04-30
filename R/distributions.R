#' Distribution mean and variance to par1 and par2 (and vice versa)
#'
#' Given parameter 1 and 2 (e.g., mean and variance or shape and scale), finds parameter 1 and parameter 2 (or mean and variance) of the specified distribution. Supported distributions are: gamma, beta, beta-binomial, gumbel, geometric
#' @param dist character name of the distribution. one of "gamma","beta","beta-binomial","gumbel", "geometric"
#' @param par1 parameter 1 of the distribution (usually shape or alpha)
#' @param par2 parameter 2 of the distribution (usually scale or beta)
#' @param par3 OPTIONAL, but may be the 3rd distribution parameter e.g., N of the beta-binomial
#' @param mean_par mean of the distribution
#' @param var_par variance of the distribution
#' @return a named vector of the converted parameters
#' @export
convert_dist_pars <- function(dist="gamma",par1=NULL,par2=NULL,par3=NULL, mean_par=NULL,var_par=NULL){
    ## Figure out which direction
    get_mean <- !is.null(par1) & !is.null(par2) & is.null(mean_par) & is.null(var_par)
    get_pars <- is.null(par1) & is.null(par2) & !is.null(mean_par) & !is.null(var_par)

    if(!get_mean & !get_pars){
        message("Error - you have not passed enough parameters to convert either way")
      return(NULL)
    }

    ## Choose appropriate distribution and convert parameters
     if(dist=="gamma"){
       solve_fr <- rgamma
       solve_fd <- dgamma
        if(get_pars){
            ret_pars <- gamma_mean_var_to_shape_scale(mean_par,var_par)
        } else {
            ret_pars <- gamma_shape_scale_to_mean_var(par1,par2)
        }
     } else if (dist == "gumbel"){
       solve_fr <- evd::rgumbel
       solve_fd <- evd::dgumbel
         if(get_pars){
             ret_pars <- gumbel_mean_var_to_mu_beta(mean_par,var_par)
         } else {
             ret_pars <- gumbel_mu_beta_to_mean_var(par1,par2)
         }
     } else if(dist == "beta"){
       solve_fr <- rbeta
       solve_fd <- dbeta
         if(get_pars){
             ret_pars <- beta_mean_var_to_alpha_beta(mean_par,var_par)
         } else {
             ret_pars <- beta_alpha_beta_to_mean_var(par1,par2)
         }
     } else if(dist=="geometric"){
       solve_fr <- rgeom
       solve_fd <- dgeom
         if(get_pars){
             ret_pars <- geometric_mean_to_p(mean_par)
         } else {
             ret_pars <- geometric_p_to_mean_var(par1)
         }
     } else if(dist=="beta-binomial"){
       solve_fr <- rbb
       solve_fd <- dbb
         if(is.null(par3)){
             message("Error - par3 should be the N of the beta-binomial distribution")
             return(NULL)
         }
        if(get_pars){
            ret_pars <- bbinomial_mean_var_n_to_alpha_beta(mean_par, var_par, par3)
        } else {
            ret_pars <- bbinomial_alpha_beta_n_to_mean_var(par1, par2, par3)
        }
     } else {
       solve_fr <- rnorm
       solve_dr <- dnorm
        if(get_pars){
            ret_pars <- c(par1, par2)
        } else {
            ret_pars <- c(mean_par, var_par)
        }
     }

    ## Fill in NULL parameters for plot
    if(get_pars){
      par1 <- ret_pars[1]
      par2 <- ret_pars[2]
    } else {
      mean_par <- ret_pars[1]
      var_par <- ret_pars[2]
    }

    ## Generate plot data
    if(dist != "beta-binomial"){
      if(dist == "gamma"){
        plot_dat <- data.frame(x=solve_fr(10000,par1, scale=par2))
        x_range <- seq(min(plot_dat$x), max(plot_dat$x), by=(max(plot_dat$x)-min(plot_dat$x))/100)
        binwidth1 <- (max(plot_dat$x)-min(plot_dat$x))/25
        plot_dat2 <- data.frame(x=x_range,y=solve_fd(x_range, par1, scale=par2))
      } else {
        plot_dat <- data.frame(x=solve_fr(10000,par1, par2))
        x_range <- seq(min(plot_dat$x), max(plot_dat$x), by=(max(plot_dat$x)-min(plot_dat$x))/100)
        binwidth1 <- (max(plot_dat$x)-min(plot_dat$x))/25
        plot_dat2 <- data.frame(x=x_range,y=solve_fd(x_range, par1, par2))
      }
    } else {
      plot_dat <- data.frame(x=solve_fr(10000,par1, par2, par3))
      x_range <- floor(seq(min(plot_dat$x), max(plot_dat$x), by=(max(plot_dat$x)-min(plot_dat$x))/100))
      binwidth1 <- floor((max(plot_dat$x)-min(plot_dat$x))/25)
      plot_dat2 <- data.frame(x=x_range,y=solve_fd(x_range, par1, par2, par3))
    }

    ## Create histogram and density plot
    p <- ggplot(data=plot_dat) +
      theme_overall() +
      xlab("x") + ylab("Count") +
      ggtitle(paste0("10,000 random draws from ", dist, " distribution with mean=",
                     signif(mean_par,3),", var=",signif(var_par,3),
                     ";\nParameter 1=", signif(par1,3), " and parameter 2=", signif(par2,3)))

    if(dist %in% c("gamma","beta","gamma", "geometric")){
      p <- p +  geom_histogram(aes(x=x),fill="blue",col="black",alpha=0.25,boundary=0,bins=20)
    } else {
      p <- p +  geom_histogram(aes(x=x),fill="blue",col="black",alpha=0.25,bins=20)
    }
    if(dist != "beta-binomial"){
      p2 <- ggplot(data=plot_dat2) +
        geom_ribbon(aes(x=x,ymax=y,ymin=0),fill="blue",col="black",alpha=0.25) +
        theme_overall() +
        xlab("x") + ylab("Density") +
        ggtitle(paste0("PDF of ", dist, " distribution with mean=",signif(mean_par,3),", var=",signif(var_par,3), ";\nParameter 1=", signif(par1,3), " and parameter 2=", signif(par2,3)))
      p_main <- p/p2
    } else {
      p_main <- p
    }

    list(plot=p_main,pars=ret_pars)
}

#' Fit distribution to data vector
#'
#' Given parameter 1 and 2 (e.g., mean and variance or shape and scale), finds parameter 1 and parameter 2 (or mean and variance) of the specified distribution. A key use case is to fit to posterior draws to create priors for downstream MCMC. This uses optim (default Nelder-Mead). Supported distributions are: normal, gamma, beta, log-normal, geometric, poisson
#' @param dist character name of the distribution
#' @param x vector of data to be fitted to
#' @param ver either "draws" or "density". If "draws", estimates the parameters directly. If "density", fits the distribution to the density of x.
#' @param lower_bound OPTIONAL lower bound of optimization (uses Brent method)
#' @param upper_bound OPTIONAL upper bound of optimization (uses Brent method)
#' @param control list of controls for optim
#' @return a vector of the estimated parameters and a ggplot2 object showing the fit
#' @export
fit_distribution_quick <- function(dist="normal",x,ver="draws",lower_bound=NULL,upper_bound=NULL,control=list(abstol = 1e-8, reltol = 1e-8)){
    if(is.null(lower_bound) != is.null(upper_bound)){
        message("Warning - if bounds are required, then both upper and lower bounds must be specified")
    }
    use_method <- "Nelder-Mead"
    if(!is.null(lower_bound) & !is.null(lower_bound)){
        use_method<-"Brent"
    } else {
        lower_bound <- -Inf
        upper_bound <- Inf
    }

    ## Remove NA from the data to be fitted
    x <- x[!is.na(x)]

    if(ver=="draws"){
        dat_use <- x
        if(dist=="normal"){
            use_f <- fit_normal
            solve_f <- dnorm
            start_par <- c(mean(x,na.rm=TRUE),sd(x,na.rm=TRUE))
        } else if(dist=="gamma"){
            use_f <- fit_gamma
            solve_f <- dgamma
            start_par <- convert_dist_pars("gamma",mean_par=mean(x,na.rm=TRUE),var_par=var(x,na.rm=TRUE))
        } else if(dist=="log-normal"){
            use_f <- fit_log_normal
            solve_f <- dlnorm
            start_par <- c(mean(log(x),na.rm=TRUE),sd(log(x),na.rm=TRUE))
        } else if(dist == "poisson"){
            use_f <- fit_poisson
            solve_f <- dpois
            start_par <- c(mean(x,na.rm=TRUE))
        } else if(dist=="geometric"){
            use_f <- fit_geometric
            solve_f <- dgeom
            start_par <- c(0.25)
        } else if(dist=="beta"){
            use_f <- fit_beta
            solve_f <- dbeta
            start_par <- convert_dist_pars("beta",mean_par=mean(x,na.rm=TRUE),var_par=var(x,na.rm=TRUE))
        } else {
            use_f <- fit_normal
            solve_f <- dnorm
            start_par <- c(mean(x,na.rm=TRUE),sd(x,na.rm=TRUE))
        }
    } else {
        dat_use <- density(x)
        if(dist=="normal"){
            use_f <- fit_normal_density
            solve_f <- dnorm
            start_par <- c(mean(x,na.rm=TRUE),sd(x,na.rm=TRUE))
        } else if(dist=="gamma"){
            use_f <- fit_gamma_density
            solve_f <- dgamma
            start_par <- convert_dist_pars("gamma",mean_par=mean(x,na.rm=TRUE),var_par=var(x,na.rm=TRUE))
        } else if(dist=="log-normal"){
            use_f <- fit_log_normal_density
            solve_f <- dlnorm
            start_par <- c(mean(log(x),na.rm=TRUE),sd(log(x),na.rm=TRUE))
        } else if(dist=="beta"){
            use_f <- fit_beta_density
            solve_f <- dbeta
            start_par <- convert_dist_pars("beta",mean_par=mean(x,na.rm=TRUE),var_par=var(x,na.rm=TRUE))
        } else {
            use_f <- fit_normal_density
            solve_f <- dnorm
            start_par <- c(mean(x,na.rm=TRUE),sd(x,na.rm=TRUE))
        }
    }

    message("Starting parameters used: ")
    print(signif(start_par,3))
    ## Find the parameters of the distribution which best fit the vector
    res <- optim(par=start_par,fn=use_f, dat=dat_use,lower=lower_bound,upper=upper_bound,method=use_method)

    best_pars <- res$par
    min_x <- min(x,na.rm=TRUE)
    max_x <- max(x,na.rm=TRUE)
    by_x <- (max_x-min_x)/100
    solve_x <- seq(min_x,max_x,by=by_x)

    if(length(best_pars) > 1){
        plot_data <- data.frame(x=solve_x,y=solve_f(solve_x, best_pars[1],best_pars[2]))
    } else {
        plot_data <- data.frame(x=solve_x,y=solve_f(solve_x, best_pars[1]))
    }
    plot_data$ynorm <- plot_data$y/max(plot_data$y)

    best_pars1 <- best_pars
    if(length(best_pars1) == 1){
        best_pars[2] <- NA
    }

    p1 <- ggplot(plot_data) +
        geom_histogram(data=data.frame(x=x),aes(x=x,y=..ncount..,col="Data"),fill="grey70") +
        geom_line(aes(x=x,y=ynorm,col="Fitted"),size=0.75) +
        scale_color_manual(name="",values=c("Fitted"="blue","Data"="black")) +
        ylab("Normalized density") +
        xlab("x") +
        ggtitle(paste0("Fitted ",dist," distribution to provided data (raw)\n Parameter estimates: ",
                       signif(best_pars[1],3),", ",signif(best_pars[2],3))) +
        theme_overall()

    p2 <- ggplot(plot_data) +
        geom_density(data=data.frame(x=x),aes(x=x,col="Data"),fill="grey40",alpha=0.1) +
        geom_line(aes(x=x,y=y,col="Fitted")) +
        scale_color_manual(name="",values=c("Fitted"="blue","Data"="black")) +
        ylab("Density") +
        xlab("x") +
        ggtitle(paste0("Fitted ",dist," distribution to provided data (density)\n Parameter estimates: ",
                       signif(best_pars[1],3),", ",signif(best_pars[2],3))) +
        theme_overall()

    ## Plot the fit against the data
    return(list(plot=p1/p2,pars=best_pars))
}

## Helpers for distribution reparameterizations
gamma_mean_var_to_shape_scale <- function(mean_par, var_par){
    scale <- var_par/mean_par
    shape <- mean_par/scale
    return(c(shape=shape,scale=scale))
}
gamma_shape_scale_to_mean_var <- function(shape, scale){
    mean_par <- shape*scale
    var_par <- shape*(scale^2)
    return(c(mean=mean_par,var=var_par))
}

beta_mean_var_to_alpha_beta <- function(mean_par, var_par){
    beta_alpha <- mean_par * ((mean_par * (1-mean_par) / var_par) - 1)
    beta_beta <- beta_alpha*(1/mean_par - 1)
    return(c(shape1=beta_alpha,shape2=beta_beta))
}

beta_alpha_beta_to_mean_var <- function(alpha, beta){
    mean_par <- alpha/(alpha+beta)
    var_par <- alpha*beta / (((alpha+beta)^2)*(alpha+beta+1))
    return(c(mean=mean_par,var=var_par))
}

bbinomial_alpha_beta_n_to_mean_var <- function(alpha, beta, n){
    mean_par <- bb_mean(n,alpha,beta)
    var_par <- bb_var(n,alpha,beta)
    return(c(mean=mean_par, var=var_par))
}

bbinomial_mean_var_n_to_alpha_beta <- function(mean, var, n) {
    y <- mean
    z <- var
    x <- (n - y) / y
    top <- z * (1 + x)^2 - (n^2) * x
    bot <- n * x * (1 + x) - z * ((1 + x)^3)
    a2 <- top / bot
    b2 <- x * a2
    return(c("a" = a2, "b" = b2))
}

gumbel_mu_beta_to_mean_var <- function(mu_par, beta_par){
    mean_par <- mu_par + beta_par*exp(1)
    var_par <- pi^2 / 6 * beta_par^2
    return(c(mean=mean_par,var=var_par))
}
gumbel_mean_var_to_mu_beta <- function(mean_par, var_par){
    beta_par <- sqrt(var_par*6/pi^2)
    mu_par <- mean_par - beta_par*exp(1)
    return(c(location=mu_par,scale=beta_par))
}

geometric_mean_to_p <- function(mean_par){
    return(1/mean_par)
}

geometric_p_to_mean_var <- function(mean_par){
    return(c(mean=1/p,var=(1-p)/p^2))
}

## Cost functions to pass to optim
#' Function for optim to fit poisson distribution
fit_poisson <- function(pars, dat){
    mean <- pars[1]
    -sum(dpois(dat, mean, log=TRUE))
}

#' Function for optim to fit a gamma distribution
fit_gamma <- function(pars, dat){
    -sum(dgamma(dat, shape=pars[1], scale=pars[2], log=TRUE))
}
fit_gamma_density <- function(pars, dat) {
    shape1 <- pars[1]
    scale <- pars[2]
    out <- dgamma(dat$x, shape=shape1, scale=scale)
    return(sum((out - dat$y)^2))
}

#' Function for optim to fit a geometric distribution
fit_geometric <- function(prob, dat){
    -sum(dgeom(x=dat, prob, log=TRUE))
}

fit_beta <- function(pars, dat){
    -sum(dbeta(dat, pars[1], pars[2], log=TRUE))
}

fit_beta_density <- function(pars, dat){
    shape1 <- pars[1]
    shape2 <- pars[2]
    out <- dbeta(dat$x, shape1, shape2)
    return(sum((out - dat$y)^2))
}

#' Discretised gamma fit
dgamma_discrete_mean <- function(dat, mean, var, use_log=TRUE){
    scale <- var/mean
    shape <- mean/scale
    ddgamma(dat, shape=shape,scale=scale,log=use_log)
}

fit_gamma_discrete <- function(pars,dat){
    scale <- pars[1]
    shape <- pars[2]
    -sum(ddgamma(dat, shape=shape, scale=scale, log=TRUE))
}

fit_normal <- function(pars, dat){
    -sum(dnorm(dat, pars[1],pars[2],log=TRUE))
}

fit_normal_density <- function(pars, dat){
    shape1 <- pars[1]
    shape2 <- pars[2]
    out <- dnorm(dat$x, mean = shape1, sd = shape2)
    return(sum((out - dat$y)^2))
}

fit_log_normal <- function(pars, dat){
    -sum(dlnorm(dat, pars[1],pars[2],log=TRUE))
}

fit_log_normal_density <- function(pars, dat){
    shape1 <- pars[1]
    shape2 <- pars[2]
    out <- dlnorm(dat$x, mean = shape1, sd = shape2)
    return(sum((out - dat$y)^2))
}

#' Estimate vector mode
#'
#' @param x the vector to be estimated
#' @return the estimated mode of the given vector of values
#' @examples
#' x <- runif(1000)
#' y <- estimate_mode(x)
#' @export
estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

#' Creates a covariance matrix from a vector of standard deviations and a correlation matrix
#' 
#' @param sds vector of standard deviations
#' @param cor_mat a correlation matrix with nrow and ncol the same as length(sds)
#' @return the covariance matrix
#' @export
create_cov_mat <- function(sds, cor_mat){
  diag(sds) %*% cor_mat %*% diag(sds)
}