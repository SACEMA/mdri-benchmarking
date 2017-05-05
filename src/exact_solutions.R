library(cubature)
library(mnormt)
library(pracma)
library(R2Cuba)
library(reshape2)
library(optparse)
#library(RMySQL) # <- This can be a massive pain to setup in windows and its only needed for the database part

# This script computes MDRI for a biomarker signal in the form given by the inner6 function.
# To change the biomarker, change that function
# There are two integration steps:
# An inner integral over T
# An outer integral over beta, gamma and delta
# Multiple techniques are implemented for each integral.
# Using integrate and adaptIntegrate for the inner and outer integrals respectively has worked pretty well for me - Note that some of the methods are very slow, expecially when applied to the outer integral
# it was originally written to work with a database. The params variable was populated from the database, you have to do it in some other way. I hard coded a couple of values into it so that you can see how it works.
# It was also originally written to be called from the commandline so you will see some weird stuff at the start and end of the script that you can just delete if you dont want to use it from the commandline.
# full_run is the function that actually call all the other functions and computes the integral
# Sorry about the line length - best viewed on a full hd screen with small font size :)

# Command line arguments
option_list <- list(
    make_option(c("-t", "--threshold"), action = "store", help = "Threshold to use", dest = 'threshold', default = 40, type = "integer"),
    make_option(c("-T", "--bigT"), action = "store", help = "bigT", dest = 'bigT', default = 365, type = "integer"),
    make_option(c("-i", "--inner-method"), action = "store", help = "Method for the inner part of the integral (1d). integrate, quad, quadl, integral_kronrod, integral_richardson, integral_clenshaw, integral_simpson, integral_romberg, adaptIntegrate, suave, vegas, cuhre", default = "integrate", dest = 'innerMethod'),
    make_option(c("-o", "--outer-method"), action = "store", help = "Method for the outer part of the integral (3d). adaptIntegrate, suave, vegas", dest = 'outerMethod', default = 'adaptIntegrate'),
    make_option(c("-c", "--calculate"), action = "store_true", help = "Specify this option if you want to calculate an integral", dest = "calculate")
    )

opt <- parse_args(OptionParser(option_list = option_list))
threshold <- opt$threshold
bigT <- opt$bigT
innerMethod <- opt$innerMethod
outerMethod <- opt$outerMethod

# This is where I read the parameters from the database. You have to setup the params variable some other way

# con <- dbConnect(MySQL(), user = 'root', password = "xxxxxxxx", host = "localhost", db = "assay_calib_sims")
# parsFromDB <- dbGetQuery(con, paste("SELECT 'a', variable, value FROM biology_parameter_sets where ps = '", par_set, "' AND biol_id = ", biol_id, sep = ""))
# parsFromDB[,3]<-as.numeric(parsFromDB[,3])
# params <- dcast(parsFromDB, a ~ variable, mean)

# Temporary stop gap to illustrate how the parameters should look
# NOTE CAREFULLY: this is a data.frame
params <- structure(list(a = "a", alpha_beta_sd = 0, alpha_delta_sd = 0, 
    alpha_gamma_sd = 0, alpha_mu = 0, alpha_sd = 0, beta_delta_sd = 3.15, 
    beta_gamma_sd = -28, beta_mu = 4, beta_sd = 1.4, delta_mu = 85, 
    delta_sd = 7.5, e0 = 2, e1 = 0, e2 = 0.3, e3 = 0.5, gamma_delta_sd = -45, 
    gamma_mu = 190, gamma_sd = 50), .Names = c("a", "alpha_beta_sd", 
"alpha_delta_sd", "alpha_gamma_sd", "alpha_mu", "alpha_sd", "beta_delta_sd", 
"beta_gamma_sd", "beta_mu", "beta_sd", "delta_mu", "delta_sd", 
"e0", "e1", "e2", "e3", "gamma_delta_sd", "gamma_mu", "gamma_sd"
), row.names = c(NA, -1L), class = "data.frame")

stopifnot(outerMethod %in% c("adaptIntegrate", "suave", "vegas"))
stopifnot(innerMethod %in% c("integrate", "quad", "quadl", "integral_kronrod", "integral_richardson", "integral_clenshaw", "integral_simpson", "integral_romberg", "adaptIntegrate", "suave", "vegas", "cuhre"))

inner6 <- function(t, beta, gamma, delta, e0, e1, e2, e3, threshold)
{
#    print(c(t, beta, gamma, delta, e0, e1, e2, e3, threshold)) # Useful for debugging
    alpha <- 0
    seroconversion_date <- 0
    Z <- ((alpha - delta)/(1+(((t-seroconversion_date)/gamma)**beta)))+delta
    error_term <- e2*(Z**e3) + e1*(Z) + e0
    pnorm((threshold - Z)/error_term)
}

inner_int6 <- function(x, e0, e1, e2, e3, threshold, bigT, innerMethod)
{
    beta = x[1]
    gamma = x[2]
    delta = x[3]
    if (innerMethod == "integrate")
        {answ <- try(integrate(inner6, 0, bigT, beta = beta, gamma = gamma, delta = delta, e0 = e0, e1 = e1, e2 = e2, e3 = e3, threshold = threshold)$value)}
    if (innerMethod == "quad")
        {answ <- try(quad(inner6, 0, bigT, beta = beta, gamma = gamma, delta = delta, e0 = e0, e1 = e1, e2 = e2, e3 = e3, threshold = threshold))}
    if (innerMethod == "quadl")
        {answ <- try(quadl(inner6, 0, bigT, beta = beta, gamma = gamma, delta = delta, e0 = e0, e1 = e1, e2 = e2, e3 = e3, threshold = threshold))}
    if (innerMethod == "integral_kronrod")
        {answ <- try(integral(inner6, 0, bigT, method = "Kronrod", beta = beta, gamma = gamma, delta = delta, e0 = e0, e1 = e1, e2 = e2, e3 = e3, threshold = threshold))}
    if (innerMethod == "integral_richardson")
        {answ <- try(integral(inner6, 0, bigT, method = "Richardson", beta = beta, gamma = gamma, delta = delta, e0 = e0, e1 = e1, e2 = e2, e3 = e3, threshold = threshold))}
    if (innerMethod == "integral_clenshaw")
        {answ <- try(integral(inner6, 0, bigT, method = "Clenshaw", beta = beta, gamma = gamma, delta = delta, e0 = e0, e1 = e1, e2 = e2, e3 = e3, threshold = threshold))}
    if (innerMethod == "integral_simpson")
        {answ <- try(integral(inner6, 0, bigT, method = "Simpson", beta = beta, gamma = gamma, delta = delta, e0 = e0, e1 = e1, e2 = e2, e3 = e3, threshold = threshold))}
    if (innerMethod == "integral_romberg")
        {answ <- try(integral(inner6, 0, bigT, method = "Romberg", beta = beta, gamma = gamma, delta = delta, e0 = e0, e1 = e1, e2 = e2, e3 = e3, threshold = threshold))}
    if (innerMethod == "adaptIntegrate")
        {answ <- try(adaptIntegrate(inner6, lowerLimit = 0, upperLimit = bigT, beta = beta, gamma = gamma, delta = delta, e0 = e0, e1 = e1, e2 = e2, e3 = e3, threshold = threshold)$integral)}
    if (innerMethod == "suave")
        {answ <- try(suave(ndim = 1, ncomp = 1, integrand = inner6, beta = beta, gamma = gamma, delta = delta, 
                           e0 = e0, e1 = e1, e2 = e2, e3 = e3, threshold = threshold, lower = 0, upper = bigT,
                           flags=list(verbose=0, final=1, pseudo.random=0, smooth=0, mersenne.seed=NULL))$value)}
    if (innerMethod == "vegas")
        {answ <- try(vegas(ndim = 1, ncomp = 1, integrand = inner6, beta = beta, gamma = gamma, delta = delta, 
                           e0 = e0, e1 = e1, e2 = e2, e3 = e3, threshold = threshold, lower = 0, upper = bigT,
                           flags=list(verbose=0, final=1, pseudo.random=0, smooth=0, mersenne.seed=NULL))$value)}
    if(class(answ) == "try-error")
    {
        answ <- romberg(inner6, 0, bigT, beta = beta, gamma = gamma, delta = delta, e0 = e0, e1 = e1, e2 = e2, e3 = e3, threshold = threshold)
        cat("non-convergence")
    }
    answ
}

inner_int_dens6only <- function(x, e0, e1, e2, e3, threshold, bigT, mean_vec, cov_mat, innerMethod)
{
    y <- dmnorm(x, mean_vec, cov_mat)
    y
}

inner_int_dens6 <- function(x, e0, e1, e2, e3, threshold, bigT, mean_vec, cov_mat, innerMethod)
{
    z <- inner_int6(x, e0, e1, e2, e3, threshold, bigT, innerMethod)
    y <- dmnorm(x, mean_vec, cov_mat)
    z*y
}

outer6 <- function(params, threshold, bigT, innerMethod, outerMethod)
{
    # really should have done this with 'with' instead
    e0 <- params$e0
    e1 <- params$e1
    e2 <- params$e2
    e3 <- params$e3
    beta_mu <- params$beta_mu
    gamma_mu <- params$gamma_mu
    delta_mu <- params$delta_mu
    beta_sd <- params$beta_sd
    gamma_sd <- params$gamma_sd
    delta_sd <- params$delta_sd
    beta_gamma_sd <- params$beta_gamma_sd
    beta_delta_sd <- params$beta_delta_sd
    gamma_delta_sd <- params$gamma_delta_sd
    
    mean_vec <- c(beta_mu, gamma_mu, delta_mu)
    cov_mat <- matrix(c(beta_sd**2, beta_gamma_sd, beta_delta_sd, 
                      beta_gamma_sd, gamma_sd**2, gamma_delta_sd,
                      beta_delta_sd, gamma_delta_sd, delta_sd**2), 
                      nrow = 3, ncol = 3)

    if (outerMethod == "adaptIntegrate")
    {
    answ <- adaptIntegrate(inner_int_dens6, 
                   lowerLimit = c(max(0, beta_mu - 4*beta_sd), max(0, gamma_mu - 4*gamma_sd), max(0, delta_mu - 4*delta_sd)),
                   upperLimit = c(beta_mu + 4*beta_sd, gamma_mu + 4*gamma_sd, delta_mu + 4*delta_sd),
                   e0 = e0, e1 = e1, e2 = e2, e3 = e3,
                   threshold = threshold, bigT=bigT,
                   innerMethod = innerMethod,
                   mean_vec = mean_vec, cov_mat = cov_mat, maxEval = 1000000
                   )
    normFactor <- adaptIntegrate(inner_int_dens6only, 
                   lowerLimit = c(max(0, beta_mu - 4*beta_sd), max(0, gamma_mu - 4*gamma_sd), max(0, delta_mu - 4*delta_sd)),
                   upperLimit = c(beta_mu + 4*beta_sd, gamma_mu + 4*gamma_sd, delta_mu + 4*delta_sd),
                   e0 = e0, e1 = e1, e2 = e2, e3 = e3,
                   threshold = threshold, bigT=bigT,
                   innerMethod = innerMethod,
                   mean_vec = mean_vec, cov_mat = cov_mat, maxEval = 1000000
                   )
    error <- answ$error/normFactor$integral
    answ <- answ$integral/normFactor$integral
    error <- error*answ
    }
    if (outerMethod == "suave")
    {
    answ <- suave(ndim = 3, ncom = 1, inner_int_dens6,  rel.tol = 0.001,
                   e0 = e0, e1 = e1, e2 = e2, e3 = e3,
                   threshold = threshold, bigT=bigT,
                   innerMethod = innerMethod,
                   mean_vec = mean_vec, cov_mat = cov_mat,
                   lower = c(max(0, beta_mu - 4*beta_sd), max(0, gamma_mu - 4*gamma_sd), max(0, delta_mu - 4*delta_sd)),
                   upper = c(beta_mu + 4*beta_sd, gamma_mu + 4*gamma_sd, delta_mu + 4*delta_sd), max.eval = 1000000
                   )
    normFactor <- suave(ndim = 3, ncom = 1, inner_int_dens6only,  rel.tol = 0.001,
                   e0 = e0, e1 = e1, e2 = e2, e3 = e3,
                   threshold = threshold, bigT=bigT,
                   innerMethod = innerMethod,
                   mean_vec = mean_vec, cov_mat = cov_mat,
                   lower = c(max(0, beta_mu - 4*beta_sd), max(0, gamma_mu - 4*gamma_sd), max(0, delta_mu - 4*delta_sd)),
                   upper = c(beta_mu + 4*beta_sd, gamma_mu + 4*gamma_sd, delta_mu + 4*delta_sd), max.eval = 1000000
                   )
    error <- answ$abs.error/normFactor$value
    answ <- answ$value/normFactor$value
    }
    if (outerMethod == "vegas")
    {
    answ <- vegas(ndim = 3, ncom = 1, inner_int_dens6,  rel.tol = 0.001,
                   e0 = e0, e1 = e1, e2 = e2, e3 = e3,
                   threshold = threshold, bigT=bigT,
                   innerMethod = innerMethod,
                   mean_vec = mean_vec, cov_mat = cov_mat,
                   lower = c(max(0, beta_mu - 4*beta_sd), max(0, gamma_mu - 4*gamma_sd), max(0, delta_mu - 4*delta_sd)),
                   upper = c(beta_mu + 4*beta_sd, gamma_mu + 4*gamma_sd, delta_mu + 4*delta_sd), max.eval = 1000000
                   )
    normFactor <- vegas(ndim = 3, ncom = 1, inner_int_dens6only,  rel.tol = 0.001,
                   e0 = e0, e1 = e1, e2 = e2, e3 = e3,
                   threshold = threshold, bigT=bigT,
                   innerMethod = innerMethod,
                   mean_vec = mean_vec, cov_mat = cov_mat,
                   lower = c(max(0, beta_mu - 4*beta_sd), max(0, gamma_mu - 4*gamma_sd), max(0, delta_mu - 4*delta_sd)),
                   upper = c(beta_mu + 4*beta_sd, gamma_mu + 4*gamma_sd, delta_mu + 4*delta_sd), max.eval = 1000000
                   )
    error <- answ$abs.error/normFactor$value
    answ <- answ$value/normFactor$value
    }
    return(list(answ = answ, err = error))
}

# The function that you should actually call.
full_run <- function(params, bigT, threshold, innerMethod, outerMethod)
{
    results <- NULL
    for (i in 1:nrow(params)) # If the params data.frame has multiple rows, then run multiple integrals
    {
        integralValue <- outer6(params[i,], threshold = threshold, bigT = bigT, innerMethod = innerMethod, outerMethod = outerMethod)
        results <- rbind(results, data.frame(i = i, innerMethod = innerMethod, outerMethod = outerMethod, 
                                             threshold = threshold, bigT = bigT,
                                             value = integralValue$answ, error = integralValue$err
                                             )
        )
    }
    # Stash the results in the database
    # dbWriteTable(con, "zExactSolutions", results, append = TRUE)
    return (results)
}

if (opt$calculate){
    print (opt)
    print (full_run(params, bigT, threshold, innerMethod, outerMethod))
}
