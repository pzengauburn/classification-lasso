######################################################################
# monte carlo  
######################################################################

rm(list = ls()); 
library(glmnet); 

data_output_dir = file.path("..", "results-Rdata");

case = "ind";                         # which correlation structure
rho = 0.8;                            # rho in correlation matrix  
alpha = 0.5;                          # ratio of n and p 
sparsity = 0.05;                      # proportion of non-zeros 

iter = 500; 

loglambda.seq = seq(-4, 2, by = 0.5); 

######################################################################
# load model setting 
######################################################################

p = 1000; n = p * alpha; 

Rmat = diag(rep(1.0, p));

sigma2 = 2.0; 
Sigma = sigma2 * Rmat; 

mu0 = rep(0, p); 
mu0[seq(1, p, by = 1 / sparsity)] = 1; 
mu0 = Sigma %*% mu0; 
mu.dir = mu0 / sqrt(sum(mu0 * mu0)); 

mu.length = 2.0;
mu = mu.length * mu.dir; 

######################################################################
# compute Sigma.half and Sigma.jj 
######################################################################

eig = eigen(Sigma); 
Sigma.half = eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors); 

Sigma.inv = eig$vectors %*% diag(1/eig$values) %*% t(eig$vectors); 
Sigma.jj = diag(Sigma.inv); 

######################################################################
# simulate samples from multivariate normal distribution 
# mean = mu and Sigma^{1/2} = Sigma.half
######################################################################

rmvnorm = function(n, mu, Sigma.half)
{
    stopifnot(length(mu) == nrow(Sigma.half), length(mu) == ncol(Sigma.half));
    p = length(mu); 
    sweep(matrix(rnorm(n * p), nrow = n) %*% Sigma.half, 2, mu, "+"); 
}

NN = length(loglambda.seq); 
precision.rate.mc = array(0, dim = c(iter, NN));  

for(i in 1:iter)
{
    ######################################################################
    # simulate data and fit model 
    ######################################################################

    y = c(rep(1, n/2), rep(-1, n/2));  
    x = matrix(0, nrow = length(y), ncol = p); 
    x[y == 1, ] = rmvnorm(sum(y == 1),  mu, Sigma.half);   
    x[y == -1, ] = rmvnorm(sum(y == -1), -mu, Sigma.half); 
    fit = glmnet(x, y, family = "binomial", alpha = 1, standardize = FALSE, intercept = FALSE);
    b.hat = coef(fit, s = exp(loglambda.seq) / sqrt(p) / alpha, exact = TRUE, x = x, y = y)[-1, ];
    b.hat = b.hat * sqrt(p); 

    ######################################################################
    # precision rate 
    ######################################################################

    # new data set to calculate the predicted precision rate
    n0 = 500; 
    y0 = c(rep(1, n0/2), rep(-1, n0/2)); 
    x0 = matrix(0, nrow = length(y0), ncol = p); 
    x0[y0 == 1, ]  = rmvnorm(sum(y0 == 1),  mu, Sigma.half); 
    x0[y0 == -1, ] = rmvnorm(sum(y0 == -1), -mu, Sigma.half);
    precision.rate.mc[i, ] = colMeans(y0 * (x0 %*% b.hat) > 0);
}

######################################################################
# save results 
######################################################################

if(!dir.exists(data_output_dir)) dir.create(data_output_dir); 

file_output = paste0("mc_", case, "_", rho, "_", alpha, "_", sparsity, ".Rdata"); 

save(case, rho, alpha, sparsity, n, p, mu, Sigma, Sigma.jj, iter, 
     precision.rate.mc, loglambda.seq, 
     file = file.path(data_output_dir, file_output));

cat("Stop:", date(), "\n\n");

######################################################################
# THE END 
######################################################################
