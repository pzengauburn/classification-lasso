######################################################################
# compute saddle points  
######################################################################

rm(list = ls()); 
source("saddle-point.R"); 

output.filename = "saddle-list-iid.Rdata";

######################################################################
# set parameters 
######################################################################

case.seq = c("ind");
alpha.seq = seq(0.5, 1.5, by = 0.5); 
sparsity.seq = c(0.01, 0.05, 0.1); 
loglambda.seq = seq(-4, 2, by = 2); 
rho = 0.8; 

p = 1000; 

penalty.option = "phi-lasso"; 

######################################################################
# set saddle.list.iid
######################################################################

saddle.list.iid = vector("list", 
        length(case.seq) * length(alpha.seq) * length(sparsity.seq) * length(loglambda.seq));

count = 0; 
for(case in case.seq) 
for(alpha in alpha.seq)
for(sparsity in sparsity.seq)
for(loglambda in loglambda.seq)
{
    count = count + 1; 
    n = p * alpha; 

    saddle.list.iid[[count]]$case = case; 
    saddle.list.iid[[count]]$alpha = alpha; 
    saddle.list.iid[[count]]$sparsity = sparsity; 
    saddle.list.iid[[count]]$loglambda = loglambda; 

    saddle.list.iid[[count]]$rho = rho; 
    saddle.list.iid[[count]]$p = p; 
    saddle.list.iid[[count]]$n = n;  
    saddle.list.iid[[count]]$mu.length = 2.0;
    saddle.list.iid[[count]]$saddle = NULL; 
}

######################################################################
# compute parameters 
######################################################################

start.time = date(); 

for(i in 1:length(saddle.list.iid))
{
    with(saddle.list.iid[[i]],
        cat("\ncase =", case, "alpha =", alpha, 
            "sparsity =", sparsity, "loglambda =", loglambda, "\n"));

    ##################################################################
    # load model setting 
    ##################################################################

    Rmat = diag(rep(1.0, saddle.list.iid[[i]]$p));

    sigma2 = 2.0; 
    Sigma = sigma2 * Rmat; 

    mu0 = rep(0, saddle.list.iid[[i]]$p); 
    mu0[seq(1, saddle.list.iid[[i]]$p, by = 1 / saddle.list.iid[[i]]$sparsity)] = 1; 
    mu0 = Sigma %*% mu0; 
    mu.dir = mu0 / sqrt(sum(mu0 * mu0)); 

    mu = saddle.list.iid[[i]]$mu.length * mu.dir; 

    ##################################################################
    # compute saddle points 
    ##################################################################

    saddle.ini = NULL; 

    try({
        saddle = saddle.point(alpha = saddle.list.iid[[i]]$alpha, mu = mu, Sigma = Sigma, 
                lambda = exp(saddle.list.iid[[i]]$loglambda), ini = saddle.ini, 
                ncores = 2, loss = "logistic", penalty = penalty.option,
                maxiter = 1000, tol = 1e-4, mc.sample = 1e+5, ratio = 0.5, 
                trace.params = TRUE, details = FALSE)

        cat("iter =", saddle$iter, ", conv =", saddle$conv, "\n\n")
        saddle.list.iid[[i]]$saddle = saddle; 
    }); 
}

stop.time = date(); 

######################################################################
# compute saddle.list.info 
######################################################################

saddle.list.info = NULL; 
for(i in 1:length(saddle.list.iid))
if(!is.null(saddle.list.iid[[i]]$saddle))
{
    tmp1 = as.data.frame(saddle.list.iid[[i]][
            c("case", "rho", "n", "p", "alpha", "sparsity", "loglambda", "mu.length")]); 
    tmp2 = as.data.frame(saddle.list.iid[[i]]$saddle[
            c("iter", "conv", "q", "R", "q0", "zeta.0", "zeta", "R0")]); 
    saddle.list.info = rbind(saddle.list.info, cbind(tmp1, tmp2)); 
}

saddle.list.info$precision.theory = 
    pnorm(saddle.list.info$mu.length * saddle.list.info$R / sqrt(saddle.list.info$q0));

######################################################################
# save saddle.list.iid and saddle.list.info 
######################################################################

save(start.time, stop.time, saddle.list.iid, saddle.list.info, 
    file = output.filename);

######################################################################
# THE END 
######################################################################
