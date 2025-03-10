######################################################################
# draw all figures 
# 2025-03-09 
######################################################################

######################################################################
# Figure 1. precision rate vs log(lambda) 
######################################################################

rm(list = ls()); 
library(ggplot2);
figdf = read.csv(file.path("data", "df-precision-power.csv"));

for(icase in c("ind", "blockdiag", "ar1", "banded"))
{
    fig = ggplot(subset(figdf, case == icase), 
            aes(x = loglambda, y = precision.theory, col = as.factor(sparsity), linetype = as.factor(sparsity))) + 
                geom_line(linewidth = 1) + geom_point() + theme_bw() + 
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
                labs(x = expression(log(lambda)), y = "precision rate") + 
                labs(color = "sparsity", linetype = "sparsity") + 
                geom_errorbar(aes(ymin = precision.low, ymax = precision.upp), linetype = 1, width = 0.05); 

    ggsave(paste0("fig-1-precision-", icase, ".pdf"), plot = fig, 
        width = 6.4, height = 4.8);
}

######################################################################
# Figure 2. optimal precision rate vs alpha 
######################################################################

rm(list = ls()); 
library(ggplot2);
figdf = read.csv(file.path("data", "df-optimal-precision.csv"));

fig = ggplot(figdf, 
        aes(x = alpha, y = precision.theory, col = as.factor(sparsity), linetype = as.factor(sparsity))) + 
                geom_line(linewidth = 1) + geom_point() + theme_bw() + 
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
                labs(x = expression(alpha), y = "precision rate") + 
                labs(color = "sparsity", linetype = "sparsity") + 
                geom_errorbar(aes(ymin = CI.low, ymax = CI.upp), linetype = 1, width = 0.02); 

ggsave(paste0("fig-2-optimal-precision.pdf"), plot = fig, 
        width = 6.4, height = 4.8);

######################################################################
# Figure 3. histograms for w and debiased-w
######################################################################

rm(list = ls()); 
load(file.path("data", "data-histogram.Rdata"));

pdf("fig-3-w.pdf", paper = "special", width = 6.4, height = 4.8);
par(mar = c(4, 4, 0.5, 0.5));
hist(b.hat, breaks = 20, prob = TRUE, main = "", xlab = "estimated w"); 
dev.off(); 

index = which(abs(w0) > 1e-5); 
p1 = hist(w.bar[index]);
p2 = hist(w.bar[-index]);
p1$density = p1$density * sparsity; 
p2$density = p2$density * (1 - sparsity); 
xgrid = seq(-30, 30, length = 101);

pdf("fig-3-w-debiased.pdf", paper = "special", width = 6.4, height = 4.8);
par(mar = c(4, 4, 0.5, 0.5));
plot(p2, col = rgb(0, 0, 1, 1/4), xlim = c(-15, 20), freq = FALSE, main = "", xlab = "debiased w");
plot(p1, col = rgb(1, 0, 0, 1/4), xlim = c(-15, 20), freq = FALSE, add = TRUE);
lines(xgrid, (1 - sparsity) * dnorm(xgrid, sd = w.sd), col = "red");
lines(xgrid, sparsity * dnorm(xgrid, mean = mu.true, sd = w.sd), col = "blue");
dev.off(); 

######################################################################
# Figure 4. boxplots for confidence levels  
######################################################################

rm(list = ls()); 
library(ggplot2);
figdf = read.csv(file.path("data", "df-CI.csv"));

for(icase in c("ind", "blockdiag", "ar1", "banded"))
{

    fig = ggplot(subset(figdf, case == icase), 
            aes(x = factor(loglambda), y = CIprob2)) + 
               geom_boxplot() + geom_hline(yintercept = 0.95, col = "red") + 
               theme_bw() + xlab(expression(log(lambda))) + ylab("confidence level") +
               theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
               facet_grid(sparsity ~ .);

    ggsave(paste0("fig-4-CI-", icase, ".pdf"), plot = fig, 
        width = 6.4, height = 4.8);
}

######################################################################
# Figure 5. power vs log(lambda)
######################################################################

rm(list = ls()); 
library(ggplot2);
figdf = read.csv(file.path("data", "df-precision-power.csv"));

for(icase in c("ind", "blockdiag", "ar1", "banded"))
{
   fig = ggplot(subset(figdf, case == icase), 
            aes(x = loglambda, y = power.theory, col = as.factor(sparsity), linetype = as.factor(sparsity))) + 
                geom_line(linewidth = 1) + geom_point() + theme_bw() + 
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
                xlab(expression(log(lambda))) + ylab("power") + 
                labs(color = "sparsity", linetype = "sparsity") + 
                geom_errorbar(aes(ymin = power.low, ymax = power.upp), linetype = 1, width = 0.05); 

    ggsave(paste0("fig-5-power-", icase, ".pdf"), plot = fig, 
        width = 6.4, height = 4.8);
}

######################################################################
# Figure 6. ROC, case = "ind" and "ar1", log(lambda) = -3, -1.5, 0.5 
######################################################################

rm(list = ls()); 
library(ggplot2);
figdf = read.csv(file.path("data", "df-ROC.csv"));

for(icase in c("ind", "ar1"))
for(isparsity in c(0.01, 0.05, 0.1))
{
    fig = ggplot(subset(figdf, (case == icase) & (abs(sparsity - isparsity) < 1e-5)), 
            aes(x = false.pos, y = true.pos, col = as.factor(loglambda), linetype = as.factor(loglambda))) + 
            geom_line(linewidth = 1) + theme_bw() + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
            labs(x = "false positive rate", y = "true positive rate") + 
            labs(color = expression(log(lambda)), linetype = expression(log(lambda))); 

    ggsave(paste0("fig-6-ROC-", icase, "-", isparsity * 100, ".pdf"), plot = fig, 
        width = 6.4, height = 4.8);
}

######################################################################
# Figure 7. AUC vs log(lambda)
######################################################################

rm(list = ls()); 
library(ggplot2);
figdf = read.csv(file.path("data", "df-AUC-tau.csv"));

for(icase in c("ind", "blockdiag", "ar1", "banded"))
{
    fig = ggplot(subset(figdf, icase == case), 
            aes(x = loglambda, y = auc.mean, col = as.factor(sparsity), linetype = as.factor(sparsity))) + 
                geom_point() + theme_bw() + 
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
                xlab(expression(log(lambda))) + ylab("AUC") + 
                labs(color = "sparsity", linetype = "sparsity") + 
                geom_errorbar(aes(ymin = auc.low, ymax = auc.upp), linetype = 1, width = 0.05); 

    ggsave(paste0("fig-7-AUC-", icase, ".pdf"), plot = fig, 
        width = 6.4, height = 4.8);
}

######################################################################
# Figure 8. tau vs log(lambda)
######################################################################

rm(list = ls()); 
library(ggplot2);
figdf = read.csv(file.path("data", "df-AUC-tau.csv"));

for(icase in c("ind", "blockdiag", "ar1", "banded"))
{
    fig = ggplot(subset(figdf, case == icase), 
            aes(x = loglambda, y = tau, col = as.factor(sparsity), linetype = as.factor(sparsity))) + 
                geom_line() + geom_point() + theme_bw() + 
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
                xlab(expression(log(lambda))) + ylab(expression(tau)) + 
                labs(color = "sparsity", linetype = "sparsity"); 

    ggsave(paste0("fig-8-tau-", icase, ".pdf"), plot = fig, 
        width = 6.4, height = 4.8);
}

######################################################################
# Figure 9. precision rate at tau vs alpha
######################################################################

rm(list = ls()); 
library(ggplot2);
figdf = read.csv(file.path("data", "df-optimal-tau.csv"));

fig = ggplot(figdf, 
        aes(x = alpha, y = precision.theory, col = as.factor(sparsity), linetype = as.factor(sparsity))) + 
                geom_line(linewidth = 1) + geom_point() + theme_bw() + 
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
                labs(x = expression(alpha), y = "precision rate") + 
                labs(color = "sparsity", linetype = "sparsity") +  
                geom_errorbar(aes(ymin = CI.low, ymax = CI.upp), linetype = 1, width = 0.02); 

ggsave("fig-9-optimal-tau.pdf", plot = fig, 
        width = 6.4, height = 4.8);

######################################################################
# THE END
######################################################################
