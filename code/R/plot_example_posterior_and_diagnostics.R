library(foreign)
library(ggplot2)
library(rstan)
library(plyr)
library(reshape2)
library(gridExtra)

load("ck_data.RData")

axis_line_color <- "#222222"
theme_set(theme_classic())
theme_update(
  axis.line.x = element_line(size = 3, color = axis_line_color),
  axis.line.y = element_blank(), 
  axis.ticks.y = element_blank(), 
  axis.text.y = element_blank(),
  axis.title = element_text(face = "bold", size = 13),
  strip.background = element_rect(fill = "gray35", color = "gray35"),
  strip.text = element_text(color = "white", face = "bold", size = 15),
  legend.position = "none"
)


# plot example posterior of bias for a single congress ----------------------------
bias <- extract(stanfit, pars = "bias")[[1]]
congress <- sort(unique(ck_data$congress))
t <- sample(length(congress), size = 1)

bias <- bias[,t]
mean_bias <- mean(bias)
sd_bias <- sd(bias)

graph <- ggplot(data.frame(x = bias), aes(x)) + 
  geom_density(fill = "#222222") + 
  stat_function(size = 1, color = "skyblue", fun = function(x) {
    dnorm(x, mean_bias, sd_bias)
    }) +
  labs(x = paste0("Bias: ", congress[t],"th congress"), y = "") 
graph
  



# plot diagnostics --------------------------------------------------------
Ndraws <- length(extract(stanfit,pars="lp__")[[1]])
stats <- rstan::summary(stanfit)$summary[,c("Rhat","n_eff","sd","se_mean")]
stats <- as.data.frame(stats)
stats$neff_ratio <- stats$n_eff/Ndraws
stats$mcse_ratio <- stats$se_mean / stats$sd
graph_diagostics <- ggplot(stats)

graph_rhat <- graph_diagostics + 
  stat_bin(aes(x = Rhat, y=..count../sum(..count..)), color = "black", fill = "skyblue", alpha = 1, size = 0.25) + 
  labs(x = bquote(paste(hat(R)," statistic")), y = NULL)
graph_rhat
ggsave("final_models/rhat.pdf", w = 3, h = 2)  


graph_neff <- graph_diagostics + 
  stat_bin(aes(x = neff_ratio, y=..count../sum(..count..)), color = "black", fill = "skyblue", alpha = 1, size = 0.25) + 
  labs(x = bquote(n[eff]/N), y = NULL)
graph_neff


graph_mcse <- graph_diagostics + 
  stat_bin(aes(x = mcse_ratio, y=..count../sum(..count..)), color = "black", fill = "skyblue", alpha = 1, size = 0.25) + 
  labs(x = bquote(mcse/sd), y = NULL)
graph_mcse


  


# distributions of predictive test statistics ------------------------------------
obs <- ddply(ck_data, "congress", summarise, y = sum(majps))
obs$congress <- as.numeric(obs$congress)

yrep <- extract(stanfit, pars = "y_rep")[[1]]
yreps <- matrix(NA, nrow = nrow(yrep), ncol = length(congress))
for (c in congress) {
  i <- c - min(congress) + 1
  yreps[,i] <- rowSums(yrep[, ck_data$congress == c])
}

rowmeans_yrep <- rowMeans(yreps)
rowsds_yrep <- apply(yreps, 1, sd)
test_stats <- data.frame(means = rowmeans_yrep, sds = rowsds_yrep)
mean_obs <- mean(obs$y)
sd_obs <- sd(obs$y)

graph_test_stats <- ggplot(data = test_stats)  

graph_test_stats_mean <- graph_test_stats + 
  stat_bin(aes(x = means, y=..count../sum(..count..)), color = "black", fill = "skyblue", alpha = 1, size = 0.25) +
  geom_segment(x = mean_obs, xend = mean_obs, y = 0, yend = Inf, size = 2, color = "#222222") + 
  labs(x = bquote(T(y[rep]) == mean(y[rep])), y = NULL) 

graph_test_stats_mean


graph_test_stats_sd <- graph_test_stats + 
  stat_bin(aes(x = sds, y=..count../sum(..count..)), color = "black", fill = "skyblue", alpha = 1, size = 0.25) +
  geom_segment(x = sd_obs, xend = sd_obs, y = 0, yend = Inf, size = 2, color = "#222222") + 
  labs(x = bquote(T(y[rep]) == sd(y[rep])), y = NULL) 

graph_test_stats_sd



# distributions of posterior predictive simulations ----------------------------------------
graphs_yrep <- list()
ids <- sample(nrow(yreps), 11)

graphs_yrep[[1]] <- ggplot(data = data.frame(x = obs$y)) +
  stat_bin(aes(x = x, y=..count../sum(..count..)), color = "black", fill = "darkgray", alpha = 1, size = 0.25) +
  labs(x = bquote(y), y = NULL) 

for (i in seq_along(ids)) {
  graphs_yrep[[i+1]] <- ggplot(data = data.frame(x = yreps[ids[i],])) +
    stat_bin(aes(x = x, y=..count../sum(..count..)), color = "black", fill = "skyblue", alpha = 1, size = 0.25) +
    labs(x = bquote(y[rep]), y = NULL) 
}
do.call(grid.arrange, c(graphs_yrep,nrow = 3))

