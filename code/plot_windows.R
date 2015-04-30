# Plot results from Cox and Katz's method for different sized rolling windows

library(ggplot2)
library(plyr)
library(reshape2)
library(mvtnorm)

load("ck_data.RData")


ests <- function(data, t, N = 1e4) {
  congresses <- unique(data$congress)
  out <- sapply(congresses, simplify = FALSE, FUN = function(C) {
    model <- glm(cbind(majps, nvotes - majps) ~ lnmajvavg, 
                 data = data, family = binomial,
                 subset = congress >= C - t & congress <= C + t)
    betas <- rmvnorm(N, mean = coef(model), sigma = vcov(model))
    betas[,1] <- plogis(betas[,1]) - 0.5
    return(betas)
  })
  
  meds <- rep(NA, length(congresses))
  for(i in seq_along(out)) meds[i] <- median(out[[i]][,"(Intercept)"])
  
  meds
}

windows <- 0:5
results <- lapply(windows, function(t) ests(ck_data, t))
bias <- do.call(c, results)

congresses <- sort(unique(ck_data$congress))
nCong <- length(congresses)
window <- rep(c("t %+-% 0", "t %+-% 1", "t %+-% 2", "t %+-% 3", "t %+-% 4", "t %+-% 5"), each = nCong)
dfm <- melt(data.frame(Congress = congresses, Bias = bias, Window = window), 
            id.vars = c("Congress","Window"))

axis_line_color <- "#222222"
theme_set(theme_classic())
theme_update(
  axis.line.x = element_line(size = 3, color = axis_line_color),
  axis.line.y = element_line(size = 0.5, color = axis_line_color),
  axis.title = element_text(face = "bold", size = 13),
  strip.background = element_rect(fill = "gray35", color = "gray35"),
  strip.text = element_text(color = "white", face = "bold", size = 15),
  legend.position = "none"
)

base <- ggplot(dfm, aes(x = Congress, y = value, group = Window))
graph <- base + 
  stat_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "gray", size = .25, alpha = 0.5) + 
  geom_line(aes(color = Window), size = 1, show_guide = FALSE) + 
  geom_point(color = "gray35", size = 1.5) +
  facet_grid(Window~., labeller = label_parsed) + 
  ylab("Bias toward majority") 

graph



