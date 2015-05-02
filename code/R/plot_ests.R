library(rstan)
library(ggplot2)
library(plyr)
library(reshape2)
library(gridExtra)

axis_line_color <- "#222222"
theme_set(theme_classic())
theme_update(
  axis.line.x = element_line(size = 3, color = axis_line_color),
  axis.line.y = element_line(size = 0.5, color = axis_line_color),
  legend.title = element_text(face = "bold"),
  legend.text = element_text(size = 18),
  plot.title = element_text(size = 11),
  strip.background = element_blank(), 
  strip.text = element_blank()
)


posterior <- extract(stanfit, pars = c("lambda", "rho"))
probs <- c(0.025, .25, 0.5, 0.75, 0.975)
quants <- function(posterior, param) t(apply(posterior[[param]], 2, quantile, probs))
lambda <- quants(posterior, "lambda")
rho <- quants(posterior, "rho")
colnames(lambda) <- colnames(rho) <- c("lb95", "lb50", "est", "ub50", "ub95")

lambda <- data.frame(cbind(Congress = 46:106, lambda))
rho <- data.frame(cbind(Congress = 46:106, rho))

df <- data.frame(Congress = 46:106, rbind(lambda, rho), 
                 parameter = rep(c("lambda", "rho"), each = nrow(lambda)))
df$parameter <- factor(df$parameter, levels = c("rho", "lambda"))
param_labels <- list(bquote(rho), bquote(lambda))

base <- ggplot(df, aes(x = Congress, ymin = lb95, ymax = ub95, y = est, fill = parameter))
graph <- base + geom_ribbon(alpha = 0.25) + 
  geom_ribbon(aes(ymin = lb50, ymax = ub50), alpha = 0.65) +
  geom_line(aes(color = parameter)) + 
  geom_point(aes(color = parameter)) + 
  scale_fill_brewer(name = "", type = "qual", labels = param_labels) + 
  scale_color_brewer(name = "", type = "qual", labels = param_labels) +
  ylab("Parameter value") 
  
graph + facet_grid(parameter ~., scales = "free")   

