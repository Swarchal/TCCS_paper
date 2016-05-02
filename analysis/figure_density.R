library(ggplot2)

df_delta_theta <- readRDS("data/df_delta_theta")

# plot density plot of differences
# vast majaority of differences between cell lines are small
ggplot(data = df_delta_theta,
       aes(x = difference)) + 
    geom_density(fill = "white") + 
    xlab(expression(theta))
ggsave("figures/density_delta_theta.eps", width = 5, height = 4)