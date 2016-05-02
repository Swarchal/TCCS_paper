library(ggplot2)
library(Smisc)
library(viridis)

df_delta_theta <- readRDS("data/df_delta_theta")

# plot heatmap of differences between compounds
ggplot(data = df_delta_theta,
       aes(x = A, y = B, fill = difference)) +
    geom_tile() + 
    coord_fixed() + 
    xlab_rotate() +
    facet_wrap(~drug, ncol = 3) + 
    scale_fill_viridis(bquote(theta)) +
    xlab("") +
    ylab("")
ggsave("figures/delta_theta_heatmap.eps", width = 6, height = 9)