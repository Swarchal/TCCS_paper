library(phenoScreen)
library(dplyr)
library(ggplot2)

# create circular histogram for MCF7, KPL4, SN38 data
# used in figure 4.
# hopefully looks different between the two as they have theta of 179

df <- readRDS("data/df_delta_theta")

# subset SN38 data

df_sn38 <- df %>% filter(drug == "sn38",
                         B == "KPL4",
                         A == "MCF7")

# turns out I have two values because I've summarised all the concentration data
# and replicates to a single value

# need to use the original data and calculate a theta value for each replicate
# and plot that on the rose plot

# load original dataframe
df_new <- read.csv("data/df_cell_subclass.csv")

# calculate PCA on all the data before subsetting
pca_out <- prcomp(df_new[, get_featuredata(df_new)])$x[,1:2]

pca_df <- data.frame(pca_out,
                     df_new[, -get_featuredata(df_new)])


# subset SN38 at 1 uM
cell_lines <- c("MCF7", "KPL4")
df_want <- pca_df %>% filter(Metadata_compound == "sn38",
                             Metadata_CellLine %in% cell_lines,
                             Metadata_concentration == 1000)

# calculate theta against the reference vector
df_want$theta <- NA
for (i in 1:nrow(df_want)){
    df_want$theta[i] <- theta0(c(df_want$PC1[i], df_want$PC2[i]))
}



# plot data
ggplot(data = df_want,
       aes(x = theta,
           group = Metadata_CellLine)) +
    geom_histogram(binwidth = 15,
                   aes(fill = Metadata_CellLine)) +
    coord_polar(start = -1.57, direction = -1) +
    scale_x_continuous(breaks = seq(0, 360, by = 45), expand = c(0,0), lim = c(0, 360)) +
    scale_size_area() +
    geom_text(data = NULL, size = 4, x = 175, y = 15,
              label = paste("theta =", format(round(theta_out, 2), nsmall = 2))) +
    scale_fill_brewer(name = "Cell line", palette = "Dark2") + 
    xlab("") + ylab("") + 
    ggsave("figures/directional_sn38.eps", height = 7, width = 5)