library(dplyr)
library(phenoDist)
library(ggplot2)
library(gtools)

################################################################################
# This script demonstrates using the normalised feature values rather than
# the principal components.
# This may be beneficial when we only have a small number of features, and want
# to improve interpretability. Although if you're using high-dimensional
# datasets (20+) features, then it is recommended you use principal components
# as a method to reduce dimensionality. Otherwise you may have issues when
# calculating distances from the negative control.
################################################################################

# load data
# N.B this data is already normalised against the DMSO controls
df <- read.csv("../../data/df_cell_subclass.csv")

# just a single cell line in this example to make things easier
df <- filter(df, Metadata_CellLine == "MDA231")

# select just rows with concentration >= 300 and the DMSO values
df_high_conc <- filter(df, Metadata_concentration >= 300)
df_dmso <- filter(df, Metadata_compound == "DMSO")
df_both <- rbind(df_high_conc, df_dmso)


# center data so that DMSO centroid is at the origin
df_centered <- centre_control(df_both,
                              cols = 1:260, # TODO check this!
                              cmpd_col = "Metadata_compound",
                              cmpd = "DMSO")

# function to compute the distance (L1 norm) from the centroid to each point
l1_dist_from_origin <- function(x) {
    origin <- rep(0, length(x))
    sum(abs(x - origin))
}

# Now we want to determine the distance from the origin for the negative control
# values to determine a cutoff for inactive compounds

dmso_dist <- apply(filter(df_centered, Metadata_compound == "DMSO")[, 1:260], 1,
                   l1_dist_from_origin)

cutoff <- mean(dmso_dist) + sd(dmso_dist)

df_centered$Metadata_dist_from_origin <- apply(df_centered[, 1:260], 1,
                                               l1_dist_from_origin)

controls <- c("DMSO", "STS")
`%notin%` <- function(x, y) !(x %in% y)

df_cmpd <- filter(df_centered,
                  Metadata_compound %notin% controls &
                  Metadata_dist_from_origin > cutoff)

# refactor to remove empty classes (i.e STS)
df_cmpd$Metadata_compound <- factor(df_cmpd$Metadata_compound)


# select just feature data and compound labels
df_cmpd <- df_cmpd[, c(1:260, 266)]

head(df_cmpd)

df_agg <- group_by(df_cmpd, Metadata_compound) %>%
    summarise_each(funs(mean))

# function to calculate theta in high-dimensions
cosine_pairs <- function(x, cols, degrees = TRUE) {

    if (!is.list(x) || is.data.frame(x)) {
        stop(paste(substitute(x), "should be a list"),
             call. = FALSE)
    }

    # initialise empty vectors
    # TODO pre-allocate vector lengths
    vals <- numeric()
    A <- character()
    B <- character()

    # get pairs of compounds
    # want to get compounds against itself
    name <- names(x)
    pairs_names <- combinations(n = length(name), 2, name,
                                repeats.allowed = TRUE)

    for (i in 1:nrow(pairs_names)) {

        tmp1 <- x[[pairs_names[i, 1]]]
        tmp2 <- x[[pairs_names[i, 2]]]

        # store values in a matrix
        mat <- matrix(nrow = nrow(tmp1),
                      ncol = nrow(tmp2))

        # loop through rows in cmpd A and cmpd B
        # calculate the cosine similarity between the two vectors
        # store values as co-ordinates in the matrix
        for (j in 1:nrow(tmp1)) {
            for (k in 1:nrow(tmp2)) {
                mat[j, k] <- cosine_sim_vector(tmp1[j, cols],
                                               tmp2[k, cols])
            }

        }

        val <- as.vector(mat)
        vals <- append(vals, val)
        A <- append(A, rep(pairs_names[i, 1], length(val)))
        B <- append(B, rep(pairs_names[i, 2], length(val)))

    }

    # convert cosine similarity to 0 -> 180 degrees (biologist friendly)
    if (degrees) {
        vals <- (1 - cossim_to_angsim(vals)) / (1/180)
    }

    data.frame(A, B, vals)
}

df_split <- split(df_agg, df_agg$Metadata_compound)

cosine_out <- cosine_pairs(df_split, cols = 1:260)

ggplot(data = cosine_out,
      aes(x = A, y = B,
          fill = vals)) +
    geom_tile() +
    coord_fixed() +
    viridis::scale_fill_viridis("theta") +
    Smisc::xlab_rotate() +
    xlab("") + ylab("") +
    ggtitle("Theta values between compounds, MDA-MB-231, using feature values")
ggsave("../../figures/feautures_plot.pdf")
