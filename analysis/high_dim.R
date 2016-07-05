library(dplyr)
library(phenoDist)
library(phenoScreen)
library(gtools)
library(ggplot2)


###############################################################################
# Example workflow showing how to use the TCCS method with high-dimensional data.
# This may be needed when the dataset contains many measured features, so that
# two principal components do not represent a large enough proportion of the
# variance contained in the data, or when screening a large number of compounds
# that require greater discriminating power to separate them from each other
###############################################################################


# load data
df <- read.csv("data/df_cell_subclass.csv")

# just a single cell line in this example to make things more simple
df <- filter(df, Metadata_CellLine == "MDA231")

# select just those rows with a concentration > 300nM, and the DMSO values
df_high_conc <- df %>% filter(Metadata_concentration >= 300)
df_dmso <- df %>% filter(Metadata_compound == "DMSO")
df_both <- rbind(df_high_conc, df_dmso)

# calculate principal components
pca <- prcomp(df_both[, get_featuredata(df_both)])


# create dataframe of principal components and metadata
df2 <- data.frame(pca$x,
                  df_both[, grep("Metadata", colnames(df_both))])


# center principal components so that DMSO centroid is at (0,0)
df_centered <- centre_control(df2,
                              cols = 1:260,
                              cmpd_col = "Metadata_compound",
                              cmpd = "DMSO")


#-----------------------------------------------------------------------------
# For a more unbiased approach at selecting the dimensionality of the dataset,
# it's common to use the number of principal components that capture a specified
# proportion of the variance in the data.
#
# In this case we demonstrate using the minumum number of principal components
# that capture 70% of the variance in out dataset


# calculate number of principal components to capture 70% of the variance
n_pc <- min(which(cumsum(pca$sdev^2) / sum(pca$sdev^2) >= 0.7))

pdf("figures/cumulative_PC_variance.pdf")
plot(cumsum(pca$sdev^2) / sum(pca$sdev^2),
     type = "l", lwd = 3, col = "cornflowerblue",
     main = "Cumulative variance\ndescribed by principal components",
     xlab = "Principal component",
     ylab = "Prop. variance")
dev.off()
#---------------------------
# Need to remove points within a minimum distance of the DMSO centroid
# distance of each point from the DMSO centroid?
# points are now centered at dmso centroid due to `centre_control()`, so we can
# calculate the distance from the origin to each point. Though in high-dimensional
# space l1 norm (Manhattan distance) might be more sensible than Euclidean distance
#---------------------------

l1_dist_from_origin <- function(x) {
    origin <- rep(0, length(x))
    sum(abs(x - origin))
}

# l1 norm for DMSO controls, find SD and define cutoff
dmso_dist <- apply(filter(df_centered, Metadata_compound == "DMSO")[, 1:260], 1,
                   l1_dist_from_origin)


# calculate the minimum cutoff distance from the negative control
cutoff <- mean(dmso_dist) + sd(dmso_dist)

# add column of distance from origin to subset by
df_centered$Metadata_dist_from_origin <- apply(df_centered[, 1:260], 1,
                                               l1_dist_from_origin)

controls <- c("DMSO", "STS")
`%notin%` = function(x,y) !(x %in% y)
df_cmpd <- filter(df_centered,
                  Metadata_compound %notin% controls & Metadata_dist_from_origin > cutoff)

# don't subset yet, want to keep all dimensionality to run tests
df_sub <- data.frame(df_cmpd[, 1:260],
                     df_cmpd[, grep("Metadata", colnames(df_cmpd))])


# refactor to remove empty classes (i.e STS)
df_sub$Metadata_compound <- factor(df_sub$Metadata_compound)


df_agg <- df_sub %>% group_by(Metadata_compound) %>%
    summarise_each(funs(mean))

# re-order columns (such a pain in R)
df_agg <- df_agg %>% select(-Metadata_compound, everything())


# function to calculate cosine similarity between pairs of vectors
# in N-dimensional space
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

    # convert cosine similarity to 0 -> 180 degrees
    if (degrees) {
        vals <- (1 - cossim_to_angsim(vals)) / (1/180)
    }

    data.frame(A, B, vals)
}


# split df_sub into a list, each element containing a compound
df_split <- split(df_agg, df_agg$Metadata_compound)

# output for 70% variance dimensionality
cosine_out <- cosine_pairs(df_split, cols = 1:n_pc)
# output for dimensionality = 2
cosine_out_n2 <- cosine_pairs(df_split, cols = 1:2)
# output for dimensionality = 10
cosine_out_n10 <- cosine_pairs(df_split, cols = 1:10)
# ouput for dimensionality = 260
cosine_out_nall <- cosine_pairs(df_split, cols = 1:260)


ggplot(data = cosine_out,
      aes(x = A, y = B,
          fill = vals)) +
    geom_tile() +
    coord_fixed() +
    viridis::scale_fill_viridis("theta") +
    Smisc::xlab_rotate() +
    xlab("") + ylab("") +
    ggtitle("n PC = 70% variance")
ggsave("figures/PC70.pdf")

ggplot(data = cosine_out_n2,
      aes(x = A, y = B,
          fill = vals)) +
    geom_tile() +
    coord_fixed() +
    viridis::scale_fill_viridis("theta") +
    Smisc::xlab_rotate() +
    xlab("") + ylab("") +
    ggtitle("n PC = 2")
ggsave("figures/PC2.pdf")

ggplot(data = cosine_out_n10,
      aes(x = A, y = B,
          fill = vals)) +
    geom_tile() +
    coord_fixed() +
    viridis::scale_fill_viridis("theta") +
    Smisc::xlab_rotate() +
    xlab("") + ylab("") +
    ggtitle("n PC = 10")
ggsave("figures/PC10.pdf")

ggplot(data = cosine_out_nall,
      aes(x = A, y = B,
          fill = vals)) +
    geom_tile() +
    coord_fixed() +
    viridis::scale_fill_viridis("theta") +
    Smisc::xlab_rotate() +
    xlab("") + ylab("") +
    ggtitle("n PC = all")
ggsave("figures/PCall.pdf")