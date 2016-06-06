library(ggplot2)
library(dplyr)
library(phenoScreen)
library(phenoDist)
library(Smisc)
library(caret)
library(reshape2)
library(viridis)
library(dplyr)

# load data
df <- read.csv("data/df_cell_subclass.csv")

# principal components of the feature data columns
pca <- prcomp(df[, get_featuredata(df)])

# figure of cumulative variance of principal components
setEPS()
postscript("figures/pca_variance_explained.eps", width = 8, height = 6)
var_expl <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
prop_var <- 0.85
n_components <- min(which(var_expl >= prop_var))
plot(var_expl, type = "l",
     xlab = "principal components",
     ylab = "cumulative variance")
segments(x0 = n_components,
         y0 = -10,
         y1 = prop_var, lty = 2)
segments(x0 = -10,
         x1 = n_components,
         y0 = prop_var,
         lty = 2)
dev.off()


# create dataframe of the first 2 principal components and metadata
pca_df <- data.frame(pca$x[,0:2], # first 2 prin comps
                     df[, grep("Metadata_", colnames(df))]) # metadata


# calculate the multivariate z-factor for each cell lines
pca_df_z <- data.frame(pca$x,
                       df[, grep("Metadata_", colnames(df))])

cl_z_factor <- sapply(split(pca_df_z, pca_df_z$Metadata_CellLine),
                      function(x){
                          multi_z(x,
                                  feature_cols = get_featuredata(x),
                                  cmpd_col = "Metadata_compound",
                                  pos = "STS",
                                  neg = "DMSO")})

# dataframe of cell lines and z-factor values
cl_z_df <- data.frame(cell_line = rownames(data.frame(cl_z_factor)),
                      z_prime = cl_z_factor)


# sort by values of z_prime
cl_z_df <- transform(cl_z_df, cell_line = reorder(cell_line, - z_prime))

# dotchart of z-prime values
ggplot(data = cl_z_df,
       aes()) +
    geom_segment(aes(x = 0,
                     xend = z_prime,
                     y = cell_line,
                     yend = cell_line),
                 col = "gray40") +
    geom_point(aes(z_prime, cell_line), size = 2.5) +
    xlab("multivariate Z'") +
    ylab("") +
    theme(axis.text.y = element_text(face = "bold"))
ggsave("figures/z_factor.eps", width = 6, height = 4)
#########################################################################



# centre principal components so that the DMSO centroid is centered
# at co-ordinates 0,0
pca_df <- centre_control(pca_df,
                         cols = get_featuredata(pca_df),
                         cmpd_col = "Metadata_compound",
                         cmpd = "DMSO")



# euclidean distance function
distance <- function(x, y){
    dist <- sqrt(x^2 + y^2)
    return(dist)
}

# calculate norm (length) of each vector
pca_df$dist <- NA
for (row in 1:nrow(pca_df)){
    pca_df$dist[row] <- distance(pca_df$PC1[row], pca_df$PC2[row])
}


# select a single cell line
df_mda231 <- filter(pca_df, Metadata_CellLine == "MDA231")

# select single compound data within that cell line
df_mda231_barasertib <- filter(df_mda231, Metadata_compound == "barasertib")


# scatter plot of first 2 principal components
# barasertib datapoints coloured by concentration
ggplot() +
    geom_point(data = df_mda231,
               colour = "gray50",
               aes(x = PC1,
                   y = PC2)) +
    geom_point(size = 3,
               data = df_mda231_barasertib,
               aes(x = PC1,
                   y = PC2,
                   colour = Metadata_concentration)) +
    geom_line(data = df_mda231_barasertib,
              size = 1,
              aes(x = PC1,
                  y = PC2,
                  colour = Metadata_concentration)) +
    scale_color_viridis(name = "Concentration (nM)",
                        trans = "log10")
ggsave("figures/increasing_barasertib_mda231.eps", width = 8, height = 6)


compound <- "cycloheximide"
# select a single cell line
df_mda231 <- filter(pca_df, Metadata_CellLine == "MDA231")

# select single compound data within that cell line
df_mda231_compound <- filter(df_mda231, Metadata_compound == compound)

# scatter plot of first 2 principal components of cycloheximide
# points coloured by concentration
ggplot() +
    geom_point(data = df_mda231,
               colour = "gray50",
               aes(x = PC1,
                   y = PC2)) +
    geom_point(size = 3,
               data = df_mda231_compound,
               aes(x = PC1,
                   y = PC2,
                   colour = Metadata_concentration)) +
    geom_line(data = df_mda231_compound,
              size = 1,
              aes(x = PC1,
                  y = PC2,
                  colour = Metadata_concentration)) +
    scale_color_viridis(name = "Concentration (nM)",
                        trans = "log10")
ggsave("figures/increasing_cycloheximide_mda231.eps", width = 8, height = 6)



pca_df$theta <- NA # initialise empty column for loop

# loop through rows of data calculating theta for each vector(PC1, PC2)
for (i in 1:nrow(pca_df)){
    pca_df$theta[i] <- theta0(c(pca_df$PC1[i], pca_df$PC2[i]))
}


# filter just barasertib data
df_barasertib <- filter(pca_df, Metadata_compound == "barasertib")

# circular hisotgram of batasertib theta values
ggplot(data = df_barasertib,
       aes(x = theta,
           group = Metadata_concentration)) +
    geom_histogram(binwidth = 15,
                   aes(fill = Metadata_concentration)) +
    coord_polar(start = -1.57, direction = -1) +
    scale_x_continuous(breaks = seq(0, 360, by = 45), expand = c(0,0), lim = c(0, 360)) +
    scale_size_area() +
    xlab("") + ylab("") +
    scale_fill_viridis(name = "concentration (nM)",
                       trans = "log10")
ggsave("figures/directional_histogram_barasertib.eps", width = 8, height = 6)


# circular histogram of batasertib theta values
# small plot for each cell line
ggplot(data = df_barasertib,
       aes(x = theta,
           group = Metadata_concentration)) +
    geom_histogram(binwidth = 15,
                   aes(fill = Metadata_concentration)) +
    coord_polar(start = -1.57, direction = -1) +
    scale_x_continuous(breaks = seq(0, 360, by = 45), expand = c(0,0), lim = c(0, 360)) +
    scale_size_area() +
    scale_fill_viridis(name = "concentration (nM)",
                       trans = "log10") +
    xlab("") + ylab("") +
    facet_wrap(~Metadata_CellLine, ncol = 2) +
    theme(axis.text.x = element_text(size = 6))
ggsave("figures/directional_histogram_barasertib_split.eps", width = 8, height = 12)


# filter barasertib and monastrol data
wanted_compounds <- c("monastrol", "barasertib")
df_two <- filter(pca_df, Metadata_compound %in% wanted_compounds, Metadata_CellLine == "MDA231")
df_two$Metadata_compound <- relevel(df_two$Metadata_compound, "cycloheximide")

# circular histogram of barasetib and monastrol data
ggplot(data = df_two,
       aes(x = theta,
           group = Metadata_compound)) +
    geom_histogram(binwidth = 15,
                   aes(fill = Metadata_compound)) +
    coord_polar(start = -1.57, direction = -1) +
    scale_x_continuous(breaks = seq(0, 360, by = 45), expand = c(0,0), lim = c(0, 360)) +
    scale_size_area() +
    scale_fill_brewer(name = "compound", palette = "Set2")
ggsave("figures/barasertib_monastrol_hist.eps", width = 8, height = 6)


# function to calculate the average vector from replicates
average_vector <- function(dat){
    # data will be in the form of rows = vectors, columns = components
    means <- as.vector(colMedians(dat))
    return(means)
}


# filter just barasertib data from MDA-231 cell line
barasertib_data <- filter(pca_df, Metadata_compound == "barasertib" &
                              Metadata_CellLine == "MDA231")

# calculate the average vector from the replicates of baraserib in MDA231
# from the vector(PC1, PC2)
vector_info_barasertib <- matrix(c(barasertib_data$PC1, barasertib_data$PC2), ncol = 2)
vector_barasertib <- average_vector(vector_info_barasertib)

# filter monastrol data from MDA-231 cell line
monastrol_data <- filter(pca_df, Metadata_compound == "monastrol" &
                             Metadata_CellLine == "MDA231")

# calculate the average vector from the replicates of monastrol in MDA231
# from the vector(PC1, PC2)
vector_info_monastrol <- matrix(c(monastrol_data$PC1, monastrol_data$PC2), ncol = 2)
vector_monastrol <- average_vector(vector_info_monastrol)

# calculate theta between two the averaged vectors of baraserib and monastrol
theta_out <- theta(vector_barasertib, vector_monastrol)


# circular histogram of theta values from monastrol and barasertib
# this time labelled with the average vector and calcualted theta value
# between the two vectors
ggplot(data = df_two,
       aes(x = theta,
           group = Metadata_compound)) +
    geom_histogram(binwidth = 15,
                   aes(fill = Metadata_compound)) +
    coord_polar(start = -1.57, direction = -1) +
    scale_x_continuous(breaks = seq(0, 360, by = 45), expand = c(0,0), lim = c(0, 360)) +
    scale_size_area() +
    geom_vline(xintercept = theta0(vector_monastrol)) +
    geom_vline(xintercept = theta0(vector_barasertib)) +
    geom_text(data = NULL, size = 4, x = 225, y = 10,
              label = paste("theta =", format(round(theta_out, 2), nsmall = 2))) +
    scale_fill_brewer(name = "compound", palette = "Set2")
ggsave("figures/barasertib_monastrol_hist_ann.eps", width = 8, height = 6)


# calculate theta value between two cell lines (MDA231 & HCC1569) from their
# average PC1/2 vectors
# select barasertib data for MDA231 and HCC1569 lines:
data_comp_cells <- filter(pca_df, Metadata_compound == "barasertib",
                          Metadata_CellLine == "MDA231" | Metadata_CellLine == "HCC1569")

# mean vector MDA231
just_mda <- filter(data_comp_cells, Metadata_CellLine == "MDA231")
vector_mda <- average_vector(matrix(c(just_mda$PC1, just_mda$PC2), ncol = 2))

# mean vector HCC1569
just_hcc <- filter(data_comp_cells, Metadata_CellLine == "HCC1569")
vector_hcc <- average_vector(matrix(c(just_hcc$PC1, just_hcc$PC2), ncol = 2))

# theta value between the 2 cell line's averaged vectors
theta_out <- theta(vector_mda, vector_hcc)

# circular histogram of MDA-231 and HCC1569 treated with barasertib, with
# labelled average vectors and theta value between the two cell lines
ggplot(data = data_comp_cells,
       aes(x = theta,
           group = Metadata_CellLine)) +
    geom_histogram(binwidth = 15,
                   aes(fill = Metadata_CellLine)) +
    coord_polar(start = -1.57, direction = -1) +
    scale_x_continuous(breaks = seq(0, 360, by = 45), expand = c(0,0), lim = c(0, 360)) +
    scale_size_area() +
    geom_vline(xintercept = theta0(vector_mda)) +
    geom_vline(xintercept = theta0(vector_hcc)) +
    geom_text(data = NULL, size = 4, x = 175, y = 15,
              label = paste("theta =", format(round(theta_out, 2), nsmall = 2))) +
    scale_fill_brewer(name = "Cell line", palette = "Pastel1")
ggsave("hcc1569_231_hist_ann.eps", width = 8, height = 6)

##########################################
#---          Cosine analysis         ---#
##########################################


# filter single concentration (100nM)
concentration <- 1000
df_1000 <- filter(df, Metadata_concentration == concentration)
controls <- c("DMSO", "STS")
# get control data (DMSO & staurosporine)
df_dmso <- filter(df, Metadata_compound %in% controls)

# row-bind 100nM and control data into a single dataframe
df_new <- rbind(df_1000, df_dmso)

# calculate the first two principal components of the featuredata
pca_out <- prcomp(df_new[, get_featuredata(df_new)])$x[,1:2]
# create dataframe of first two principal components and metadata
pca_df <- data.frame(pca_out,
                     df_new[, grep("Metadata_", colnames(df))])


# calculate theta and vector distances
# initialise empty columns for loop
pca_df$theta <- NA
pca_df$vector_norm <- NA
# loop through rows of principal components, calculating the theta value
# against a place-holder vector (1, 0), and the norm from the origin
for (i in 1:nrow(pca_df)){
    pca_df$theta[i] <- theta0(c(pca_df$PC1[i], pca_df$PC2[i]))
    pca_df$vector_norm[i] <- norm_vector(c(pca_df$PC1[i], pca_df$PC2[i]))
}

# create %notin% function
`%notin%` <- function(x, y) !(x %in% y)

# create cutoff constants
cutoff_n <- 1
max_cutoff_n <- 100

# calculate cutoff from standard deviations of the norms from the origin
cutoff <- cutoff_n * sd(pca_df$vector_norm)
max_cutoff <- max_cutoff_n * sd(pca_df$vector_norm)

pca_df$cutoff <- paste("<", cutoff_n)
pca_df$cutoff[pca_df$vector_norm > cutoff] <- paste(max_cutoff_n, "> x >", cutoff_n)
pca_df$cutoff[pca_df$vector_norm > max_cutoff] <- paste(">", max_cutoff_n)

# scatter plot of first two principal components, coloured by whether they are
# beyond the cut-off value or not
ggplot(data = pca_df,
       aes(x = PC1,
           y = PC2,
           col = as.factor(cutoff))) +
    geom_point() +
    coord_fixed() +
    scale_color_brewer("Standard deviations", palette = "Set1")
ggsave("figures/cutoff.eps", width = 8, height = 6)


# unwanted control data labels
unwanted <- c("DMSO", "STS")

cell_lines <- c("MDA231",
                "SKBR3",
                "MDA157",
                "T47D",
                "KPL4",
                "MCF7",
                "HCC1569",
                "HCC1954")


# function to filter rows that are beyond the cutoff for each cell line
# in the compound data
cell_line_cutoff <- function(x){
    filter(pca_df, Metadata_CellLine == x,
           vector_norm > cutoff,
           Metadata_compound %notin% unwanted) %>%
        distinct(Metadata_compound)
}

# convert cell-line names to lower case for eval() to match variable names
for (i in cell_lines){
    assign(tolower(i),
           cell_line_cutoff(i))
}

# as compounds may be within the cut-off in some cell lines and beyond the
# cut-off in other cell lines, only the compounds that are beyond the cutoff in
# all eight cell lines are used in the futher analyses
common_compounds <- Reduce(intersect, list(mda231$Metadata_compound,
                                           skbr3$Metadata_compound,
                                           mda157$Metadata_compound,
                                           t47d$Metadata_compound,
                                           kpl4$Metadata_compound,
                                           mcf7$Metadata_compound,
                                           hcc1569$Metadata_compound,
                                           hcc1954$Metadata_compound))

filter_common_compounds <- function(x){
    filter(eval(parse(text = x)), Metadata_compound %in% common_compounds)
}

# convert cell lines names back to lower-case (again)
for (i in cell_lines){
    assign(tolower(i),
           filter_common_compounds(tolower(i)))
}


th_A <- mda231$theta
th_B <- kpl4$theta

out_test <- sapply(th_A, function(x, y = th_B){abs(x - y)})

out_test <- apply(out_test, 1:2, fold_180)
dimnames(out_test) <- list(mda231$Metadata_compound, mda231$Metadata_compound)

# can use diag() to extract the diagonal of the matrix, which returns the angle
# between the drugs between the two cell-lines

diag_out <- as.data.frame(diag(out_test))
diag_out <- cbind(drug = rownames(diag_out), diag_out)
rownames(diag_out) <- NULL
names(diag_out)[2] <- "difference"

diag_out$drug <- with(diag_out, reorder(drug, difference))


cell_lines <- c("mda231",
                "skbr3",
                "mda157",
                "t47d",
                "kpl4",
                "mcf7",
                "hcc1569",
                "hcc1954")


# dataframe of all combinations of cell-lines:
clb <- expand.grid(cell_lines, cell_lines)

# PITA factors
clb <- sapply(clb, as.character)

# start empty data frame to place results into
# will contain all combinations of cell-lines, drugs and their dt values
df_delta_theta <- data.frame(A = NA,
                             B = NA,
                             drug = NA,
                             difference = NA)

# function to find delta-theta values between two cell-lines
find_delta_theta <- function(a, b){
    a_ <- get(a)
    b_ <- get(b)
    th_A <- a_$theta
    th_B <- b_$theta

    out_test <- sapply(th_A, function(x, y = th_B){abs(x - y)})
    out_test <- apply(out_test, 1:2, fold_180)
    # compound vectors are identical across all cell lines
    # use any one of them (mda231 in this case)
    dimnames(out_test) <- list(mda231$Metadata_compound, mda231$Metadata_compound)

    # can use diag() to extract the diagonal of the matrix, which returns the angle
    # between the drugs between the two cell-lines
    diag_out <- as.data.frame(diag(out_test))
    diag_out <- cbind(drug = rownames(diag_out), diag_out)
    rownames(diag_out) <- NULL
    names(diag_out)[2] <- "difference"

    diag_out$A <- eval(substitute(a))   # add cell-line name
    diag_out$B <- eval(substitute(b))   # add cell-line name

    # refactor 'drug' so in numerical order according to difference
    diag_out$drug <- with(diag_out, reorder(drug, difference))
    diag_out
}

# loop through all possible combinations of cell-lines and store as a list of
# dfs
list_of_df <- list()
for (row in 1:nrow(clb)){
    list_of_df[row] <- list(find_delta_theta(clb[row, 1], clb[row, 2]))
}

# row-wise bind of list into single df
df_delta_theta <- do.call(rbind, list_of_df)

# make cell lines uppercase for figures
df_delta_theta[, 3:4] <- apply(df_delta_theta[, 3:4], 2, toupper)

saveRDS(df_delta_theta, file = "data/df_delta_theta")
