library(viridis)
library(dplyr)
library(gplots)
library(phenoDist)
library(phenoScreen)

df <- read.csv("data/df_cell_subclass.csv")

# first two principal components
pca <- prcomp(df[, get_featuredata(df)])$x[,1:2]

df2 <- data.frame(pca,
                  df[, grep("Metadata", colnames(df))])

pca_subset <- df2

# calculate theta value for each group
# group by compound & cell-line and get average PC1:2 values
pca_vectors <- group_by(pca_subset, Metadata_compound, Metadata_CellLine) %>% 
    summarise_each(funs(mean), PC1, PC2)

pca_vectors$theta <- NA

# calculate a theta value for each row
# each row will have a unique compound -- celll-ine combination
for (row in 1:nrow(pca_vectors)){
    pca_vectors$theta[row] <- theta0(c(pca_vectors$PC1[row],
                                       pca_vectors$PC2[row]))
}


pca_vectors <- filter(pca_vectors, Metadata_compound == "saracatinib")

data = pca_vectors$theta
x = pca_vectors$Metadata_CellLine

out <- sapply(data, function(x, y = data){abs(x - y)})

# add compound levels to the matrix
row.names(out) <- x
colnames(out) <- x

# if values are greater than 180, then they become more similar again
# 180 degrees should be the maximum values as this signifies completely opposite directions
# values close to 360 will be very similar
# horrible hack to get around this:
fold_180 <- function(x){
    out <- ifelse(x > 180, 360 - x, x)
    return(out)
}

# apply function to every cell in the matrix
out_mess <- apply(out, 1:2, fold_180)

# plot and save figure
setEPS()
postscript("figures/saracatinib_diff.eps", width = 7, height = 7)
heatmap.2(out_mess,
          margins = c(8, 7),
          col = viridis,
          trace = "none",
          keysize = 1.1,
          key.xlab = expression(theta),
          density.info = "none",
          key.title = "none")
dev.off()
