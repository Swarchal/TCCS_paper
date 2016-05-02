library(dplyr)
library(ggplot2)

df <- read.csv("data/df_cell_subclass.csv")

df_sara <- filter(df, Metadata_compound == "saracatinib", Metadata_concentration >= 100)
df_sub <- filter(df_sara, Metadata_CellLine %in% c("SKBR3", "KPL4", "MDA231"))
df_sub$Metadata_CellLine <- droplevels(df_sub$Metadata_CellLine)
# reorder factor levels
df_sub$Metadata_CellLine <- factor(df_sub$Metadata_CellLine, levels(df_sub$Metadata_CellLine)[c(3,1,2)])

ggplot(data = df_sub,
       aes(x = Metadata_CellLine,
           y = cell_Neighbors_NumberOfNeighbors_2)) + 
    geom_boxplot(fill = "darkorange") + 
    geom_jitter(width = 0.4) + 
    xlab("Cell Line") +
    ylab("Z-score Number of Neighbours (cells)") 
ggsave("figures/saracatinib_features.eps", width = 5, height = 4)