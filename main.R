files <- list.files("analysis")
source("load_packages.R")
for (file in files){
    cat(paste("- Running", file, "\n"))
    source(paste("analysis", file, sep = "/"))
}
cat(paste(" - DONE", "\n"))