# install necessary packages

load_packages <- function(){
    # checks if packages are already installed,
    # installs packages if needed
    
    cran_packages <- c(
        "devtools",
        "dplyr",
        "reshape",
        "reshape2",
        "readr",
        "caret",
        "viridis",
        "gplots",
        "ggplot2"
    )
    
    github_packages <- c(
        "phenoScreen",
        "phenoDist",
        "focus",
        "Smisc"
    )
    
    for (package in cran_packages){
        if (!require(package, character.only = TRUE)){
            install.packages(package)
        }
    }
    
    for (g_package in github_packages){
        if (!require(g_package, character.only = TRUE)){
            devtools::install_github(paste("swarchal", g_package, sep = "/"))
        }
    }
}

load_packages()