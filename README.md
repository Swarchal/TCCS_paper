The R scripts assume the use of RStudio and the .rproj file to correctly set the working directory.

To use without RStudio, the working directory has to be set to `<users_path_to_file>/cosine_paper2/` for the file paths to work correctly.

To run analysis use the `main.R` function, either within RStudio or from the command line with `Rscript main.R`  which first checks and installs the necessary packages from CRAN and GitHub, and then runs the .R files in `analysis/` to generate the plots used in the paper, which are saved as .eps files in `figures/`.

NB. The CellProfiler pipelines will not run, due to the lack of hosted images.
