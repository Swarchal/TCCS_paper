# TCCS (Theta Comparative Cell Scoring)

Data and code to generate figures from the paper: *Warchal, Dawson and Carragher -  'Development of the Theta-Comparative-Cell-Scoring (TCCS) Method to Quanitfy Diverse Phenotypic Responses between Distinct Cell Types'*

To run the analysis use the `main.R` function, either within R or from the command line with `Rscript main.R`  which first checks and installs the necessary packages from CRAN and Github, and then runs the .R files in analysis/ to generate the plots used in the paper, which are saved as .eps files in figures/.

An R package containing the functions used in this analysis is available at [github.com/swarchal/tccs](http://www.github.com/swarchal/tccs).

**:warning: Warning:** Running the `load_packages()` or `main()` function will download and install current R packages from CRAN and GitHub if they are not already installed.

**Note:** The CellProfiler pipelines are provided without access to the raw images.
