# influenza-seasonality-shanghai
## System Requirements
- R version 4.3.0
- R packages (included in the code)
- GNU Scientific Library (GSL-2.7)

Before running the code, please unzip the following files, which are too large to be stored on GitHub without LFS:

- modeling_analysis/files/net_reproduction_number.csv.zip
- modeling_analysis/files/estimated_ili_pos.csv.zip
- modeling_analysis/files/calibrated_param_distributions.csv.zip

After that operation, you should have the same files without the `.zip` extension.

## Installation 
- To download and install R, visit [this website](https://cran.r-project.org/) and follow the "Download and Install R" instructions.
- To download and install GSL-2.7, follow the instructions described [here](https://www.gnu.org/software/gsl/)

- Set working directory to your root directory
- To run statistical analysis:
  - Go into the statistical_analysis folder and open the file statistical_analysis.R 
    - Install the R libraries needed for the analysis by setting INSTALL = T 
    - Run the analysis
- To run the modeling analysis:
  - Compile the code by opening a terminal in the main project directory and typing: *modeling_analysis/seasonal_influenza_model/MyCompiler.sh*
  - Go into the modeling_analysis folder and open the file modeling_analysis.R
    - Install the R libraries needed for the analysis by setting INSTALL = T 
    - Run the analysis in small sections separated by a header
    - To run the model, set RUN_SYSTEM = T


## Project directories
### Statistical analysis
Contains all the files needed to run the statistical analysis. 

**statistical_analysis.R** – Contains the code for running the descriptive analysis, the regression analysis, and estimating the predicted number of contacts. 

### Modeling analysis
Contains all the files needed to run the modeling analysis. 

**modeling_analysis.R** - Contains the code for running the modeling analysis and estimating the potential reproduction numbers, incidence of ILI<sup>+</sup> per 10,000, peak week, peak week incidence, and final infection attack rate.

**seasonal_influenza_model/Differential_main.c** - Contains the code for the influenza transmission model. 

**seasonal_influenza_model/myfn.c** - Contains the functions used for the influenza transmission model code.

**files/positive_ili.csv** - Contains estimates calculated using influenza-like-illness data collected by the CDC in Shanghai, China. 
