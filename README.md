# influenza-seasonality-shanghai
## System Requirements
- R version 4.3.0
- R packages (included in the code)
- GNU Scientific Library (GSL-2.7)
- C compiler (gcc or clang)

Before running the code, please unzip the following files, which are too large to be stored on GitHub without LFS:

- modeling_analysis/files/net_reproduction_number.csv.zip
- modeling_analysis/files/estimated_ili_pos.csv.zip
- modeling_analysis/files/calibrated_param_distributions.csv.zip

After that operation, you should have the same files without the `.zip` extension.

## Installation 
- To download and install R, visit [this website](https://cran.r-project.org/) and follow the "Download and Install R" instructions. This will take only a few minutes. 
- To download and install GSL-2.7, follow the instructions described [here](https://www.gnu.org/software/gsl/). This may take around 10 minutes to install. 

## Instructions for Use
- Set working directory to the project root
- To run statistical analysis:
  - Go into the statistical_analysis folder and open the file statistical_analysis.R 
    - Install the R libraries needed for the analysis by setting INSTALL = T 
    - Run the analysis by highlighting small sections of the code and clicking the "run" command. 
    - The full analysis will take around 10 minutes to run. 
    - To make sure the code is working, check that the plots produced from the code are similar to those in the main text. 
- To run the modeling analysis:
  - Compile the code by opening a terminal in the main project directory and typing: *modeling_analysis/seasonal_influenza_model/MyCompiler.sh*
    - The ./Mycompiler.sh script will only work when running from the seasonal_influenza_model folder.
  - Go into the modeling_analysis folder and open the file modeling_analysis.R
    - Install the R libraries needed for the analysis by setting INSTALL = T 
    - Run the analysis in small sections by highlighting the code you want to run and then clicking the "run" command. Each section of the code is separated by a header.
    - To run the influenza transmission model, set RUN_SYSTEM = T
    - Warning: the full analysis may take a few hours to run with the default setup. This time can be reduced by reducing the number of simulations or commenting the sections that estimate the reproduction number (lines 86-114 and lines 276-331).
    - - To make sure the code is working, check that the plots produced from the code are similar to those in the main text. 

## Project directories
### Statistical analysis
Contains all the files needed to run the statistical analysis. 

**statistical_analysis.R** – Contains the code for running the descriptive analysis, the regression analysis, and estimating the predicted number of contacts. 

***Expected Estimates:***
- Descriptive statistics of the variables included in the analysis. (Table 1)
- Regression coefficients, p-values, and confidence intervals of all variables included in the regression analysis. (Table 1)
- Average daily number of contacts by week from October 1, 2017, to September 30, 2018, as predicted using the regression results. 

***Expected Outputs:***
- Figure 1A -  Daily maximum temperature for each day from October 1, 2017, to September 30, 2018, the 505 seasonal trend of the temperatures, daily variation between the maximum temperature, and seasonal trend.
- Figure 1B - Estimated daily number of total contacts for each week from October 1, 2017, to September 30, 2018.

### Modeling analysis
Contains all the files needed to run the modeling analysis. 

**modeling_analysis.R** - Contains the code for running the modeling analysis and estimating the potential reproduction numbers, incidence of ILI<sup>+</sup> per 10,000, peak week, peak week incidence, and final infection attack rate.

**seasonal_influenza_model/Differential_main.c** - Contains the code for the influenza transmission model. 

**seasonal_influenza_model/myfn.c** - Contains the functions used for the influenza transmission model code.

**files/positive_ili.csv** - Contains estimates calculated using influenza-like-illness data collected by the CDC in Shanghai, China. 

***Expected Estimates:***
- Potential reproduction number for each week from October 1, 2017, to September 30, 2018.
- Incidence of ILI<sup>+</sup> per 10,000 
- Peak week 
- Peak week incidence
- Final infection attack rate.

***Expected Output:***
- Figure 2A - Estimated weekly incidence of ILI+ infections during the 2017-2018 influenza season in Shanghai, China.
- Figure 2B - Estimated posterior distribution of weekly net reproduction number.
- Figure 2C - Posterior distribution of the final infection attack rate for the 2017-2018 influenza season.
- Figure 3A - Estimated weekly contact-dependent basic reproduction number from October 1, 2017, to September 30, 2018, using the estimated posterior distribution of the transmission risk (𝛽) for the 2017-2018 influenza season.
- Figure 3B - Infection attack rate, peak week incidence of ILI+, and peak week for one year of simulations for different epidemic starting times.

