# influenza-seasonality-shanghai

**If you do not have R installed on your computer, visit [this website](https://cran.r-project.org/) and follow the "Download and Install R" instructions.**

## Instructions
- Open statistical_analysis.R 
  - Check that the working directory contains the following: functions.R, XX, temperature_shanghai_2017_2018.csv, and regression_output.csv
  - Check that the source file is functions.R
  - Install the R libraries needed for the analysis by setting INSTALL = T 
  - Run the analysis
- Download and install the GNU Scientific Library by following the instructions described [here](https://www.gnu.org/software/gsl/)
- Compile the code by opening a terminal in the main project directory and typing: XX
- Open model_script.R
- Check that the working directory contains the following: XX.csv and XX

## Project directories
### Statistical analysis
Contains all the files needed to run the statistical analysis. 

**Files:**

**statistical_analysis.R** – Contains the code for running the descriptive analysis, the regression analysis, and estimating the predicted number of contacts. Requires input from: functions.R, temperature_shanghai_2017_2018.csv, XX.csv, and regression_output.csv.

**functions.R** - Contains functions needed for the statistical analysis that are not included in an R library. Serves as a source file for statistical_analysis.R.

**temperature_shanghai_2017_2018.csv** - Contains the maximum daily temperature data in Shanghai, China from October 1, 2017, to September 30, 2018. This data was obtained from wunderground.com. 

**XX.csv** - Contains a cleaned and modified subset dataset of the diary-based contact survey data that includes the variables of interest in the current study. The complete contact survey dataset can be found on Zenodo - include link here / reference. 

### Modeling analysis
Contains all the files needed to run the statistical analysis. 

**Files:**

**model_script.R** - Contains the code for running the modeling analysis and estimating the incidence of ILI<sup>+</sup> per 10,000, peak week, peak week incidence, and final infection attack rate. Produces Figures 2 and 3. Requires input from: estimated_contacts_shanghai.csv, number_of_ili_.csv, mcmc_results.csv, and model_output.R.

**estimated_contacts_shanghai.csv** - Contains the estimated mean number of contacts produced from the statistical analysis. 

**number_of_ili_.csv** - Contains the number of observed ILI<sup>+</sup>. REF.

**mcmc_results.csv** - Contains the posterior distributions for the initial number of infected individuals, the per-contact transmission risk, the reporting rate, and the over-dispersion parameter.

**model_output.R** - Contains the resulting transmission dynamics produced from the model in a format that can easily be read for estimating the incidence of ILI<sup>+</sup> per 10,000, peak week, peak week incidence, and final infection attack rate. 




- Data license stuff:
  - Clean data (The CC 4.0 Attrib Zenodo)
  - Wunderground:
  - CDC China
