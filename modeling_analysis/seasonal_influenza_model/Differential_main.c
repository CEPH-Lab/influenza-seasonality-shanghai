#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "myfn.h"
#include "mcmc-io-main/mcmc_io.h"

// Input files
#define PARAMFILE "Parameters_2.csv" // "Parameters_3.csv" // "Parameters_4.csv" //   "Parameters_1.csv" // 
#define ILIFILE "ILIcombineddataH1pdm.csv" // "ILIcombineddataBY.csv" //
#define CONTACTS "contactsHolidayExtended.csv" // "contactsHoliday.csv" // "contacts40reduced.csv" // "contacts.csv" //    

// Output files
#define PARAMOUTPUT "Param_Output.csv" // Outuput file for parameters
#define WKLYINC "EstWklyIncidence.csv" // Output - Est. Inc
#define POT_REP "Potential_Reproduction.csv" // Output - Potential Rep
#define NET_REP "Net_Reproduction.csv" // Output - Net Rep

int main(int argc, char* argv[]){ 
  char *paramfile = NULL;  
  char *contacts = NULL;

  char *output = "";
  char *parameters = NULL;
  char *wklyinc = NULL;
  char *cumulative = NULL;
  char *pot_rep = NULL;
  char *net_rep = NULL;

  switch(argc){

    case 1:
      paramfile = PARAMFILE;
      contacts = CONTACTS; 
      break;
    
    case 2:
      paramfile = argv[1];
      contacts = CONTACTS; 
      break;
    
    case 3:
      paramfile = argv[1];
      contacts = argv[2];
      break;
    
    default:
      paramfile = argv[1];
      contacts = argv[2];
      output = argv[3];
      break;
  }
  //printf("%s\n", paramfile);
  //printf("%s\n", contacts);

  outputFile(output, &parameters, PARAMOUTPUT);
  outputFile(output, &wklyinc, WKLYINC);
  outputFile(output, &cumulative, "CumulativeIncidence.csv");
  outputFile(output, &pot_rep, POT_REP);
  outputFile(output, &net_rep, NET_REP);

  char *data = NULL;
  data = ILIFILE;

  FILE *fParameters = NULL; // File pointer for Parameter Input
  FILE *fWklyInc = NULL; // File pointer to Weekly Incidence Output file
  FILE *fCumulative = NULL; // File pointer to Weekly Incidence Output file
  FILE *fParamOutput = NULL; // File pointer to Parameter Output file
  FILE *fContacts = NULL; // File pointer to contact data
  FILE *fPotentialRep = NULL; // File pointer to potential reproduction output file
  FILE *fNetRep = NULL; // File pointer to net reproduction output file

 Param param = {0., 0, 0, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 , 0, 0, 0, 0, 0, 0}; 
  ILIinput ili;
  Output results, *pResults = NULL;
  pResults = &results; 
  
  double loglike; // Accepted Log Likelihood
  
  int *fluSeasonBegins = NULL;
  int *fluSeasonEnds = NULL;
  bool accept; // Boolean for Accepting and Rejecting Likelihoods. Initialized later.
  int count = 0; // Number of weeks during flu season
  int i;

  // Set seed for generator
  gsl_rng *gsl;
  gsl_rng_env_setup();
  gsl = gsl_rng_alloc (gsl_rng_default); // fprintf(stdout, "Generator name: %s\n", gsl_rng_name(gsl)); 


  /* Temporary variables for storing "candidate" values */
  Param param_tmp = {0., 0, 0, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 , 0, 0, 0, 0, 0, 0}; 
  Output results_tmp, *pResults_tmp = NULL;
  pResults_tmp = &results_tmp;
  double loglike_tmp; // "Candidate" Log Likelihood

  /* Reading Data and Parameters into program */
  read_ili_csv(data, &ili); // Reads in ILI data
  if(read_ili_csv(data, &ili)){
    free_ili_input(&ili);
    return EXIT_FAILURE;
  }
  /*for (int i = 0; i < ili.size; i++){
    printf("Year: %d, Week: %d, Weekly Incidence: %d\n", ili.year[i], ili.week[i], ili.estInc[i]);
  }
  printf("Data has %lu entries.\n", ili.size); */

  read_csv_double_vector(contacts, &(param.contacts), &(param.t));
  if(read_csv_double_vector(contacts, &(param.contacts), &(param.t))){
    free(&(param.contacts));
    free(&(param.t));
    return EXIT_FAILURE;
  }
  /*for(int i = 0; i < param.t; i++){
    fprintf(stdout, "%d, %f\n", i, param.contacts[i]);
  }*/

  read_csv_double_vector(contacts, &(param_tmp.contacts), &(param_tmp.t));
  if(read_csv_double_vector(contacts, &(param_tmp.contacts), &(param_tmp.t))){
    free(&(param_tmp.contacts));
    free(&(param_tmp.t));
    return EXIT_FAILURE;
  }

  /* Input the Parameters from file and Align it with the Contact Data */
  inputParameters(&param, fParameters, paramfile ); inputParameters(&param_tmp, fParameters, paramfile ); 
  alignContactData(&param, &ili); alignContactData(&param_tmp, &ili); 
  //printParams(&param); printParams(&param_tmp);

  /* Allocate memory for Results and Temporary Results */
  memoryAllocation(&param, pResults);
  memoryAllocation(&param_tmp, pResults_tmp);

  /* Identifying Beginning and End of Flu Season */
  for(i = 0; i <= ili.size; i++){
    if(ili.year[i] == param.fluStartYear && ili.week[i] == param.fluStartWeek){
      fluSeasonBegins = &ili.estInc[i];
      // fprintf(stderr, "Flu Start: week %d of %d.\n Incidence = %d\n", ili.week[i], ili.year[i], ili.estInc[i]);
    }
    if(ili.year[i] == param.fluEndYear && ili.week[i] == param.fluEndWeek){
      fluSeasonEnds = 1 + &ili.estInc[i];
      // fprintf(stderr, "Flu End: week %d of %d.\n Incidence = %d\n", ili.week[i], ili.year[i], ili.estInc[i]);
    }
  }

  count = fluSeasonEnds - fluSeasonBegins;  
  if(count <= 0){
    fprintf(stderr, "Error: Not enough weeks.\n");
    exit(0);
  }
  //fprintf(stdout, "Weeks = %d\n", count);

  /* Estimates Initial Results with Initial Parameters */
  systemSolver(&param, pResults);
  weeklyInc(&param, pResults); // fprintf(stdout, "Flu: %d Model: %d\n", count, pResults->weeks);
  if(count > pResults->weeks){
    fprintf(stderr, "Error: Number of weeks in the data exceeds those in the model.\n");
    exit(0);
  }
  /*for(i = 0; i < pResults->weeks; i++){
    fprintf(stdout, "Estimated Weekly Incidence: %d %.5f\n", i, pResults->wklyInc[i]);
    //fflush(stdout);
  }*/
  //fprintf(stdout, "AR: %f\n", param.attackRate);
  ReproductionNumber(&param, pResults);
  
  writeWklyInc(fWklyInc, wklyinc, &param, pResults, 1, count);

  fCumulative =  fopen(cumulative, "w+");
  if(fCumulative == NULL){
      fprintf(stderr, "Error: cannot write Weekly Incidence file. Check file directory. %s\n", cumulative);
      exit(-1);
  }
  writeReproductionNumber(fCumulative, &param, pResults->cInc);
  
  fPotentialRep =  fopen(pot_rep, "w");
  if(fPotentialRep == NULL){
    fprintf(stderr, "Error: cannot write Potential Reproduction Number file. Check file directory. %s\n", pot_rep);
    exit(-1);
  }
  writeReproductionNumber(fPotentialRep, &param, pResults->Rpotential);

  fNetRep =  fopen(net_rep, "w");
  if(fNetRep == NULL){
    fprintf(stderr, "Error: cannot read Net Reproduction Number file. Check file directory. %s\n", net_rep);
    exit(-1);
  }
  writeReproductionNumber(fNetRep, &param, pResults->Rnet);
  

  /**************************/
  /*          MCMC          */
  /**************************/

  if(param.mcmc == 1){
    //fprintf(stdout, "MCMC is running...\n");

    loglike = nbinom_logLL(&param, pResults, fluSeasonBegins, count);
    if(loglike == -INFINITY){
      fprintf(stderr, "Log Likelihood is -inf.\n");
      exit(0);
    }

    /* Open Files to Write Output */
    fParamOutput =  fopen(parameters, "w");
    if(fParamOutput == NULL){
      fprintf(stderr, "Error: cannot write Parameter Output file. Check file directory. %s\n", parameters);
      exit(-1);
    }
    fprintf(fParamOutput, "Simulation, AR, Beta_0, Beta_1, I0_0, I0_1, RR_0, RR_1, RRholiday_0, RRholiday_1, Size_0, Size_1, LL_0, LL_1, Accept\n");

    /* Estimated Weekly Incidence */
    fWklyInc =  fopen(wklyinc, "w");
    if(fWklyInc == NULL){
      fprintf(stderr, "Error: cannot write Weekly Incidence File. Check file directory. %s\n", wklyinc);
      exit(-1);
    }

    fprintf(fWklyInc, ",");
    for(i = 0; i < pResults->weeks; i++){
      fprintf(fWklyInc, "Week %d,", i + 1);
    }
    fprintf(fWklyInc, "\n");

    for(int j = 0; j < 1; j++){
      fprintf(fWklyInc, "%d,", j);
      for(i = 0; i < pResults->weeks; i++){
        if(i >= (param.contactStart + param.start_holiday) && i <= (param.contactStart + param.end_holiday)){
          fprintf(fWklyInc, "%.4f,", pResults->wklyInc[i] * param.RRholiday);
        }
        else{
          fprintf(fWklyInc, "%.4f,", pResults->wklyInc[i] * param.RR);
        }
      }
      fprintf(fWklyInc, "\n");
    }

    /* MCMC Iterations for Log Likelihood Comparisons */
    for(int j = 0; j < param.Nsim; j++){

      accept = false;
      param_tmp.beta = param.beta + gsl_ran_gaussian(gsl, param.sigma_Beta); 
      param_tmp.It0 = param.It0 + gsl_ran_gaussian(gsl, param.sigma_I0); 
      param_tmp.RR = param.RR + gsl_ran_gaussian(gsl, param.sigma_RR); 
      param_tmp.dispersion = param.dispersion + gsl_ran_gaussian(gsl, param.sigma_disp); 
      param_tmp.RRholiday = param.RRholiday + gsl_ran_gaussian(gsl, param.sigma_RR); // param_tmp.RR; // 
      // printParams(&param_tmp);
    
      if(param_tmp.beta < 0 | param_tmp.It0 < 0 | param_tmp.RR < 0 | param_tmp.RR > 1 | param_tmp.dispersion < 0 
      |  param_tmp.RRholiday < 0 | param_tmp.RRholiday > 1){
        // | param_tmp.beta < param_tmp.gamma){ 
        //n Add maximum contacts
        goto PRINT;
      }
    
      systemSolver(&param_tmp, pResults_tmp);
      weeklyInc(&param_tmp, pResults_tmp);
      /*printf("%d\n", j);
      for(i = 0; i < pResults_tmp->weeks; i++){
        fprintf(stdout, "Weekly Inc: %d %e\n", i , pResults_tmp->wklyInc[i]);
        fflush(stdout);
      }*/
      ReproductionNumber(&param_tmp, pResults_tmp);

      loglike_tmp = nbinom_logLL(&param_tmp, pResults_tmp, fluSeasonBegins, count);

      //printf(stdout, "%f\n", loglike_tmp); fflush(stdout);

      if(gsl_ran_flat(gsl, 0, 1) < exp(loglike_tmp - loglike)){
        accept = true;
	    } else{
        accept = false;
      }

      PRINT: fprintf(fParamOutput, "%d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %d\n", 
      j, param.attackRate, param.beta, param_tmp.beta, param.It0, param_tmp.It0, param.RR, param_tmp.RR, 
      param.RRholiday, param_tmp.RRholiday, param.dispersion, param_tmp.dispersion, loglike, loglike_tmp, accept);

      if(accept == true){
        loglike = loglike_tmp;
	      param.beta = param_tmp.beta;
        param.RR = param_tmp.RR;
        param.dispersion = param_tmp.dispersion;
        param.It0 = param_tmp.It0;
        param.attackRate = param_tmp.attackRate;
        param.RRholiday = param_tmp.RRholiday;
        for(i = 0; i < pResults_tmp->weeks; i++){
          pResults->wklyInc[i] = pResults_tmp->wklyInc[i];
        }
        for(i = 0; i < param.t; i++){
          pResults->Rpotential[i] = pResults_tmp->Rpotential[i];
        }
        for(i = 0; i < param.t; i++){
          pResults->Rnet[i] = pResults_tmp->Rnet[i];
        }
	    }

      /* Prints Weekly Incidence */
      fprintf(fWklyInc, "%d", j);
      for(i = 0; i < pResults->weeks; i++){
       if(i >= (param.contactStart + param.start_holiday) && i <= (param.contactStart + param.end_holiday)){
          fprintf(fWklyInc, ",%.4f", pResults->wklyInc[i] * param.RRholiday);
       }
        else{
         fprintf(fWklyInc, ",%.4f", pResults->wklyInc[i] * param.RR);
       }
      }
     fprintf(fWklyInc, "\n");


     fprintf(fCumulative, "0,");
      for(int k = 0; k < param.t; k++){
        fprintf(fCumulative, "%.4f,", pResults->cInc[k]);
      }
      fprintf(fCumulative, "\n");

      /* Prints Daily Potential Reproduction Number */
      fprintf(fPotentialRep, "%d", j);
      for(i = 0; i < param.t; i++){
       fprintf(fPotentialRep, ",%.4f", pResults->Rpotential[i]);
      }
      fprintf(fPotentialRep, "\n");

      /* Prints Daily Net Reproduction Number */
      fprintf(fNetRep, "%d", j);
      for(i = 0; i < param.t; i++){
        fprintf(fNetRep, ",%.4f", pResults->Rnet[i]);
      }
      fprintf(fNetRep, "\n"); 

      fprintf(stdout, "MCMC simulation %d is finished.\n", j);
    }

    fprintf(stdout, "MCMC is finished.");
  }

  else {
    fprintf(stdout, "MCMC did not run.\n");
  }
  
  fclose(fParameters); fParameters = NULL;
  fclose(fParamOutput); fParamOutput = NULL;
  fclose(fWklyInc); fWklyInc = NULL;
    
  // Free memory
  free_ili_input(&ili);
  free_sys_output(pResults);
  // free_sys_output(pResults_tmp);
  free(parameters);
  free(wklyinc);
  free(cumulative);
  free(pot_rep);
  free(net_rep);

  return 0;
}
