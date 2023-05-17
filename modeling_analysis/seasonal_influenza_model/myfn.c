#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "mcmc-io-main/mcmc_io.h"
#include "myfn.h"

/* Read parameters from file and saves them into a struct and 
calculates the other needed paramaters */
void inputParameters(Param *p, FILE *fp, const char* fname){
    int row = 1;
    char ch;
    const int Dec24 = 84;

    fp = fopen(fname, "r");
    if(fp == NULL){
        fprintf(stderr, "Error: cannot read Parameters file. Check file directory. %s\n", fname);
        fprintf(stderr, "File Path: %s", fname);
        exit(-1);
    }

    while(!feof(fp)){
        if(fgetc(fp) == '\n'){
            row++;
        }
        else{
            fscanf(fp, "%d, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %d, %d, %d, %d, %lf, %lf, %lf, %lf, %d, %lf, %lf, %lf, %d, %d, %lf, %d, %d, %d", 
            &p->N, &p->gt, &p->R0, &p->It0, &p->Rem0, &p->sigma_Beta, &p->sigma_I0, &p->sigma_RR, 
            &p->sigma_disp, &p->fluStartWeek, &p->fluStartYear, &p->fluEndWeek, &p->fluEndYear, &p->scalingFactor, 
            &p->RR, &p->RRholiday, &p->dispersion, &p->Nsim, &p->hh_cnt, &p->work_cnt, &p->school_cnt, &p->fluStartDay, &p->mcmc, &p->beta,
            &p->first_week, &p->start_holiday, &p->end_holiday); 
        }
    }
    //fprintf(stderr, "Rows = %d\n", row);

    if(p->sigma_Beta < 0 | p->sigma_disp < 0 | p->sigma_I0 < 0 | p->sigma_RR < 0){
        fprintf(stderr, "Error: Sigma is less than 0.\n");
        exit(0);
    }
    
    p->gamma = 1 / p->gt;
    if(p->gamma < 0 | p->gamma > 1){
        fprintf(stderr, "Error: Gamma is not between 0 and 1. %f\n", p->gt);
        exit(0);
    }

    if(p->beta <= 0){
        p->beta = (p->R0 * p->gamma) / (p->contacts[Dec24] + p->hh_cnt + p->work_cnt + p->school_cnt);
    }
    //fprintf(stdout, "Beta = %f", p->beta);
    //fprintf(stdout, "Contacts on Dec. 24: %f\n", p->contacts[Dec24]);

    p->St0 = p->N - p->It0;
}

/* Function that allocates memory for a character string and concatonates
the string into the output directory */
void outputFile(char *src, char **dest, char *fname){
    int x = 0;
    int y = 0;

    x = strlen(src);
    y = strlen(fname);

    *dest = (char *) malloc((x + y + 2 ) * sizeof(char));

    strcpy(*dest, src);
    strcat(*dest, fname);
    //fprintf(stdout, "%s\n", *dest);
}

/* Function that prints Paramters */
void printParams(Param *p){
    printf("Total population: %d\n", p->N); 
    printf("Generation time: %.2f\n", p->gt);
    printf("Basic Reproduction Number: %.2f\n", p->R0); 
    printf("Days: %d\n", p->t);
    printf("Gamma: %.3f\n", p->gamma); 
    printf("Beta: %.3f\n", p->beta);
    printf("S0: %.2f\n", p->St0); 
    printf("I0: %.2f\n", p->It0);
    printf("Rem0: %.2f\n", p->Rem0); 
    printf("Scaling Factor: %.2f\n", p->scalingFactor);
    printf("Beta Standard Deviation: %.2f\n", p->sigma_Beta); 
    printf("I0 Standard Deviation: %.2f\n", p->sigma_I0); 
    printf("Reporting Rate Standard Deviation: %.2f\n", p->sigma_RR); 
    printf("Flu season starts week %d of %d; Day: %d\n", p->fluStartWeek, p->fluStartYear, p->fluStartDay); 
    printf("Flu season ends week %d of %d\n", p->fluEndWeek, p->fluEndYear);
    printf("Size-negbinom parameter: %f\n", p->dispersion);
    printf("Total Simulations: %d\n", p->Nsim);
    printf("Contact Start: %d\n", p->contactStart);
    printf("Household Contacts: %f\n", p->hh_cnt);
    printf("Work Contacts: %f\n", p->work_cnt);
    printf("School Contacts: %f\n", p->school_cnt);
    printf("MCMC is %d\n", p->mcmc);
}

/* Allocates memory for pointers included in structs */ 
void memoryAllocation(Param *p, Output *r){
    r->cInc = (double *) calloc(p->t, sizeof(double));
    r->susceptible = (double *) calloc(p->t, sizeof(double));
    r->Rpotential = (double *) calloc(p->t, sizeof(double));
    r->Rnet = (double *) calloc(p->t, sizeof(double));
} 

/* Frees output pointer and sets it to NULL*/
void free_sys_output(Output *p){
    free(p->cInc); p->cInc = NULL;
    free(p->wklyInc); p->wklyInc = NULL;
    free(p->susceptible); p->susceptible = NULL;
    free(p->Rpotential); p->Rpotential = NULL;
    free(p->Rnet); p->Rnet = NULL;
    p->weeks = 0;
}

/* Estimates weekly Incidence from Cumulative Incidence */
void weeklyInc(Param *p, Output *r){
    int x = 1;
    r->weeks = 0; 

    // Calculating number of weeks 
    if(p->t % 7 == 0)
        r->weeks = p->t / 7; 
    else
        r->weeks = (p->t / 7) + 1; 

    //fprintf(stdout, "Number of Weeks: %d\n", r->weeks);
    //fprintf(stdout, "%d\n", p->t % 7);

    int *beginWeek = NULL;
    beginWeek = calloc((int) r->weeks, sizeof(int));

    int *endWeek =  NULL;
    endWeek = calloc((int) r->weeks, sizeof(int));

    //beginWeek[0] = 0;
    //endWeek[0] = 0;
    //fprintf(stderr, "%d %d %d\n", 0, beginWeek[0], endWeek[0]);

    for(int k = 0; k < r->weeks; k++){
        if(k == 0 && x > 0 | k == 0 && x < 0 ){
            fprintf(stderr, "Error: Weekly Incidence does not start at [0].");
            exit(0);
        }
        beginWeek[k] = x;

        if((x + 6) > p->t)
            endWeek[k] = p->t;
        else
            endWeek[k] = x + 6;

        x = x + 7;
        //fprintf(stderr, "%d %d %d\n", k, beginWeek[k], endWeek[k]);
    }

    r->wklyInc = (double *) calloc((int) r->weeks, sizeof(double));
    
    for(int i = 0; i < r->weeks; i++){
        for(int j = 0; j < 6; j++){
            r->wklyInc[i] = r->cInc[endWeek[i]] - r->cInc[beginWeek[i]];
        }

        if(r->wklyInc[i] < 0){
            r->wklyInc[i] = 0.; 
        }

        //fprintf(stdout, "Weekly Inc: %d %.5f %.5f %.5f\n", i, r->cInc[endWeek[i]], r->cInc[beginWeek[i]], r->wklyInc[i]);
    }

    free(beginWeek); 
    free(endWeek);
}

/* Aligns Contact Data Weeks */
void alignContactData(Param *p, ILIinput *d){
    const int YearWk = 52;

    if(p->fluStartYear == 2017 && p->fluStartWeek < p->first_week | p->fluEndYear == 2018 && p->fluEndWeek > p->first_week){
        fprintf(stderr, "Contact data does not align with ILI data. %d %d %d\n", p->fluStartWeek, p->fluEndWeek, p->first_week);
        exit (0);
    } 
    else if(p->fluStartYear == 2017 && p->fluStartWeek == p->first_week){
        p->contactStart = 0;
    }
    else{
        p->contactStart = ((p->fluStartWeek - (p->first_week)));
        
    }

    //fprintf(stdout, "%d\n", p->contactStart);
}

/* Function for opening a file to write and save data */
void writeWklyInc(FILE *fp, const char* fname, Param *p, Output *r, int sim, int n){
    int i;
    int j;

    fp =  fopen(fname, "w");
    if(fp == NULL){
        fprintf(stderr, "Error: cannot write Weekly Incidence file. Check file directory. %s\n", fname);
        exit(-1);
    }

    fprintf(fp, ",");
    for(i = 0; i < r->weeks; i++){
        fprintf(fp, "Week %d,", i + 1);
    }
    fprintf(fp, "\n");

    for(j = 0; j < sim; j++){
        fprintf(fp, "%d,", j);
        for(i = 0; i < r->weeks; i++){
            if(i >= (p->contactStart + p->start_holiday) && i <= (p->contactStart + p->end_holiday)){
                fprintf(fp, "%.4f,", r->wklyInc[i] * p->RRholiday);
                //fprintf(stdout, "%d, %.4f\n", i, r->wklyInc[i] * p->RRholiday);
            }
            else{
                fprintf(fp, "%.4f,", r->wklyInc[i] * p->RR);
                //fprintf(stdout, "%d, %.4f\n", i, r->wklyInc[i] * p->RR);
            }
        }
        fprintf(fp, "\n");
    }

}

void writeReproductionNumber(FILE *fp, Param *p, double *rep){
    int i;
    
    fprintf(fp, ",");
    for(i = 0; i < p->t; i++){
     fprintf(fp, "Time %d,", i + 1);
    }
    fprintf(fp, "\n");

    fprintf(fp, "0,");
    for(i = 0; i < p->t; i++){
        fprintf(fp, "%.4f,", rep[i]);
    }
    fprintf(fp, "\n");
}

/* Function for estimating differential SIR model and calculates the 
cumulative incidence using parameters beta, gamma, I0,, S0, Rem0, N and t*/
int func (double t, const double y[], double f[], void *p)
{
  (void)(t); // avoid unused parameter warning 
  int j = t; 
  Param *pParam = p;
  double transmission = pParam->beta;
  double population = pParam->N;
  double recovery = pParam->gamma;
  double scf = pParam->scalingFactor;
  double additionalContacts = pParam->hh_cnt + pParam->work_cnt + pParam->school_cnt;
  
  /* Print parameters passed into function
  fprintf(stdout, "N: %.2f\n", population);
  fprintf(stdout, "Recovery rate: %.3f\n", recovery);
  fprintf(stdout,"Per-contact transmission: %.2f\n", transmission);*/
  //fprintf(stdout,"Scaling Factor: %.2f\n", scf); 

    //fprintf(stdout,"S0: %f\n", y[0]);
    //fprintf(stdout,"I0: %f\n", y[1]);
    //fprintf(stdout,"Rem0: %f\n", y[2]);
    if(y[0] >= 0){
        f[0] = -((transmission * (pParam->contacts[j + pParam->fluStartDay] + additionalContacts)) * (y[1] / population) * y[0]); // 
        f[1] = ((transmission * (pParam->contacts[j+ pParam->fluStartDay] + additionalContacts)) * (y[1] / population) * y[0]) - (y[1] * recovery);
        f[2] = (y[1] * recovery);
        f[3] = ((transmission * (pParam->contacts[j+ pParam->fluStartDay] + additionalContacts)) * (y[1] / population) * y[0]);
        if(f[3] < 0){
            fprintf(stderr, "Error: Cumulative Incidence is less than 0.\n");
            exit(0);
        } 
    }
    
    // fprintf(stderr, "Cumulative Incidence at time %f: %f\n", t, f[3]);

    return GSL_SUCCESS;
}

/* Solves the system for the differential SIR model using function func (above)*/
void systemSolver(Param *p, Output *r){
    int i;
    double t = 0.0, t1 = 100;
    double y[4] = {p->St0, p->It0, p->Rem0, 0.};

    gsl_odeiv2_system sys = {func, NULL, 4, p};
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 0.001, 0.001, 0.0); 

    for(i = 0; i < p->t; i++){
        double ti = i * t1 / 100.0; 
        int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
          
        if (status != GSL_SUCCESS){
            fprintf (stderr, "Error, return value = %d\n", status);
            break;
            }

        r->susceptible[i] = (y[0] / p->N);
        if(r->susceptible[i] < 0 | r->susceptible[i] > 1){
          fprintf(stderr, "Proportion of susceptible is invalid.");
            exit(0);
        }
        r->cInc[i] = y[3];
       
        //fprintf (stdout, "%.0f %.5f %.5f %.5f %.5f %.5f\n", t, y[0], y[1], y[2], y[0] + y[1] + y[2], y[3]);
        // fprintf(stdout, "%.5f %.5f\n", pResults->cInc[i], y[3]);
        // fprintf(stdout, "Susceptible Proportion: %f\n", r->susceptible[i]);
    } 

   p->attackRate = r->cInc[p->t - 1] / p->N;
   //fprintf(stdout, "AR: %f\n", p->attackRate);

  gsl_odeiv2_driver_reset(d); 
  gsl_odeiv2_driver_free (d); 
}

/* Calculates the Potential and Net Reproduction Numbers */
void ReproductionNumber(Param *p, Output *r){
    int i; 
    double additionalContacts = p->hh_cnt + p->work_cnt + p->school_cnt;

    for(i = 0; i < p->t; i++){
        r->Rpotential[i] = (p->beta * (p->contacts[i + p->fluStartDay] + additionalContacts)) / p->gamma;
        r->Rnet[i] = r->Rpotential[i] * r->susceptible[i];
        //fprintf(stdout, " %d, Susceptible: %f, Potential: %f, Net: %f\n", i, r->susceptible[i], r->Rpotential[i], r->Rnet[i]);

        if(r->Rpotential < 0 | r->Rnet < 0){
            fprintf(stderr, "Error: Reproduction numbers are less than 0.\n");
            exit(0);
        }
        if(r->Rpotential[i] < r->Rnet[i]){
            fprintf(stderr, "Error: Net reproduction number is greater than potential reproduction number.\n");
            exit(0);
        }
    }
}

/* Estimates the log likelihood for a poisson distribution from
estimated weekly incidence from the SIR model and estimated weekly
incidence produced from ILI data */
double poisson_logLL(Param *p, Output *r, int *ili, int n){ // , Param *p
    double loglike = 0.;

    for(int i = 0; i <= n; i++){
        //fprintf(stdout, "%d %f %d\n", i, r->wklyInc[i], ili[i]);
        // fflush(stdout);
        if(r->wklyInc[i] > 0 && ili[i] > 0){
            if(i >= p->contactStart && i <= n){
                loglike += log(gsl_ran_poisson_pdf(ili[i], r->wklyInc[i] * p->RRholiday)); 
                fprintf(stdout, "Week: %d ILI: %d Model: %f PDF: %e\n", 
                i, ili[i], r->wklyInc[i] * p->RRholiday, gsl_ran_poisson_pdf(ili[i], r->wklyInc[i] * p->RRholiday));
            }
            else{
                loglike += log(gsl_ran_poisson_pdf(ili[i], r->wklyInc[i] * p->RR)); 
                fprintf(stdout, "Week: %d ILI: %d Model: %f PDF: %e\n", 
                i, ili[i], r->wklyInc[i] * p->RR, gsl_ran_poisson_pdf(ili[i], r->wklyInc[i] * p->RR));
            }
            // fflush(stdout);
        }
    }
    // printf("Log likelihood: %f\n", loglike);

    return loglike;
} 

long double poisson_logLL_long(Param *p, Output *r, int *ili, int n){ // , Param *p
    long double loglike = 0.;

    for(int i = 0; i <= n; i++){
        //fprintf(stdout, "%d %f %d\n", i, r->wklyInc[i], ili[i]);
        // fflush(stdout);
        if(r->wklyInc[i] > 0 && ili[i] > 0){
            if(i >= p->contactStart && i <= n){
                loglike += log(gsl_ran_poisson_pdf(ili[i], r->wklyInc[i] * p->RRholiday)); 
                fprintf(stdout, "Week: %d ILI: %d Model: %f PDF: %e\n", 
                i, ili[i], r->wklyInc[i] * p->RRholiday, gsl_ran_poisson_pdf(ili[i], r->wklyInc[i] * p->RRholiday));
            }
            else{
                loglike += log(gsl_ran_poisson_pdf(ili[i], r->wklyInc[i] * p->RR)); 
                fprintf(stdout, "Week: %d ILI: %d Model: %f PDF: %e\n", 
                i, ili[i], r->wklyInc[i] * p->RR, gsl_ran_poisson_pdf(ili[i], r->wklyInc[i] * p->RR));
            }
            // fflush(stdout);
        }
    }
    // printf("Log likelihood: %f\n", loglike);

    return loglike;
} 

/* Estimates the log likelihood for a negative binomial distribution from
estimated weekly incidence from the SIR model and estimated weekly
incidence produced from ILI data */

double nbinom_logLL(Param *p, Output *r, int *ili, int n){ 
    double loglike = 0.;
    double* a; 
    double *probability = NULL;
    bool* holiday = 0;

    probability = (double *) calloc(r->weeks, sizeof(double));
    a = (double *) calloc(r->weeks, sizeof(double));
    holiday = (bool *) calloc(r->weeks, sizeof(bool));

    for(int i = 0; i < r->weeks; i++){
        if(i >= (p->contactStart + p->start_holiday) && i <= (p->contactStart + p->end_holiday)){
            // if(i >= p->contactStart && i <= n){
            a[i] = r->wklyInc[i]*p->RRholiday; 
            
            holiday[i] = 1;
            //fprintf(stdout, "%d, %f, %f, %d\n", i, r->wklyInc[i], a, holiday);
        } else{
            a[i] = r->wklyInc[i]*p->RR; 

            holiday[i] = 0;
            //fprintf(stdout, "%d, %f, %f, %d\n",i, r->wklyInc[i], a, holiday);
        }
        //fprintf(stdout, "%d %d %d %d\n", i, p->contactStart, n, ili[i]);
        //fprintf(stdout, "%f\n", r->wklyInc[i]);

        probability[i] = (p->dispersion) /(a[i] + p->dispersion);
        //fprintf(stdout, "%d, INC: %f, p: %f size: %f\n", i, a, probability[i], p->dispersion);
        if(probability[i] < 0 | probability[i] > 1){
            fprintf(stderr, "Probability is not between 0 and 1.");
            exit(-1);
        }
    }

    //fprintf(stdout, "%d\n", p->contactStart);
    for(int i = 0; i < n; i++){
        //fprintf(stdout, "%d, %d, %f, %f, %d\n", i, ili[i], r->wklyInc[i], a[i], holiday[i]);

        if(r->wklyInc[i + p->contactStart] > 0 && ili[i] > 0){
            loglike += log(gsl_ran_negative_binomial_pdf(ili[i], probability[i + p->contactStart], p->dispersion)); 
            /*fprintf(stdout, "Week: %d ILI: %d Model: %f PDF: %e\n", 
              i, ili[i], r->wklyInc[i] * p->RR, gsl_ran_poisson_pdf(ili[i], r->wklyInc[i] * p->RR));
             fflush(stdout);*/
        }
    }
    //fprintf(stdout, "Log likelihood: %f\n", loglike);

    free(probability);
    free(a);

    return loglike;
} 

/*void nbinom_sizeParameter(FILE *fp, const char *fname, Param *p){
    p->lengthOfSize = 1;
    char buffer[1024];
    char *tok;
    int i;

    fp = fopen(fname, "r");
    if(fp == NULL){
        fprintf(stderr, "Error: cannot read file. Check file directory. %s \n", fname);
        exit(-1);
    }

    while(!feof(fp)){
        if(fgetc(fp) == '\n'){
             p->lengthOfSize++;
        }
    } // fprintf(stdout, "Row: %d\n", p->lengthOfSize);

    p->sizeIndex = (int *) calloc(p->lengthOfSize, sizeof(int));
    p->size = (double *) calloc(p->lengthOfSize, sizeof(double));

    fgets(buffer, sizeof(buffer), fp);
    while(!feof(fp)){ 
        tok = strtok(buffer, ",");
        while(tok != NULL){
            p->sizeIndex[i] = strtock(NULL, ",");
        }
        
    }
   


    fclose(fp);
    fp = NULL;

    // Ignore first row - loop over rows
    // Ignore first column 
    // Save to Param pointer
}*/
