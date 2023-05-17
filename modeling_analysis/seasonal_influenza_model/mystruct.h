#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#ifndef MYSTRUCT_H_
#define MYSTRUCT_H_

// Struct that stores parameters needed for differential SIR model
typedef struct{
    double gt; // Generation time
    int N; // Total population
    int t; // Size of contacts
    double gamma; // Recovery rate
    double R0; // Basic Reproduction Number
    double beta; // Per-contact transmission rate
    double St0; // Susceptible at time 0
    double It0; // Infectious at time 0
    double Rem0; // Removed at time 0
    double RR; // Reporting rate during the regular period
    double RRholiday; // Reporting rate during the holiday season 
    double scalingFactor; // Scaling factor for final output
    double sigma_Beta; // standard deviation for Beta
    double sigma_I0; // standard deviation for initial Infectious
    double sigma_RR; // standard deviation for reporting rate
    int fluStartDay; // Day flu season starts
    int fluStartWeek; // Week of the year flu season begins
    int fluStartYear; // Year of interest for flu season beginning
    int fluEndWeek; // Week of the year flu season ends
    int fluEndYear; // Year of interest for flu season beginning
    double dispersion; // size parameter for negative binomial distribution
    int Nsim; // Number of simulations
    double sigma_disp; // standard deviation for dispersion
    double *contacts; // Pointer to store data for mean number of contacts for each week
    int contactStart; // First day of contact data used in model
    double attackRate; // Infection Attack Rate
    double hh_cnt; // Household Contacts
    double work_cnt; // Work Contacts
    double school_cnt; // School Contacts
    int mcmc; // MCMC running. 1 = True, 0 = False
    int first_week; // First week to run
    int start_holiday; // Number of weeks from first_week until the holiday starts
    int end_holiday; // Number of weeks from first_week until the holiday ends
} Param;

// Struct that stores output from differential SIR model
typedef struct { 
    double *cInc; // Cumulative Incidence
    double *wklyInc; // Weekly Incidence - Use as expected
    double *Rpotential; // Potential Reprodution Number
    double *Rnet; // Net Reproduction Number
    double *susceptible; // Proportion of susceptible individuals at time t
    int weeks; // Size of weekly incidence vector
} Output;

#endif