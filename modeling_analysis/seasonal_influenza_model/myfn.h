#include "mystruct.h"
#include "mcmc-io-main/mcmc_io.h"

#ifndef MYFN_H_
#define MYFN_H_

void inputParameters(Param *p, FILE *fp, const char* fname);
void printParams(Param *p);
void outputFile(char *src, char **dest, char *fname);
void nbinom_sizeParameter(FILE *fp, const char *fname, Param *p);
void writeWklyInc(FILE *fp, const char* fname, Param *p, Output *r, int sim, int n);
void writeReproductionNumber(FILE *fp, Param *p, double *rep);
int func(double t, const double y[], double f[], void *params);
void memoryAllocation(Param *p, Output *r);
void free_sys_output(Output *p);
void systemSolver(Param *p, Output *r); 
void weeklyInc(Param *p, Output *r);
double poisson_logLL(Param *p, Output *r, int *i, int n);
long double poisson_logLL_long(Param *p, Output *r, int *ili, int n);
double nbinom_logLL(Param *p, Output *r, int *ili, int n);
void alignContactData(Param *p, ILIinput *d);
void ReproductionNumber(Param *p, Output *r);

#endif