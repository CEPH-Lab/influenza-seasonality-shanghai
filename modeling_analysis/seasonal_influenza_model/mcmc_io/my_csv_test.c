/*
Read and store ILI data from a csv file.

User manual of libcsv.
https://github.com/rgamble/libcsv/blob/master/csv.pdf
*/

#include <stdio.h>
#include <stdlib.h>
#include "mcmc_io.h"


int main (int argc, char *argv[])
{
  char* fname;
  ILIinput data;

  // Check the arguments given to the program call
  if (argc < 2) {
    fprintf(stderr, "Please inform the csv file name as an arugment.\n");
    exit(EXIT_FAILURE);
  }
  fname = argv[1];

  // Read the file
  read_ili_csv(fname, &data);

  // Print its content
  for (int i = 0; i < data.size; i++){
    printf("%d, %d, %d\n", data.year[i], data.week[i], data.estInc[i]);
  }
  printf("Data has %lu entries.\n", data.size);

  // Free data with interface function.
  free_ili_input(&data);

  exit(EXIT_SUCCESS);
}
 
