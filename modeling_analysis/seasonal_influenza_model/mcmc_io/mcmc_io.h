#ifndef MCMC_IO_H
#define MCMC_IO_H

#define FILE_BUF_SIZE 1024  // Size, in bytes, of the chunks of the file that are read at each input operation.

// Struct that stores ILI data read from file.
typedef struct {
    size_t size;  // Number of elements in each array.
    int *year;    // Year data was collected
    int *week;    // Week of the year data was collected
    int *estInc;  // Estimated incidence for H1pdm
} ILIinput;

void free_ili_input(ILIinput* data_p);
int read_ili_csv(const char* fname, ILIinput* data_p);

#endif
