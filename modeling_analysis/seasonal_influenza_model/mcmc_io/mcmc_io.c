/* 
Toolset for data input/output (io) for the Influenza MCMC project.

v0.01 – Development version. Performs the basic reading without complete error checking.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "/Users/alkummer/C_PROJECTS/Exercises/V2/libcsv/csv.h"

#include "mcmc_io.h"


// ------------------------------------------------------------------------------------------------
// AUXILIARY STRUCTS AND FUNCTIONS
// ------------------------------------------------------------------------------------------------

// Auxiliary struct with extra variables to help on the file parsing.
typedef struct {
    ILIinput* data_p;

    size_t capacity;  // Assured number of allocated positions in each vector.
    size_t curr_row;  // Row index (1-based) currently being read.
    size_t curr_col;  // Column index (1-based) currently being read.

    int err_status;
    size_t err_row;
    size_t err_col;

} ILIinputAux;

/*
Frees the dynamically allocated arrays for the ILIinput struct. 
Sets its pointers to NULL and its size to 0.
*/
void free_ili_input(ILIinput* data_p){
    free(data_p->year); data_p->year = NULL;
    free(data_p->week); data_p->week = NULL;
    free(data_p->estInc); data_p->estInc = NULL;
    data_p->size = 0;
}


static int is_space(unsigned char c) {
  if (c == CSV_SPACE || c == CSV_TAB) return 1;
  return 0;
}

static int is_term(unsigned char c) {
  if (c == CSV_CR || c == CSV_LF) return 1;
  return 0;
}


// ------------------------------------------------------------------------------------------------
// PARSING CALLBACK FUNCTIONS
// ------------------------------------------------------------------------------------------------

/* 
Callback function for each field that is read from file.
*/
void cb1 (void *s_v, size_t len, void *aux_vp) { 

  // Convert void pointers to meaningful types
  char* s = (char*) s_v;
  ILIinputAux* aux_p = (ILIinputAux*) aux_vp;
  ILIinput* data_p = aux_p->data_p;

  if (aux_p->curr_row == 1) return;  // Ignore first row of the file
  if (aux_p->err_status) return;  // Do not parse if an error occurred before

  switch (aux_p->curr_col){
  case 0:
  case 1:
    // DO NOTHING. These columns are ignored.
    // OBS: case 0 should not be reached, columns are 1-based.
    break;

  case 2:
    data_p->year[data_p->size] = atoi(s);
    // TODO: error check (should be the same for other fields).
    break;

  case 3:
    data_p->week[data_p->size] = atoi(s);
    // TODO error check
    break;

  case 4:
    data_p->estInc[data_p->size] = atoi(s);
    // TODO error check
    break;
  
  default:
    // TODO: THROW ERROR: more columns than expected.
    break;
  }
  
  aux_p->curr_col++;
}


/* 
Callback function for each valid row that is read from file.
*/
void cb2 (int c, void *aux_vp) { 
  ILIinputAux* aux_p = (ILIinputAux*) aux_vp;
  ILIinput* data_p = aux_p->data_p;

  // Update cursors
  aux_p->curr_col = 1;
  if (aux_p->curr_row++ == 1) return;  // Update current row AND ignore if it's the first one.

  data_p->size++;  // If line was valid, increments the size of data containers.

  // Dynamical vector reallocation.
  if (data_p->size >= aux_p->capacity){
    // TODO: call reallocation using doubling criterion.
  }
}


// ------------------------------------------------------------------------------------------------
// HIGH-LEVEL INTERFACE FUNCTIONS
// ------------------------------------------------------------------------------------------------


/* 
Reads a csv file with ILI data. 

Assumes that the file has the following 4 columns:
"index", "year","week","est_Inc"

Where the first column ("index") is ignored. 
The first row of the file is assumed as index and is also ignored.

TODO: Function is unfinished. Error checking and reporting is not yet fully functional.


@param fname  Path for the csv file. Must be a null-terminated string.
@param data_p   Pointer to an ILIinput struct, to which the data is written.

@return An integer error code.

*/
int read_ili_csv(const char* fname, ILIinput* data_p){

    // Declarations
    // ------------
    FILE *fp;
    struct csv_parser p;
    char buf[FILE_BUF_SIZE];
    size_t bytes_read;
    unsigned char options = 0;
    const size_t reserve_size = 1024;  // Could be a parameter, but no optionals in C.

    // Initialization
    // --------------
  
    // Initialization of the auxiliary parser structure.
    ILIinputAux aux;
    aux.data_p = data_p;
    aux.capacity = reserve_size;
    aux.curr_row = aux.curr_col = 1;
    aux.err_status = 0;
    aux.err_row = aux.err_col = 0;

    // Initial allocation of the struct pointers
    data_p->year = malloc(reserve_size * sizeof(int));
    data_p->week = malloc(reserve_size * sizeof(int));
    data_p->estInc = malloc(reserve_size * sizeof(int));
    data_p->size = 0;

    // Initialiation of the parser
    if (csv_init(&p, options) != 0) {
        fprintf(stderr, "Failed to initialize csv parser\n");
        exit(EXIT_FAILURE);
    }

    // This can be ignored/commented, csvlib defines default functions for space/term inside the parser.
    csv_set_space_func(&p, is_space);
    csv_set_term_func(&p, is_term);

    // Set option to append null string terminator to each field
    options += CSV_APPEND_NULL;  
    csv_set_opts(&p, options); 


    // Execution
    // ---------

    // --- File opening
    fp = fopen(fname, "rb");
    if (!fp) {
        fprintf(stderr, "Failed to open %s: %s\n", fname, strerror(errno));
        exit(EXIT_FAILURE);
    }

    // --- Main loop for reading and parsing the file
    while ((bytes_read=fread(buf, 1, FILE_BUF_SIZE, fp)) > 0) {
        // if (csv_parse(&p, buf, bytes_read, cb1, cb2, &c) != bytes_read) {
        if (csv_parse(&p, buf, bytes_read, cb1, cb2, &aux) != bytes_read) {
        fprintf(stderr, "Error while parsing file: %s\n", csv_strerror(csv_error(&p)));
        }
    }


    // Final operations
    // ----------------
    csv_fini(&p, cb1, cb2, &aux);  // Closes the csv parser.

    if (ferror(fp)) {
        fprintf(stderr, "Error while reading file %s\n", fname);
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    fclose(fp);

    csv_free(&p);

    return EXIT_SUCCESS;
}

