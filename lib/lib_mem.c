#include <stdlib.h>
#include <stdio.h>
#include "lib_mem.h"

float *
mem_calloc_1d_float(size_t n, float v0, char *msg)
{
  if (n <= 0 ) {
    fprintf(stderr, "Error: size=%i is zero or negative (%s)!\n", n, msg);
        fflush(stderr);
    return NULL;
  }

  float *var = (float *) malloc( n * sizeof(float));
  if ( var == NULL ) {
    fprintf(stderr, "Error: can't malloc enough mem (%s)!\n", msg);
        fflush(stderr);
    exit(-1);
  }

  for (size_t i = 0; i < n; i++ ) var[i] = v0;

  return var;
}
