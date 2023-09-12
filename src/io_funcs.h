#ifndef IO_FUNCS_H
#define IO_FUNCS_H

#include "gd_t.h"

/*************************************************
 * structure
 *************************************************/
typedef struct
{
  int nx;
  int nz;
  float *var; // pointer to var
} io_quality_t;


/*************************************************
 * function prototype
 *************************************************/

int
init_io_quality(io_quality_t *io_quality, gd_t *gdcurv);

int
gd_curv_coord_export(gd_t *gdcurv, char *output_dir);

int
quality_export(io_quality_t *io_quality, char *output_dir, char *var_name);

int
io_get_nextline(FILE *fp, char *str, int length);

#endif
