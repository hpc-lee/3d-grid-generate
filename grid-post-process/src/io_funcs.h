#ifndef IO_FUNCS_H
#define IO_FUNCS_H

#include "gd_t.h"
#include "par_t.h"

/*************************************************
 * function prototype
 *************************************************/

int
gd_curv_coord_export(gd_t *gdcurv, par_t *par);

int 
read_import_coord(gd_t *gdcurv, par_t *par);

int
io_get_nextline(FILE *fp, char *str, int length);

#endif
