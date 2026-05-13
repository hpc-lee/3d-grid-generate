#ifndef CONSTANTS_H
#define CONSTANTS_H

// consts
#define CONST_MAX_STRLEN 1024
#define CONST_NDIM 3

#ifndef M_PI
#define PI 3.14159265358979323846264338327950288419716939937510
#else
#define PI M_PI
#endif

#define handle_nc_err(err)                       \
{                                                \
  if (err != NC_NOERR) {                         \
     fprintf(stderr,"nc error: %s\n", nc_strerror(err)); \
     exit(-1);                                   \
  }                                              \
}

#endif
