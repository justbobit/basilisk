#ifndef BASILISK_HEADER_19
#define BASILISK_HEADER_19
#line 1 "./../alex_functions.h"
/**
## My set of functions 


A function to rescale normals so that they are unit vectors w.r.t. the
2-norm (by default, the 1-norm is adopted for efficiency purposes). */

coord normal (Point point, scalar c) {
  coord n = mycs (point, c);
  double nn = 0.;
  foreach_dimension()
    nn += sq(n.x);
  nn = sqrt(nn);
  foreach_dimension()
    n.x /= nn;
  return n;
}

/**
A function to compute 2-norm normal in every cell. */

void compute_normal (scalar f, vector normal_vector) {
  foreach() {
    coord n = normal (point, f);
    foreach_dimension() 
      normal_vector.x[] = n.x;
  }
  boundary((scalar*){normal_vector});
}

/** 
#sign function
*/
double sign2 (double x)
{
  return(x > 0. ? 1. : x<0 ? -1. : 0.);
}

#endif
