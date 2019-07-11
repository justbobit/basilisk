/**
# An advection solver

We wish to solve the advection equations
$$
\partial_tf_i+\mathbf{u}\cdot\nabla f_i=0
$$
where $\mathbf{u}$ is the velocity field and $f_i$ are a list of
passive tracers.  This can be done with a flux-based advection scheme
such as the 2nd-order, unsplit, upwind scheme of [Bell-Collela-Glaz,
1989](references.bib#bell89).

The main time loop is defined in [run.h](). A stable timestep needs to
respect the [CFL
condition](http://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition). */

#include "run.h"
#include "timestep.h"

/**
We allocate the (face) velocity field. For compatibility with the
other solvers, we allocate it as `uf` and define an alias. The
`gradient` function is used to set the type of slope-limiting
required. The default is to not use any limiting (i.e. a purely
centered slope estimation). */

face vector uf[];
#define u uf

double limiter_general (double s0, double s1, double s2)
{
  double the; 
  #if MINMOD
     the = 1.;
  #elif SUPERBEE
     the = 2.;
  #else
     the = 0;
  #endif 
  double d1,d2,d3;
  d1 = the*(s1-s0);
  d2 = (s2-s0)/2.;
  d3 = the*(s2-s1);
  if(s0 < s1 && s1 < s2) {
    if(d2 < d1) d1 = d2; 
    return min(d1,d3);
  }
  if (s0 > s1 && s1 > s2) {
    if(d2 > d1) d1 = d2;
      return max(d1,d3);
  }
  return(0);  
}

#if MINMOD
  double (* gradient) (double, double, double) = limiter_general;
#elif SUPERBEE
  double (* gradient) (double, double, double) = limiter_general;
#else
  double (* gradient) (double, double, double) = NULL;
#endif

/**
Here we set the gradient functions for each tracer (as defined in the
user-provided `tracers` list). */

extern scalar * tracers;

event defaults (i = 0) {
  for (scalar f in tracers)
    f.gradient = gradient;
}

/**
We apply boundary conditions after user initialisation. */

event init (i = 0) {
  boundary ((scalar *){u});
  boundary (tracers);
}

/**
The timestep is set using the velocity field and the CFL
criterion. The integration itself is performed in the events of
[tracer.h](). */

event velocity (i++,last) {
  dt = dtnext (timestep (u, DT));
}

#include "tracer.h"
