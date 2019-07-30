/**
#Issue with gradient calculation and embedded boundaries

Example displaying an issue for embedded boundaries close to a mesh face when
the embedded is parallele to the interface.

*/

#include "embed.h"
#include "diffusion.h"
#include "view.h"

#define MAXLEVEL 5
#define T_eq     0.
#define TL_inf  -1.
#define T_eq     0.

#define plane(x, y, H) (y - H)

scalar TL[] ;
TL[embed]  = dirichlet(T_eq);
TL[bottom] = dirichlet(TL_inf); 

/**
The domain is 4 units long, centered vertically. */

int main() {
  L0 = 4.;
  CFL = 0.5;
  origin (-L0/2., -L0/2.);


  double T[2];
  T[0] = -(5.-0.00011)*L0/(1 << MAXLEVEL);
  T[1] = -(5.-0.5)*L0/(1 << MAXLEVEL);
  
  N = 1 << MAXLEVEL;
  init_grid (N);

  vertex scalar dist[];
  int ii, iteration = 100;
  int ij;
  double T1[iteration][2], T2[iteration][2];
  for (ii = 1; ii <= 2; ii++){
    char name[80];
    sprintf (name, "Temperature%d.mp4", ii);
  
    foreach_vertex(){
      dist[] = -plane(x,y,T[ii-1]);
    }
    boundary ({dist});
    restriction({dist});
    fractions (dist, cs, fs);
    boundary({cs,fs});
    restriction({cs,fs});

    double dt = 1., t = 0;
    foreach() {
      TL[] = T_eq;
    }

    boundary({TL});
    restriction({TL});
    for (ij=1;ij<iteration;ij++){

      diffusion(TL, dt, fs);
      t ++;

      boundary({TL});
      restriction({TL});

      draw_vof("cs");
      cells();
      squares("TL", min =-1., max=0.);
      save (name);
      stats s2 = statsf(TL);
      if(ii==1){
        T1[ij-1][0] = s2.min;
        T1[ij-1][1] = s2.max;
      }
      else{
        T2[ij-1][0] = s2.min;
        T2[ij-1][1] = s2.max;
      }
    }
  }
  for (ij=1;ij<iteration;ij++){
    fprintf(stderr, "%d %g %g %g %g\n", ij, 
      T1[ij-1][0], T1[ij-1][1], T2[ij-1][0], T2[ij-1][1]);
  }
}
/**
## Results

![Temperature field, embedded boundary close to a mesh face]
(mini_cell/Temperature1.mp4)

![Temperature field, embedded boundary at the middle of a cell]
(mini_cell/Temperature2.mp4)

~~~gnuplot Temperature profiles
set term svg
set xlabel 'Iteration count'
set ylabel 'Temperature'
set xrange [*:*]
set yrange [-1.5:4]
plot 'log' u 1:2 w l t 'T1_{min}', \
     'log' u 1:3 w l t 'T1_{max}', \
     'log' u 1:4 w l t 'T2_{min}', \
     'log' u 1:5 w l t 'T2_{max}' 
~~~
*/