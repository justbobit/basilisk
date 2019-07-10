/**
# Zalesak's notched disk

This classical test advects a notched disk from [Zalesak (1979)](https://www.sciencedirect.com/science/article/pii/0021999179900512). The interface should come back to its original position.
The difference between the initial and final shapes is 
a measure of the errors accumulated during advection.

We will need the advection solver combined with the VOF advection
scheme. */

#include "utils.h"
#include "distance.h"
#include "tag.h"
#include "advection.h"
#include "vof.h"
#include "../level_set.h"
#define Pi 3.141592653589793
#define Nb_tour 2
#define T 2*Nb_tour*3.14
#define N_display 1 << 8
// #include "alimare/elementary_body.h"

double NB_width;
int     nb_cell_NB =  1 << 2 ;  // number of cells for the NB

/**
The volume fraction is stored in scalar field `f` which is listed as
an *interface* for the VOF solver. We do not advect any tracer with
the default (diffusive) advection scheme of the advection solver. */

scalar f[], dist[];
scalar * interfaces = {f}, * tracers = {dist};
scalar * level_set = NULL;

coord * ZalesakCurve(void) {
/**
#Definition of the analytic curve of the 0-level set function

We trace a polyline with 2 points on the circular part and then add the last two
points. Of course we start and end with the same point. Thus, we have properly
defined a closed contour which can be used for the distance() function.

*/
  Array * h = array_new();
  
  double theta_i = acos(sqrt(1-sq(0.025/0.15)));
  double theta_f = 2.*Pi-theta_i;
  double delta_theta = 0.5*2.*Pi/360.;
  double theta;

  coord A = {0.525, 0.62};
  coord B = {0.475, 0.85};
  coord C = {0.525, 0.85};
  coord D = A;
  {
    coord toto = {A.x,A.y} ;
    array_append(h, &toto,sizeof(coord));
  }
  for(theta=theta_i; theta<theta_f; theta+=delta_theta){
    double xx = 0.5  + 0.15*sin(theta);
    double yy = 0.75 - 0.15*cos(theta);
    coord toto = {xx,yy};
    array_append(h, &toto,sizeof(coord));
    array_append(h, &toto,sizeof(coord));
  }
  {  
    coord toto = {B.x, B.y};
    array_append(h, &toto, sizeof(coord));
    array_append(h, &toto, sizeof(coord));
  }
  {
    coord toto = {C.x, C.y};
    array_append(h, &toto, sizeof(coord));
    array_append(h, &toto, sizeof(coord));
  }
  {
    coord toto = {D.x,D.y} ;
    array_append(h, &toto,sizeof(coord));
  }

  coord p = {nodata};
  array_append(h, &p, sizeof(coord));

  return (coord *) array_shrink(h);
}

int MAXLEVEL;

/**
We center the unit box on the origin and set a maximum timestep of 0.1 */

int main() {
  
  /**
  We then run the simulation for different levels of refinement. */

  for (MAXLEVEL = 7; MAXLEVEL <= 9; MAXLEVEL++) {
    init_grid (1 << MAXLEVEL);
    run();
  }
}

/**
We define the levelset function $\phi$ on each vertex of the grid and
compute the corresponding volume fraction field. */

event init (i = 0) {
  DT = T/(1280*Nb_tour);
  NB_width = L0*nb_cell_NB / (1<<MAXLEVEL);
  coord * p = ZalesakCurve();
  scalar d[];
  distance (d, p);
  boundary({d});

  // while (adapt_wavelet ({d}, (double[]){1e-3}, MAXLEVEL+2).nf);


  foreach_vertex()
    dist[] = (d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
  boundary ({dist});
  face vector s[];
  fractions (dist, f, s);

  foreach_vertex(){
    dist[] = -clamp(dist[], -1.02*NB_width, 1.02*NB_width);
  }
  boundary({dist});
  fprintf (stderr,"N %d %d \n", N, N_display);
  if (1) {
    scalar l[];
    foreach()
      l[] = dist[];
    output_ppm (l, file = "dist_init.png", n = 400, min = -0.1*NB_width,
     max=0.1*NB_width);
  }


}

/**
This event defines the velocity field. The disk rotates around the center of
the grid (0.5,0.5).

On trees we first adapt the grid so that the estimated error on
the volume fraction is smaller than $5\times 10^{-3}$. We limit the
resolution at `MAXLEVEL` and we only refine the volume fraction field
`f`.
*/

event velocity (i++) {

#if TREE
  adapt_wavelet ({f,dist}, (double[]){1.e-3, 1.e-3*NB_width}, MAXLEVEL, 
    list= {f,dist});
#endif

  foreach_face(x)
    u.x[] = Pi/3.14*(0.5-y);

  foreach_face(y)
    u.y[] = Pi/3.14*(x-0.5);
  boundary ((scalar *){u});
}

/**
At the start and end of the simulation we check the sum, min and max
values of the volume fraction field. The sum must be constant to
within machine precision and the volume fraction should be bounded by
zero and one. */

event logfile (t = {0,T}) {
  stats s = statsf (f);
  fprintf (stderr, "#REZ %f %.12f %.9f %g \n", t, s.sum, s.min, s.max);
  stats s2 = statsf (dist);
  fprintf (stderr, "#REZ %f %.12f %.9f %g\n", t, s2.sum, s2.min, s2.max);
}


/**
To compute the error, we reinitialise field `e` at the end of the
simulation with the initial shape and compute the difference with the
final shape. We output the norms as functions of the maximum
resolution `N`. */

event field (t = T) {
  scalar e[], e2[];
  coord * p = ZalesakCurve();
  scalar d[];
  distance (d, p);

  foreach_vertex()
    e2[] = (d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
  boundary ({e2});
  face vector s[];
  fractions (e2, e, s);

  boundary({dist});

  foreach(){
    e[]  -= f[];
    e2[] -= fabs(dist[])< NB_width/2. ? dist[] : e2[];
  }
  norm n  = normf (e);
  norm n2 = normf (e2);
  fprintf (stderr, "%d %g %g %g %g %g %g\n", N, n.avg, n.rms, n.max, 
            n2.avg, n2.rms, n2.max);
}

#if TREE
event levels (t = {T}) {
  if (N == N_display) {
    scalar l[];
    foreach()
      l[] = dist[];
    output_ppm (l, file = "dist.png", n = 400, min = -0.1*NB_width,
     max=0.1*NB_width);

  }
}
#endif

#if 0 //movie if needed
event levels2 (t += T/(100*Nb_tour)){
  if (N == N_display) {
    scalar l[];
    // foreach()
    //   l[] = level;
    // output_ppm (l, file = "levels2.gif", n = 400, min = 0, max = 7);

    foreach()
      l[] = f[];
    output_ppm (l, file = "f_reversed2.gif", n = 800, 
      opt = "--delay 1",min = 0, max = 1);

    foreach()
      l[] = fabs(dist[]) < NB_width ? 1. : 0.;
    output_ppm (l, file = "NB.gif", n = 400, 
      opt = "--delay 1",min = 0, max = 1);

    foreach()
      l[] = dist[];
    output_ppm (l, file = "dist2.gif", n = 400, 
      opt = "--delay 1",min = -0.1*NB_width, max=0.1*NB_width);
  }
}
#endif 

event LS_reinitialization(i+=4,last){
  if(i>15){
    LS_reinit2(dist,L0/(1 << MAXLEVEL), 0.9*NB_width,
      1.4*(nb_cell_NB << 1));
  }
}

/**
We output the interfaces of the vof variable.
*/
event shape (t = {0,T}) {
  if (N == N_display) output_facets (f);
}

/**
## Results

To be completed


~~~gnuplot Convergence rates for adaptive grids

ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)

f(x)=a+b*x
fit f(x) 'log' u (log($1)):(log($4)) via a,b
f2(x)=a2+b2*x
fit f2(x) 'log' u (log($1)):(log($2)) via a2,b2

f3(x)=a3+b3*x
fit f3(x) 'log' u (log($1)):(log($7)) via a3,b3
f4(x)=a4+b4*x
fit f4(x) 'log' u (log($1)):(log($5)) via a4,b4

set xlabel 'Maximum resolution'
set ylabel 'Maximum error'
set key bottom left
set logscale
set xrange [64:1024]
set xtics 64,2,1024
set grid ytics
set cbrange [1:1]

plot 'log' u 1:4 t 'VOF max (cart)'  , exp(f(log(x)))       t ftitle(a ,b ), \
     'log' u 1:2 t 'VOF norm1 (cart)', exp(f2(log(x))) lw 3 t ftitle(a2,b2), \
     'log' u 1:7 t 'LS max (cart)'   , exp(f3(log(x)))      t ftitle(a3,b3), \
     'log' u 1:5 t 'LS norm1 (cart)' , exp(f4(log(x))) lw 3 t ftitle(a4,b4)

~~~

The order of convergence is 1.5 for the VOF variable and 2 for the level set
function. WIP. Due to the initial shape of the disk ? Both field are deformed 
to the left.


~~~gnuplot Shapes of the interface - VOF variable
reset
set size ratio -1
plot [0.25:0.75][0.5:1]'out' w l t "adaptive"
~~~

Shape

Initial Level set field    |  Final Level Set field
:-------------------------:|:-------------------------:
![Initial](zalesak/dist_init.png)(width="400" height="300") |  ![Final](zalesak/dist.png)(width="400" height="300")  


 */
