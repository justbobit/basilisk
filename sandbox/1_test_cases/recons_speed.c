/**
# Reconstruction of a field off an interface


This is a test case for the reconstruction of a field only defined on an
interface off of it. It uses most of the functions related to the use of a
level set.
*/

#include "run.h"
#include "timestep.h"
#include "bcg.h"
#include "fractions.h"
#include "curvature.h"

#include "embed.h"
#include "view.h"
#include "../level_set.h"

#define MAXLEVEL 7


scalar dist[];
vector v_pc[];
face vector v_pc_f[];
scalar * level_set = {dist};

scalar * velocity = {v_pc.x,v_pc.y};

int nb_cell_NB;
double  NB_width ;              // length of the NB

/**
The initial interface defined by the level set function is a circle + a cosine
function.
*/

double geometry(double x, double y, double Radius) {

  coord center;
  center.x = 0.5;
  center.y = 0.5;

  double theta = atan2 (y-center.y, x-center.x);
  double R2  =  sq(x - center.x) + sq (y - center.y) ;
  double s = -( sqrt(R2*(1.+0.25*cos(8*theta))) - Radius);

  return s;
}


int main(){
	int N = 1 << MAXLEVEL;
	init_grid(N);
	nb_cell_NB = 1 << 3 ;
	NB_width = nb_cell_NB * L0 / (1 << MAXLEVEL);


  foreach_vertex() {
    dist[] = clamp(-geometry(x,y,L0/5.),-1.5*NB_width,1.5*NB_width);
  }
  LS_reinit2(dist,0.49*L0/(1 << MAXLEVEL), 
    1.2*NB_width,
    10);

  int iloop;

  for(iloop=1; iloop<=5000;iloop++){

  fprintf(stderr, "# %d\n", iloop);

  if(iloop%100 ==0){
	 	foreach_vertex() {
	    dist[] = clamp(dist[],-1.5*NB_width,1.5*NB_width);
	  }	
  }

  boundary ({dist});
  restriction({dist});
  fractions (dist, cs, fs);

  boundary({cs,fs});
  restriction({cs,fs});

  scalar curve[];
  curvature (cs,curve);
  boundary({curve});


  stats s = statsf(curve);
  fprintf(stderr, "%d %g %g\n", iloop, s.min, s.max);
  

  /**
  We set the velocity phase change field to be dependant on the curvature : 
  $$
  v_{pc} =  \{ -4.5 -\kappa * n_x, -4.5 -\kappa * n_y \}
  $$
  therefore the equilibrium will be reached once we reach a circle whose
  curvature is -4.5.
   */

  double maxk = max(fabs(s.min), fabs(s.max));
  foreach(){
    if(cs[] != 0. && cs[] != 1.){

      coord n = normal (point, cs);
      v_pc.x[] = (-4.5 -curve[])/maxk*n.x*L0/(1 << MAXLEVEL);
      v_pc.y[] = (-4.5 -curve[])/maxk*n.y*L0/(1 << MAXLEVEL);
     }
     else{
      v_pc.x[] = 0.;
      v_pc.y[] = 0.;
    }
  }

  
  
  boundary((scalar*){v_pc});

  double dt = L0/(1 << MAXLEVEL);

  recons_speed(dist, dt, nb_cell_NB, NB_width, velocity );


/**
Advection with the reconstructed speed.

We directly use $v_{pc}$ which is a cell-centered field, we could use a
face-centered field, but my first tries were unsuccessful.
*/   


  advection(level_set, v_pc, 0.8*L0 / (1 << MAXLEVEL));

  LS_reinit2(dist,L0/(1 << MAXLEVEL), 
  	1.2*NB_width,
    1);
  if( iloop%90==0){
  	view (fov = 20., quat = {0,0,0,1}, tx = -0.5, ty = -0.5, 
  		bg = {1,1,1}, width = 600, height = 600, samples = 1);
	  draw_vof("cs");
	  squares("curve", min =-30, max = 30);
	  save ("curve.mp4");

    draw_vof("cs");
    squares("dist", min =-NB_width, max = NB_width);
    save ("dist.mp4");
  clear();  

	  draw_vof("cs");
	  squares("v_pc.x", min =-0.2*L0/(1 << MAXLEVEL), 
	  	max = 0.2*L0/(1 << MAXLEVEL));
	  save ("vpcx.mp4");
  }
  if(iloop%500==0) output_facets (cs, stdout);	
	}
	
		dump();
}



/**
![Animation of the level set function](recons_speed/dist.mp4)(loop)

![Animation of the phase change velocity](recons_speed/vpcx.mp4)(loop)

~~~gnuplot Curvature evolution
set title 'Curvature' font 'Helvetica,20'
set key left
f(x) =  -4.5
plot 'log' u 1:2 w l t 'min',  \
     'log' u 1:3 w l  t 'max', \
     f(x) dt 2 t 'equilibrium curvature'

~~~

~~~gnuplot Evolution of the interface
set size ratio -1
plot 'out' w l t ''
~~~

## References

~~~bib

@article{peng_pde-based_1999,
  title = {A {PDE}-Based Fast Local Level Set Method},
  volume = {155},
  issn = {0021-9991},
  url = {http://www.sciencedirect.com/science/article/pii/S0021999199963453},
  doi = {10.1006/jcph.1999.6345},
  pages = {410--438},
  number = {2},
  journaltitle = {Journal of Computational Physics},
  author = {Peng, Danping and Merriman, Barry and Osher, Stanley and Zhao, Hongkai and Kang, Myungjoo},
  urldate = {2019-09-09},
  date = {1999-11}
}
~~~

*/
