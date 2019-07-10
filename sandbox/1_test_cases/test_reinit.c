#include "grid/cartesian.h"
#include "advection.h"
#include "../level_set.h"
#include "../basic_geom.h"


double NB_width;
int    nb_cell_NB =  1 << 2 ;  // number of cells for the NB

scalar dist[];
scalar * tracers = {dist};
scalar * level_set = NULL;
int MAXLEVEL = 7;


int main() {
  origin (-L0/2., -L0/2.);
  N = 1 << MAXLEVEL;
  run(); 
}

event init (i = 0)
{
	DT=1.;
	NB_width = L0*nb_cell_NB / (1<<MAXLEVEL);
	coord center_circle;
	double Radius = 0.125/L0;
  center_circle.x = 0.;
  center_circle.y = 0.;


  vertex scalar dist[];
  foreach_vertex()
    dist[] = circle(x,y,center_circle,Radius);
  boundary ({dist});

  scalar l[],l2[];
	foreach(){
		double l1 = 0;
		foreach_dimension()
		  l1 += pow((dist[]    - dist[-1,0])/Delta,2.);
		l[] = sqrt(l1);
		l2[] = dist[];
		// fprintf (stderr, "%f \n", l[]);
	}
  output_ppm (l, file = "grad_dist_init.png", n = 512,min = 0.95, max = 1.05);
  output_ppm (l2, file = "dist_init.png", n = 512,min = -0.1, max = 0.1);
		fprintf (stderr, "TOP \n");
  LS_reinit2(dist,L0/(1 << MAXLEVEL), 0.8*NB_width,
    1.4*nb_cell_NB);
	foreach(){
		double l1 = 0;
		foreach_dimension()
		  l1 += pow((dist[]    - dist[-1,0])/Delta,2.);
		l[] = sqrt(l1);
		l2[] = dist[];
		// fprintf (stderr, "%f \n", l[]);
	}
	output_ppm (l,  file = "grad_dist_reinit.png", n = 512,min = 0.95, 
		max = 1.05);
	output_ppm (l2, file = "dist_reinit.png", n = 512,min = -0.1, max =
		0.1);
}