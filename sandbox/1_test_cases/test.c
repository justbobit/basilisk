#include "embed.h"
#include "navier-stokes/centered.h"
// #include "navier-stokes/perfs.h"
#include "basic_geom.h"
#include "view.h"


face vector muv[];


double geometry (double x, double y){
	
  coord center_rectangle, size_rectangle;
  center_rectangle.x = L0/2.-0.3;
  center_rectangle.y = L0/2.;

  size_rectangle.x = 0.015;
  size_rectangle.y = 0.25; 

  double r[10];


  double s = 1.;
  int i;
  for(i = 0; i<=7; i++){
  	r[i] = rectangle(x, y, center_rectangle, size_rectangle);
  	center_rectangle.x += 0.081;
  	s = min (s,r[i]);
  }
  return s;
}


#define MAXLEVEL 10

int main() {
  int N = 1 << MAXLEVEL;
	init_grid (N);
	run();
}



event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*0.125/900.;
}

/**
The fluid is injected on the left boundary with a unit velocity. The
tracer is injected in the lower-half of the left boundary. An outflow
condition is used on the right boundary. */

u.n[left]   = dirichlet(1.);
u.t[left]   = dirichlet(1.);
u.n[bottom] = dirichlet(1.);
u.t[bottom] = dirichlet(1.);

p[left]     = neumann(0.);
pf[left]    = neumann(0.);

u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
u.n[top]   = neumann(0.);
u.t[top]   = neumann(0.);

p[right]  = dirichlet(0.);
pf[right] = dirichlet(0.);
p[top]    = dirichlet(0.);
pf[top]   = dirichlet(0.);



/**
The top and bottom walls are free-slip and the cylinder is no-slip. */

u.n[embed] = fabs(y) > 1. ? neumann(0.) : dirichlet(0.);
u.t[embed] = fabs(y) > 1. ? neumann(0.) : dirichlet(0.);

event init (t = 0)
{

  /**
  The domain is the intersection of a channel of width a fan of 7 rectangles */
  L0 = 4.;
  vertex scalar phi[];
  foreach_vertex() {
    phi[] = geometry(x,y);
  }
  boundary ({phi});
  fractions (phi, cs, fs);

  /**
  We set the initial velocity field. */
  
	view (fov = 20., quat = {0,0,0,1}, tx = -0.5, 
		ty = -0.5, bg = {1,1,1}, width = 600, height = 600, samples = 1);
	draw_vof("cs");	
  squares("phi",min = -0.001, max= 0.001);
  save ("phi.png");

  foreach_face(){
  	if(cs[]>0.)u.x[] = cs[] ? 1. : 0.;
  }
}


/**
We check the number of iterations of the Poisson and viscous
problems. */

event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

/**
We produce animations of the vorticity and tracer fields... */

event movies (t+=0.0035; t <= 6.)
{
  scalar omega[], m[];
  vorticity (u, omega);
  foreach()
    m[] = cs[] - 0.5;
  boundary ({m});
 	view (fov = 10., quat = {0,0,0,1}, tx = -0.578524, ty = -0.604008, 
 		bg = {1,1,1}, width = 600, height = 600, samples = 1);
  draw_vof("cs");	
  squares("omega",min = -200, max = 200, linear = true);
  save ("vort.mp4");  
}

/**
We adapt according to the error on the embedded geometry, velocity and
tracer fields. */

event adapt (i++) {
  adapt_wavelet ({cs,u}, (double[]){1e-2,3e-2,3e-2}, 10, 4);
}