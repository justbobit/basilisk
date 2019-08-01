/**
#Double embed boundary calculation

First test case with a double calculation inside and outside of an embedded boundary.
This is based on the example of the Bénard–von Kármán Vortex Street for flow 
around a cylinder at Re=160, I modified the geometry for fun and use a
Zalesak's notched disk. Both fluids have the same properties. We use cs for
the cell fraction in one phase and 1-cs in the other.

At the end of each iteration of the solver we swap the variables, change the fs,
the cs and do a solver iteration on the new set of variables.
We use the centered Navier-Stokes solver, with embedded boundaries and advect
the passive tracer *f*. 

I still have issues regarding the use of the same timestep for both phases. They
hare truly independant. Once this is done I will work on making both phases
fields dependant of one another and then make the boundary movement dependant on
both fields.


![Animation of cs*u.x + (1-cs)*u2.x.](moving_embed/visu.mp4)(loop)

![Animation of the velocity of the interface](moving_embed/v_pc.mp4)(loop)

*/

#define DOUBLE_EMBED  1
#define LevelSet      1

#include "embed.h"
#include "../centered_alex.h"
#include "../level_set.h"
#include "diffusion.h"
#include "tracer.h"
#include "view.h"
#include "alex_functions.h"

#define MIN_LEVEL 6
#define MAXLEVEL 6
#define latent_heat 10.
#define T_eq         0.
#define TL_inf       1.
#define TS_inf       -1.



#define H0 -L0/3.-1.001*L0/(1 << MAXLEVEL)
#define dH_refine (2.*L0/(1 << 5))
#define DT_MAX  1.

#define T_eq         0.
#define plane(x, y, H) (y - H)



scalar TL[], TS[], dist[];
scalar * tracers = {TL, TS};
scalar * level_set = {dist};
face vector v_pc[];
face vector muv[];
mgstats mgT;


double lambda[2];


int     nb_cell_NB =  1 << 3 ;  // number of cells for the NB
double  NB_width ;              // length of the NB
  

TL[embed]  = dirichlet(T_eq);
TL[top]    = dirichlet(TL_inf); 

TS[embed]  = dirichlet(T_eq);
TS[bottom] = dirichlet(TS_inf); 


/**
The domain is 4 units long, centered vertically. */

int main() {
  L0 = 4.;
  CFL = 0.5;
  origin (-L0/2., -L0/2.);
  
  N = 1 << MIN_LEVEL;
  init_grid (N);

  NB_width = L0*nb_cell_NB / (1<< MAXLEVEL);
  mu = muv;
  run();
}

/**
We set a constant viscosity corresponding to a Reynolds number of 160,
based on the cylinder diameter (0.125) and the inflow velocity (1). */

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*0.125/160.;
  boundary((scalar *) {muv});
}

event init (t = 0)
{

  DT = CFL*L0 / (1 << MAXLEVEL);


  lambda[0] = 1.;
  lambda[1] = 1.;

  /**
  The domain is the intersection of a channel of width unity and a
  circle of diameter 0.125. */
  foreach_vertex(){
    dist[] = clamp(plane(x,y,H0),-NB_width,NB_width);
  }
  boundary ({dist});
  restriction({dist});
  fractions (dist, cs, fs);
  boundary({cs,fs});
  restriction({cs,fs});
  
  foreach() {
    TL[] = T_eq;
    TS[] = T_eq;
  }

  boundary({TL,TS});
  restriction({TL,TS});

}

event stability(i++){
    double dtmax2 = DT_MAX;
    timestep (uf, dtmax2);
}


event tracer_advection(i++,last){
  if(i%2 == 0){
    double L_H       = latent_heat;  
    phase_change_velocity_LS_embed (cs, fs ,TL, TS, v_pc, dist, L_H, 1.05*NB_width,
      nb_cell_NB,lambda);

    scalar cs0[];
    foreach()
      cs0[] = cs[];
    boundary({cs0});
    restriction({cs0});

    advection_LS (level_set, v_pc, dt);
    boundary ({dist});
    restriction({dist});

    fractions (dist, cs, fs);
    boundary({cs,fs});
    restriction({cs,fs});

    event ("properties");
  }
}

event tracer_diffusion(i++){

  foreach_face()
    muv.x[] = lambda[i%2]*fs.x[];
  boundary((scalar *) {muv});

  if(i%2==0){
    mgT = diffusion(TL, dt, fs);
  }
  else{
    mgT = diffusion(TS, dt, muv);
  }
}


event LS_reinitialization(i++,last){
  if(i>0 && i%2==0){
    LS_reinit2(dist,L0/(1 << MAXLEVEL), 1.4*NB_width,
      1.4*nb_cell_NB);
  }
}


/**
We produce an animation of the tracer field. */

event movies ( i+=2,last;t<500.)
{
  if(t>0.1){
    if(i%40==0){
    boundary({TL,TS});
    scalar visu[];
    foreach(){
      visu[] = (1.-cs[])*TL[]+cs[]*TS[] ;
    }
    boundary({visu});
    
    draw_vof("cs");
    squares("visu", min =-1., max = 1.);
    save ("visu.mp4");
    
    boundary((scalar *){v_pc});
    clear();
    draw_vof("cs");
    cells();
    squares("v_pc.y", min =-3.e-3, max=3.e-3);
    save ("v_pc.mp4");
    }
    stats s2    = statsf(v_pc.y);
    Point p     = locate(-0.5*L0/(1<<MIN_LEVEL),-1.51*L0/(1<<MIN_LEVEL));

    double cap  = capteur(p, TL);
    double cap3  = capteur(p, TS);
    double cap2 = capteur(p, cs);
    double T_point = (1.-cap2)*cap + (cap2)*cap3;
    fprintf (stderr, "%.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",
      t, cap, cap2, cap3, T_point, s2.min, s2.max);
  }  
}

/**
~~~gnuplot Temperature in one cell near the interface
plot 'log' u 1:2 w l t 'Liquid Temperature',  'log' u 1:4 w l t 'Solid Temperature', \
    'log' u 1:5 w l t 'Approximated Temperature'
~~~

~~~gnuplot Phase change velocity
f(x) = a + b*x
fit f(x) 'log' u (log($1)):(log($7)) via a,b
ftitle(a,b) = sprintf("%.3f/x^{%4.2f}", exp(a), -b)
set xrange[2:500]
set yrange[0:9.e-3]
plot 'log' u 1:7 w l t 'Phase Change Velocity', exp(f(log(x))) t ftitle(a,b)
~~~

*/
 