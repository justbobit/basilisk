/**
#Double diffusion of a scalar with a moving embed boundary.

We simulate the diffusion of two tracers separated by an embedded boundary. The
interface moves using the Stefan relation :

$$
  \mathbf{v}_{pc} = \frac{1}{L_H}(\lambda_1 \nabla T_L - \lambda_1 \nabla T_S)
  $$
where $L_H$ is the latent heat of water, $T_L$ and $T_S$ are the temperature
fields on both sides of the interface.

The full algorithm is done on two iterations. It is the following :

  1. Calculate the gradients on both sides of the interface
  1. Move the interface
  1. Define the diffusion coefficient
  1. Diffuse one tracer
  1. Reinit the level_set function
  1. Exchange the inside and the outside of the embed boundary
  1. Define the diffusion coefficient
  1. Diffuse the other tracer
  1. Exchange the inside and the outside of the embed boundary

![Animation of cs*u.x + (1-cs)*u2.x.](mini_cell/visu.mp4)(loop)

![Animation of the velocity of the interface](mini_cell/v_pc.mp4)(loop)


~~~gnuplot Phase change velocity
f(x) = a + b*x
fit f(x) 'log' u (log($1)):(log($7)) via a,b
ftitle(a,b) = sprintf("%.3fexp^{{%4.2f}t}", a, b)
set xrange[5:120]
plot 'log' u 1:7 w l t 'Phase Change Velocity', exp(f(log(x))) t ftitle(a,b)
~~~

~~~gnuplot Temperature in the cell located at 
plot 'log' u 1:2 w l t 'Liquid Temperature',  'log' u 1:4 w l t 'Solid Temperature', \
    'log' u 1:5 w l t 'Approximated Temperature'
~~~

*/

#define DOUBLE_EMBED 1
#define LevelSet     1

#include "embed.h"
#include "advection.h"
#include "diffusion.h"
#include "../level_set.h"
#include "view.h"

#define MIN_LEVEL 6
#define MAXLEVEL 6
#define latent_heat 100.
#define T_eq         0.
#define TL_inf       1.
#define TS_inf       -1.

#define H0 -0.45*L0
#define DT_MAX  1.

#define T_eq         0.
#define plane(x, y, H) (y - H)

scalar TL[], TS[], dist[];
scalar * tracers = {TL, TS};
scalar * level_set = {dist};
face vector v_pc[];
face vector muv[];
mgstats mgT;
scalar grad1[], grad2[];

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
  
  N = 1 << MAXLEVEL;
  init_grid (N);
  run();
}

event init(t=0){


  DT = L0/(1 << MAXLEVEL);

  NB_width = nb_cell_NB * L0 / (1 << MAXLEVEL);

  lambda[0] = 1.;
  lambda[1] = 2.5;

  foreach_vertex(){
      dist[] = plane(x,y,H0);
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

  }
}

event properties(i++){
  

  foreach_face()
    muv.x[] = lambda[i%2]*fs.x[];
  boundary((scalar *) {muv});
}

event tracer_diffusion(i++){
  if(i%2==0){
    diffusion(TL, L0/(1 << MAXLEVEL), D = muv, theta = cs);
  }
  else{
    diffusion(TS, L0/(1 << MAXLEVEL), D = muv , theta = cs);
  }
}

event LS_reinitialization(i++,last){
  if(i>0 && i%2==0){
    LS_reinit2(dist,L0/(1 << MAXLEVEL), 1.4*NB_width,
      1.4*nb_cell_NB);
  }
}

#if DOUBLE_EMBED
event double_calculation(i++,last){
// modify the fs , cs, copy the outer fields in the partially covered cells

  foreach(){
    cs[]      = 1.-cs[];
  }
  foreach_face(){
    fs.x[]      = 1.-fs.x[];
  }

  boundary({cs,fs});
  restriction({cs,fs});
}
#endif

event movies ( i+=2,last;t<240.)
{
  if(t>0.1){
    if(i%20==0){
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
    squares("v_pc.y", min =0., max=0.06);
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
