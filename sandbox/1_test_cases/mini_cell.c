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
f(x)  = a  - x/b
f2(x) = a2 - x/b2  
f3(x) = a3 - x/b3  
fit f(x)  'log1' u ($1):(log($7)) via a ,b
fit f2(x) 'log2' u ($1):(log($7)) via a2,b2
fit f3(x) 'log3' u ($1):(log($7)) via a3,b3
ftitle(a,b) = sprintf("%.3fexp^{t/{%4.4f}}",a,b)
latent(a) = sprintf("L_H  %4.4f",a)
set grid
set xlabel "time"
set ylabel "log(v_{pc})"
#set logscale y
aa = 300.
set term pngcairo size 1024,800
plot 'log1' u 1:(log($7)) w p pointtype 2  t latent(aa/1),\
         f(x) w l lw 3 t ftitle(a,b), \
    'log2' u 1:(log($7)) w p pointtype 2  t latent(aa/2), \
        f2(x) t ftitle(a2,b2), \
    'log3'  u 1:(log($7)) w p pointtype 2 lc rgb "blue" t latent(aa/3), \
        f3(x) t ftitle(a3,b3)
~~~

~~~gnuplot Temperature in the cell located at 
plot 'log1' u 1:2 w l lw 3t 'Liquid Temperature',  'log1' u 1:4 w l lw 3 t \
'Solid Temperature',  'log' u 1:5 w l lw 3 t 'Approximated Temperature'
~~~

*/

#define DOUBLE_EMBED 1
#define LevelSet     1

#include "embed.h"
#include "advection.h"
#include "diffusion.h"
#include "../level_set.h"
#include "view.h"

#define T_eq         0.
#define TL_inf       1.
#define TS_inf       -1.

int MIN_LEVEL, MAXLEVEL; 
double H0;
double latent_heat;
char filename [100];
FILE * fp1;

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
  
mgstats mg1,mg2;

TL[embed]  = dirichlet(T_eq);
TL[top]    = dirichlet(TL_inf); 

TS[embed]  = dirichlet(T_eq);
TS[bottom] = dirichlet(TS_inf); 
int j;
/**
The domain is 4 units long, centered vertically. */

int main() {
  L0 = 4.;
  CFL = 0.5;
  origin (-L0/2., -L0/2.);
  
  j = 1;
  for (j=1;j<=3;j++){
    latent_heat  = 300./j;
    MAXLEVEL = MIN_LEVEL = 5;
    H0 = -0.3*L0;
    N = 1 << MAXLEVEL;
    snprintf(filename, 100,  "log%d", j);
    fp1 = fopen (filename,"w");

    init_grid (N);
    run();
    // fclose(fp1); 
  }
}

event init(t=0){


  DT = L0/(1 << MAXLEVEL);

  NB_width = nb_cell_NB * L0 / (1 << MAXLEVEL);

  lambda[0] = 1.;
  lambda[1] = 1.;

  foreach(){
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

  foreach_face(){
    v_pc.x[] = 0.;
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

/** 
An error is made on the estimation of the volume of the cell, only for certain
cells. Here we try a simple correction propertional to that error which is :
$$
\epsilon_V = v_{pc}*dt*\Delta - (cs^{n+1}-cs^n)*\Delta^2
$$

The best way to correct that error would be using an Arbitrary Lagrangian
Eulerian formulation that integrates volume fluxes. TBD.
*/
    // foreach(){

      // if(cs0[] ==1. && cs[] !=1.){
        // TL[] +=  9.*(v_pc.y[]*dt/Delta-(cs[]-cs0[])) ;
        // fprintf(stderr, "# 1 %g %g\n", y, v_pc.y[]*dt/Delta);
      // }
      // if((cs0[]!=0. && cs[] ==0.)){
        // TS[] +=  (v_pc.y[]*dt/Delta-cs0[]) ;
        // fprintf(stderr, "# 2 %g %g %g\n", y, 
          // cs0[], v_pc.y[]*dt/Delta) ;
      // }
    // } 
  }
}


event properties(i++){
  

  foreach_face()
    muv.x[] = lambda[i%2]*fs.x[];
  boundary((scalar *) {muv});
}

event tracer_diffusion(i++){
  if(i%2==0){
    mg1 = diffusion(TL, L0/(1 << MAXLEVEL), D = muv);
  }
  else{
    mg2 = diffusion(TS, L0/(1 << MAXLEVEL), D = muv );
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

event movies ( i+=2,last;t<130.)
{
  if(t>20.){
    if(i%20==0 && j ==3){
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
    squares("v_pc.y", min =0.01, max=0.015);
    save ("v_pc.mp4");
    }

    stats s2    = statsf(v_pc.y);
    Point p     = locate(-0.5*L0/(1<<MIN_LEVEL),-1.51*L0/(1<<MIN_LEVEL));

    double cap  = capteur(p, TL);
    double cap3 = capteur(p, TS);
    double cap2 = capteur(p, cs);

    double T_point = (1.-cap2)*cap + (cap2)*cap3;
    fprintf (fp1, "%.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",
      t, cap, cap2, cap3, T_point, s2.min, s2.max);
    fprintf(fp1, "## %g %g\n", mg1.resa, mg2.resa);
  }
}

// event dumps(t=25.){
//   dump();
// }
