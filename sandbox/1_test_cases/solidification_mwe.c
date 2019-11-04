/**
#Solidification minimum working example

We simulate the diffusion of two tracers separated by an embedded boundary. The
interface moves using the Stefan relation :

$$
  \mathbf{v}_{pc} = \frac{1}{L_H}(\lambda_1 \nabla T_L - \lambda_2 \nabla T_S)
  $$
where $L_H$ is the latent heat of water, $T_L$ and $T_S$ are the temperature
fields on both sides of the interface.

The full algorithm is done on two iterations can be found on the mini_cell test
case.

![Animation of cs*u.x + (1-cs)*u2.x.](solidification_mwe/visu.mp4)(loop)

We output only the interface at different times during the simulation.


![](solidification_mwe/v_pc_1D.png)

~~~gnuplot Phase change velocity
set term pngcairo size 1200,1200 font "Helvetica,20"
set output 'v_pc_1D.png'
set title 'v_{pc} for fixed mesh, grid is 2^6x2^6, 1 Poisson ite, tolerance varies'
f(x)  = a/sqrt(x+c) - b
f2(x)  = a2/sqrt(x+c2) - b2
f3(x)  = a3/sqrt(x+c3) - b3 
fit f(x)  'log1' u ($1):($7) via a ,b,c
fit f2(x) 'log2' u ($1):($7) via a2,b2,c2
fit f3(x) 'log3' u ($1):($7) via a3,b3,c3
ftitle(a,c,b) = sprintf('%.4f/(t+%.1f)^{0.5}-%.5f',a,c,b)
latent(a) = sprintf("tolerance = %4.4e",a)
set grid
set xlabel "time"
set ylabel "v_{pc}"
set logscale y
aa = 1e-9
plot 'log1' u 1:7 w l lw 1 t latent(aa/8),\
         f(x) w l dt 2 lw 1 t ftitle(a,c,b), \
    'log2' u 1:7 w l lw 2 t latent(aa/8/8), \
        f2(x) dt 2 lw 1 t ftitle(a2,c2,b2 ), \
    'log3'  u 1:7 w l  lw 1 lc rgb "blue" t latent(aa/8/8/8), \
        f3(x) dt 2 lw 1 t ftitle(a3,c3,b3)
set output 'erreur.png'
unset xrange
unset logscale y
set yrange [-0.0002:0.005]
plot 'log1' u 1:($7-f($1)) w l lw 3 t 'err1', \
     'log2' u 1:($7-f2($1)) w l lw 3 t 'err2',\
     'log3' u 1:($7-f3($1)) w l lw 3 t 'err3'
~~~

~~~gnuplot Evolution of the interface (zoom)
set output 'interface_1D.png'
unset yrange
plot 'out' w l t ''
~~~

*/

#define DOUBLE_EMBED 1
#define LevelSet     1
#define Gibbs_Thomson 0

#include "embed.h"
#include "advection.h"
#include "diffusion.h"

#include "fractions.h"
#include "curvature.h"

#include "view.h"
#include "../level_set.h"

#define T_eq          0.
#define TL_inf        1.
#define TS_inf       -1.

#define tstart 5.

int MIN_LEVEL, MAXLEVEL; 
double H0;
double latent_heat;
char filename [100];
FILE * fp1;

#define DT_MAX  1.

#define T_eq         0.


#define plane(x, y, H) (y  - H)
// #define plane(x, y, H) (y - H)


scalar TL[], TS[], dist[];
scalar * tracers = {TL, TS};
scalar * level_set = {dist};

vector v_pc[];
scalar * LS_speed   = {v_pc.x,v_pc.y};
face vector muv[];
mgstats mgT;
scalar grad1[], grad2[];

#if Gibbs_Thomson
double  epsK = 0.005, epsV = 0.005;
#else
double  epsK = 0.000, epsV = 0.000;
#endif
scalar curve[];



double lambda[2];

int     nb_cell_NB =  1 << 3 ;  // number of cells for the NB
double  NB_width ;              // length of the NB


  
mgstats mg1,mg2;

TL[embed] = dirichlet(T_eq + epsK*curve[]-epsV*sqrt(v_pc.x[]*v_pc.x[]+v_pc.y[]*v_pc.y[]));
TL[top]   = dirichlet(TL_inf); 

TS[embed]  = dirichlet(T_eq + epsK*curve[]-epsV*sqrt(v_pc.x[]*v_pc.x[]+v_pc.y[]*v_pc.y[]));
TS[bottom] = dirichlet(TS_inf); 

int j;
int k_loop = 0;
/**
The domain is 4 units long, centered vertically. */

double timestep_LS (const face vector u, double dtmax)
{
  static double previous = 0.;
  dtmax /= CFL;
  foreach_face(reduction(min:dtmax))
    if (u.x[] != 0.) {
      double dt = Delta/fabs(u.x[]);
      if (dt < dtmax) dtmax = dt;
    }
  dtmax *= CFL;
  if (dtmax > previous)
    dtmax = (previous + 0.1*dtmax)/1.1;
  previous = dtmax;
  return dtmax;
}


int main() {
  L0 = 4.;
  CFL = 0.5;
  origin (-L0/2., -L0/2.);
  
  j = 1;
  for (j=1;j<=1;j++){

/**
Here we set up the parameters of our simulation. The latent heat $L_H$, the
initial position of the interface $h_0$ and the resolution of the grid.
*/
    latent_heat  = 50.;
    MAXLEVEL = MIN_LEVEL = 6 ;

    H0 = -0.2*L0; 
    N = 1 << MAXLEVEL;
    snprintf(filename, 100,  "log%d", j);
    fp1 = fopen (filename,"w");

    init_grid (N);
    run();
    fclose(fp1); 
  }
}

event init(t=0){

  TOLERANCE = 1.e-6;
  DT = latent_heat/50.*L0/(1 << MAXLEVEL);

  NB_width = nb_cell_NB * L0 / (1 << MAXLEVEL);

  lambda[0] = 1.;
  lambda[1] = 1.4;

  foreach(){
      dist[] = clamp(plane(x,y,H0),-1.1*NB_width, 1.1*NB_width);
  }
  boundary ({dist});
  restriction({dist});

  fractions (dist, cs, fs);
  boundary({cs,fs});
  restriction({cs,fs});

  curvature(cs,curve);
  boundary({curve});

  double y_max = -L0;

  vector h[];
  foreach(){
    if(interfacial(point, cs)){
      heights (cs, h);
    }
  }
  boundary({h});
  double hh0 = H0;
  foreach(reduction(max:y_max)){
    if(interfacial(point, cs)){
      double yy = y+Delta*height(h.y[]);
      fprintf(stderr, "%g %g %g %g\n", yy,Delta*h.y[],L0/(1 << MAXLEVEL),hh0);
      // y_max = max(y_max,y+Delta*height(h.y[]));
    }     
  }
  exit(1);
  fprintf(stderr, "%g %g\n",t, y_max);

  // stats s2    = statsf(curve);
  // view (fov = 10.4411, quat = {0,0,0,1}, tx = -0.223746, 
  //     ty = -0.00297483, bg = {1,1,1}, width = 600, 
  //     height = 600, samples = 1);
  //   draw_vof("cs");
  //   squares("curve", min =s2.min, max = s2.max);
    // save ("curve_init.png");
    // fprintf(stderr, "%g %g\n", s2.min, s2.max);

  foreach() {
    TL[] = 1.;
    TS[] = -1.;
  }

  foreach_face(){
    v_pc.x[] = 0.;
  }

  boundary({TL,TS});
  restriction({TL,TS});
}


event tracer_diffusion(i++){
  int kk;
  foreach_face()
    muv.x[] = lambda[i%2]*fs.x[];
  boundary((scalar *) {muv});
  for (kk=1;kk<=1;kk++){
    if(i%2==0){
      mg1 = diffusion(TL, L0/(1 << MAXLEVEL), D = muv);
    }
    else{
      mg2 = diffusion(TS, L0/(1 << MAXLEVEL), D = muv );
    }
  }
}


event LS_advection(i++,last){
  if(i%2 == 1 ){
    double L_H       = latent_heat;  

    scalar cs0[];


    foreach(){
      cs0[]   = cs[];
      cs[]    = 1.-cs[];
    }
    foreach_face(){
      fs.x[]  = 1.-fs.x[];
    }

    boundary({cs,fs,cs0});
    restriction({cs,fs});

    phase_change_velocity_LS_embed (cs, fs ,TL, TS, v_pc, dist, L_H, 1.05*NB_width,
      nb_cell_NB,lambda,epsK, epsV);

    recons_speed(dist, 0.5*DT, nb_cell_NB, NB_width, LS_speed);

    stats s2 = statsf (v_pc.y);

    dt = timestep_LS (v_pc, DT);

    
    boundary((scalar *){v_pc});
    if(t>tstart){

      advection_LS (level_set, v_pc, dt);
      boundary ({dist});
      restriction({dist});

      fractions (dist, cs, fs);
      boundary({cs,fs});
      restriction({cs,fs});

      curvature(cs,curve);
      boundary({curve});

      foreach(){
        cs[]      = 1.-cs[];
      }
      foreach_face(){
        fs.x[]      = 1.-fs.x[];
      }

      boundary({cs,fs});
      restriction({cs,fs});

      k_loop = 0;
      foreach(){
        if( (cs0[] != 1. && cs[] ==1.) || (cs0[] == 0. && cs[] !=0.))k_loop = 1;
      }
    }
  }
}


event LS_reinitialization(i++,last){
  if(i>0 && i%2==1){
    LS_reinit2(dist,L0/(1 << MAXLEVEL), 1.4*NB_width,
      1);
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

event movies ( i++,last;t<120.)
{
  if(i%(40*j) == 1 && t > tstart) {
    // stats s2    = statsf(curve);

    boundary({TL,TS});
    scalar visu[];
    foreach(){
      visu[] = (cs[])*TL[]+(1.-cs[])*TS[] ;
    }
    boundary({visu});
    
    // view (fov = 20.5678, quat = {0,0,0,1}, tx = 0.0308985, 
    //   ty = 0.0325629, bg = {1,1,1}, width = 800, height = 800, samples = 1);    

    // draw_vof("cs");
    // squares("visu", min =-1., max = 1.);
    // save ("visu.mp4");

    if(i%400==1 && j == 3) {
      output_facets (cs, stdout);
    }
  }
  if(i%2 == 1 && t>tstart){
    stats s3    = statsf(v_pc.y);
    Point p     = locate(-1.51*L0/(1<<MIN_LEVEL),-3.51*L0/(1<<MIN_LEVEL));

    double cap  = capteur(p, TL);
    double cap3 = capteur(p, TS);
    double cap2 = capteur(p, cs);
    double T_point = cap2*cap + (1.-cap2)*cap3;
    fprintf (fp1, "%.9f %.9f %.9f %.9f %.9f %.9f %.9f\n",
      t, cap, cap2, cap3, T_point, s3.min, s3.max);

    fprintf(stderr, "## %g %g %d\n", mg1.resa, mg2.resa, k_loop);
  }
}