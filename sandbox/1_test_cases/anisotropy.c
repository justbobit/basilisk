/**
#Melt of a solid particle

We simulate the diffusion of two tracers separated by an embedded boundary. The
interface moves using the Stefan relation :

$$
  \mathbf{v}_{pc} = \frac{1}{L_H}(\lambda_1 \nabla T_L - \lambda_1 \nabla T_S)
  $$
where $L_H$ is the latent heat of water, $T_L$ and $T_S$ are the temperature
fields on both sides of the interface.

The boundaries of the domain are set at $T_L = 1$. The ice particle is initially
at $T_S = -1$.

*/

#define DOUBLE_EMBED 	1
#define LevelSet     	1
#define Gibbs_Thomson 1

#include "embed.h"
#include "advection.h"
#include "diffusion.h"

#include "fractions.h"
#include "curvature.h"

#include "../level_set.h"
#include "view.h"

#define T_eq         0.
#define TL_inf       -1.
#define TS_inf       -1.

/**
Setup of the numerical parameters
*/
int 	MAXLEVEL = 7; 
double 	H0;

/**
Setup of the physical parameters + level_set variables
*/
scalar TL[], TS[], dist[];
vector v_pc[];

scalar * tracers 		= {TL, TS};
scalar * level_set 	= {dist};
scalar * LS_speed 	= {v_pc.x,v_pc.y};
face vector muv[];
mgstats mgT;
scalar grad1[], grad2[];
scalar curve[];
                                    

double 	latent_heat = 1000.;
double 	lambda[2]; 		// thermal capacity of each material
#if Gibbs_Thomson
double  epsK = 0.02, epsV = 0.;

TL[embed] = dirichlet(T_eq + (epsK*curve[]-epsV*sqrt(v_pc.x[]*v_pc.x[]+v_pc.y[]*v_pc.y[]))*fac1(x,y));
TS[embed] = dirichlet(T_eq + (epsK*curve[]-epsV*sqrt(v_pc.x[]*v_pc.x[]+v_pc.y[]*v_pc.y[]))*fac1(x,y));


#else
double  epsK = 0.000, epsV = 0.000;
TL[embed]  = dirichlet(T_eq);
TS[embed]  = dirichlet(T_eq);


#endif

int 		nb_cell_NB;
double  NB_width ;    // length of the NB


TL[top]    = dirichlet(TL_inf); 
TL[bottom] = dirichlet(TL_inf); 
TL[left]   = dirichlet(TL_inf); 
TL[right]  = dirichlet(TL_inf); 

TS[top]    = dirichlet(TS_inf); 
TS[bottom] = dirichlet(TS_inf); 
TS[left]   = dirichlet(TS_inf); 
TS[right]  = dirichlet(TS_inf); 

/**
Initial geometry definition. Here the interface equation is :

$$
r\left(1+ 0.15 *cos(6\theta) \right) - R
$$
where $r = \sqrt{x^2 + y^2}$, $\theta = \arctan(x,y)$ and $R = \frac{L_0}{5}$

Notice that the initial dist[] field is not really a distance, it is modified
after a few iterations of the LS_reinit() function.
*/
double geometry(double x, double y, double Radius) {

  coord center;
  center.x = 0.5;
  center.y = 0.5;

  double theta = atan2 (y-center.y, x-center.x);
  double R2  =  sq(x - center.x) + sq (y - center.y) ;
  // double s = -( sqrt(R2)*(1.+0.15*cos(6*theta)) - Radius);
  double s = -( sqrt(R2) - Radius);

  return s;
}

/**
$k_{loop}$ is a variable used when the interface arrives in a new cell. Because
we are using embedded boundaries, the solver needs to converge a bit more when
the interface is moved to a new cell. Therefore, the Poisson solver is used
$40*k_{loop} +1$ times.
*/
int k_loop = 0;


/**
Output variables
*/
mgstats mg1,mg2;


/**
Special routines
*/
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

  TOLERANCE = 5.e-8;
  int N = 1 << MAXLEVEL;
	init_grid (N);
	run();
}


event init(t=0){
	DT         = 0.3*L0 / (1 << MAXLEVEL); 	// Delta
	nb_cell_NB = 1 << 3 ; 							// number of cell in the 
																			// narrow band 
	NB_width   = nb_cell_NB * L0 / (1 << MAXLEVEL);

  lambda[0] = 1.;
  lambda[1] = 1.;

  foreach_vertex() {
    dist[] = clamp(-geometry(x,y,L0/8.),-1.5*NB_width,1.5*NB_width);
  }

  boundary ({dist});
  restriction({dist});

  fractions (dist, cs, fs);
  boundary({cs,fs});
  restriction({cs,fs});


  foreach_face(){
    v_pc.x[] = 0.;
  }
  boundary((scalar *){v_pc});

  curvature(cs,curve);
  boundary({curve});

  foreach() {
    TL[] = 0.;
    TS[] = 0.;
    foreach_dimension()
      v_pc.x[] = 0.;
  }
  boundary({TL,TS});
  restriction({TL,TS});
  
}

event properties(i++){
  foreach_face()
    muv.x[] = lambda[i%2]*fs.x[];
  boundary((scalar *) {muv});
}

event tracer_diffusion(i++){
  int kk;
  for (kk=1;kk<=(10*k_loop+1);kk++){
    if(i%2==0){
    	boundary({TL});
      mg1 = diffusion(TL, L0/(1<< MAXLEVEL), D = muv);
    }
    else{
    	boundary({TS});
      mg2 = diffusion(TS, L0/(1<< MAXLEVEL), D = muv );
    }
  }
}

event LS_advection(i++,last){
  if(i%2 ==1 && i > 200){
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

	  recons_speed(dist, 0.9*L0/(1 << MAXLEVEL), nb_cell_NB, NB_width, LS_speed);

    dt = timestep_LS (v_pc, DT);

    stats s = statsf (v_pc.y);
		fprintf(stderr, "%g %g %g \n", t, s.min, s.max);	  

    advection_LS (level_set, v_pc, dt);
    
   	boundary ({dist});
    restriction({dist});

    fractions (dist, cs, fs);
    boundary({cs,fs});
    restriction({cs,fs});

	  curvature(cs,curve);
	  boundary({curve});

    stats s2    = statsf(curve);
    fprintf(stderr, "%g %g\n", s2.min, s2.max);

    foreach(){
      cs[]      = 1.-cs[];
    }
    foreach_face(){
      fs.x[]      = 1.-fs.x[];
    }

    boundary({cs,fs});
    restriction({cs,fs});
    // event ("properties");

    k_loop = 0;
    foreach(){
      if(cs0[] != 1. && cs[] ==1.)k_loop = 1;
    }
  }

}

event LS_reinitialization(i++,last){
  if(i>0 && i%2==1){
    LS_reinit2(dist,0.5*L0/(1 << MAXLEVEL), 
  	1.2*NB_width,
    4);
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



event movies ( i++,last;t<30.)
{
  if(i%40 == 1) {
    // if(i%20==1 && j ==1){
    boundary({TL,TS});
    scalar visu[];
    foreach(){
      visu[] = (cs[])*TL[]+(1.-cs[])*TS[] ;
    }
    boundary({visu});
    // view (fov = 4.9, quat = {0,0,0,1}, tx = -0.64, ty = -0.62, 
    // 	bg = {1,1,1}, width = 800, height = 800, samples = 1);
	  
		view (fov = 20., quat = {0,0,0,1}, tx = -0.5, ty = -0.5, 
    	bg = {1,1,1}, width = 800, height = 800, samples = 1);
	  

    draw_vof("cs");
    squares("visu", min =-0.1, max = 0.1);
    save ("visu.mp4");
    stats s = statsf(v_pc.x);
    stats s2 = statsf(v_pc.y);

    draw_vof("cs");
    squares("curve", min =-5., max = 5.);
    save ("curve.mp4");
    
    boundary((scalar *){v_pc});

    scalar TT[];

    foreach(){
    	TT[] = 0.;
    }
    boundary({TT});

    fprintf(stderr, "# %g %g %g %g %g\n", t, s.min, s.max, s2.min, s2.max);
  }
}

/**

![Animation of the approximated temperature field](anisotropy/visu.mp4)(loop)

![Animation of the curvature](anisotropy/curve.mp4)(loop)

~~~gnuplot Temperature on the boundary
set key left
set yrange [-0.003:0.003]
plot 'log' u 1:2 w l t 'min',  \
     'log' u 1:3 w l  t 'max'
~~~
*/