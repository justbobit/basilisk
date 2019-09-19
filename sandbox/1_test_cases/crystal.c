
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
int 	MAXLEVEL = 6; 
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
                                    

double 	latent_heat = 10000.;
double 	lambda[2]; 		// thermal capacity of each material
double 	epsK = 0.0005, epsV = 0.0005;
int 		nb_cell_NB;
double  NB_width ;    // length of the NB

scalar curve[];

TL[embed]  = dirichlet(T_eq -epsK*curve[]-epsV*sqrt(v_pc.x[]*v_pc.x[]+v_pc.y[]*v_pc.y[]));
TL[top]    = dirichlet(TL_inf); 
TL[bottom] = dirichlet(TL_inf); 
TL[left]   = dirichlet(TL_inf); 
TL[right]  = dirichlet(TL_inf); 

TS[embed]  = dirichlet(T_eq -epsK*curve[]-epsV*sqrt(v_pc.x[]*v_pc.x[]+v_pc.y[]*v_pc.y[]));
TS[top]    = dirichlet(TS_inf); 
TS[bottom] = dirichlet(TS_inf); 
TS[left]   = dirichlet(TS_inf); 
TS[right]  = dirichlet(TS_inf); 

/**
Initial geometry definition
*/
double geometry(double x, double y, double Radius) {

  coord center;
  center.x = 0.5;
  center.y = 0.5;

  double theta = atan2 (y-center.y, x-center.x);
  double R2  =  sq(x - center.x) + sq (y - center.y) ;
  double s = -( sqrt(R2*(1.+0.25*cos(8*theta))) - Radius);
  // double s = -( sqrt(R2) - Radius);

  return s;
}

/**
mystery variable !
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
  int N = 1 << MAXLEVEL;
	init_grid (N);
	run();
}


event init(t=0){
	DT         = L0 / (1 << MAXLEVEL); 	// Delta
	nb_cell_NB = 1 << 3 ; 							// number of cell in the 
																			// narrow band 
	NB_width   = nb_cell_NB * L0 / (1 << MAXLEVEL);

  lambda[0] = 1.;
  lambda[1] = 1.;

  foreach_vertex() {
    dist[] = clamp(-geometry(x,y,L0/5.),-1.5*NB_width,1.5*NB_width);
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
    TL[] = -1;
    TS[] = -1.;
  }
  boundary({TL,TS});
  restriction({TL,TS});

  scalar visu[];
  foreach(){
    visu[] = (1.-cs[])*TL[]+(cs[])*TS[] ;
  }
  boundary({visu});

  
}

event properties(i++){
  foreach_face()
    muv.x[] = lambda[i%2]*fs.x[];
  boundary((scalar *) {muv});
}

event tracer_diffusion(i++){
  int kk;
  for (kk=1;kk<=(40*k_loop+1);kk++){
    if(i%2==0){
      mg1 = diffusion(TL, 0.5*L0/(1 << MAXLEVEL), D = muv);
    }
    else{
      mg2 = diffusion(TS, 0.5*L0/(1 << MAXLEVEL), D = muv );
    }
  }
}

event LS_advection(i++,last){
  if(i%30 == 1){
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
      nb_cell_NB,lambda);

	  recons_speed(dist, 0.5*DT, nb_cell_NB, NB_width, LS_speed);

    dt = timestep_LS (v_pc, DT);

    stats s = statsf (v_pc.y);
		fprintf(stderr, "# VELOCITY %g %g %g %g\n", dt, s.min, 
			s.max, DT);	  
		view (fov = 20., quat = {0,0,0,1}, tx = -0.5, ty = -0.5, 
			bg = {1,1,1}, width = 600, height = 600, samples = 1);
	  draw_vof("cs");
	  squares("v_pc.y", min = s.min, max = s.max);
	  save ("vpcy.png");
		
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


event movies ( i++,last;i<1000)
{
  if(i%10 == 1) {
    // if(i%20==1 && j ==1){
    boundary({TL,TS});
    scalar visu[];
    foreach(){
      visu[] = (cs[])*TL[]+(1.-cs[])*TS[] ;
    }
    boundary({visu});
    
    draw_vof("cs");
    squares("visu", min =-0.2, max = 0.2);
    save ("visu.mp4");
    
    boundary((scalar *){v_pc});

    scalar TT[];

    foreach(){
    	TT[] = 0.;
    }
    boundary({TT});

#if Gibbs_Thomson
  double epsK = 0.005;
  double epsV = 0.005;
	  scalar curve[];
	  curvature (cs,curve);
	  boundary({curve});
	  foreach(){
	    if(curve[] != nodata) TT[] = -epsK*curve[]-
	    	-epsV*sqrt(v_pc.x[]*v_pc.x[]+v_pc.y[]*v_pc.y[]);
	  }
	  boundary({TT});
	  stats s = statsf(TT);
	  fprintf(stderr, "EQUILIBRIUM TEMP %g %g\n", s.min, s.max);
#endif
	 	draw_vof("cs");
    squares("v_pc.x", min =-NB_width, max = NB_width);
    save ("v_pcx.mp4");

  }
}
