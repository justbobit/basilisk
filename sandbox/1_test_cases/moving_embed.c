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
scalar grad1[], grad2[];


int     nb_cell_NB =  1 << 3 ;  // number of cells for the NB
double  NB_width ;              // length of the NB
  

TL[embed]  = dirichlet(T_eq);
TL[top]    = dirichlet(TL_inf); 

TS[embed]  = dirichlet(T_eq);
TS[bottom] = dirichlet(TS_inf); 

void phase_change_velocity_LS_embed (scalar cs, face vector fs, scalar tr,
 scalar tr2, face vector v_pc, scalar dist, double L_H, 
 double NB_width) {
  
 /**
  The phase change velocity $\mathbf{v}_{pc}$ is

  $$
  \mathbf{v}_{pc} = \mathrm{Pe}\, D\, \nabla tr
  $$
  
  here we use the embed_gradient_face_x defined in embed that gives a proper
  definition of the gradients with embedded boundaries. */
  vector gtr[], gtr2[];
  foreach(){
    foreach_dimension(){
      gtr.x[] = 0.;
    }
    if(fabs(dist[])<NB_width){
      if ( ( (fs.x[] !=0. && fs.x[] != 1.) || (fs.y[] !=0. && fs.y[] != 1.)   ) 
              && fabs(dist[])<=0.9*NB_width){
        coord n       = facet_normal( point, cs ,fs) , p;
        normalize(&n);
        double alpha  = plane_alpha (cs[], n);
        line_length_center (n, alpha, &p);
        double c=0.;
        double temp = 0.;
        // foreach_dimension(){
        //   temp += v_pc.x[]*v_pc.x[];
        // }
        // temp = -0.2*sqrt(temp);
        double grad = dirichlet_gradient(point, tr, cs , n, p, temp, &c);
        foreach_dimension(){
          gtr.x[] += grad*n.x;
        }
      }
    }
    grad2[] = gtr.y[];
  }

  boundary((scalar*){gtr});
  foreach(){
      cs[]      = 1.-cs[];
  }
  foreach_face()
    fs.x[]      = 1.-fs.x[];

  boundary({cs,fs});
  restriction({cs,fs});
  
  foreach(){
    foreach_dimension(){
      gtr2.x[] = 0.;
    }
    if(fabs(dist[])< 0.9*NB_width){
      if(((fs.x[] !=0. && fs.x[] != 1.) || (fs.y[] !=0. && fs.y[] != 1.))){
        coord n       = facet_normal( point, cs ,fs) , p;
        normalize(&n);
        double alpha  = plane_alpha (cs[], n);
        line_length_center (n, alpha, &p);
        double c=0.;
        double temp = 0.;
        // foreach_dimension(){
        //   temp += v_pc.x[]*v_pc.x[];
        // }
        // temp = -0.2*sqrt(temp);
        double grad = dirichlet_gradient(point, tr2, cs , n, p, temp, &c);
        foreach_dimension(){
          gtr2.x[] += grad*n.x;
        }
      }
    }
    grad1[] = gtr2.y[];
  }
  boundary((scalar*){gtr2});
  
  boundary({grad1,grad2});

  foreach(){
      cs[]      = 1.-cs[];
  }
  foreach_face()
    fs.x[]      = 1.-fs.x[];
  boundary({cs,fs});
  restriction({cs,fs});

  /**
  With the the normal vector and the gradients of the tracers we can now compute
  the phase change velocity $\mathbf{v}_{pc}$, following the lines drawn in
  [meanflow.c](/sandbox/popinet/meanflow.c). We define it as the product between
  the density ratio and the diffusive flow. Note that $\mathbf{v}_{pc}$ is 
  weighted by the face metric. */
  
  face vector v_pc2[];
  foreach_face() {

    v_pc.x[]  = 0.;
    v_pc2.x[] = 0.;
  }


  foreach_face(){
    if ( (  fs.y[] !=0. && fs.y[] != 1.   ) 
              && fabs(dist[])<=0.9*NB_width){
      v_pc.x[]     =  (1.4*gtr2.x[] - gtr.x[])*0.02;
      v_pc2.x[]    = v_pc.x[];
    }
  }

  boundary((scalar *){v_pc});
  boundary((scalar *){v_pc2});

  int ii;
  for (ii=1; ii<=nb_cell_NB; ii++){
    foreach(){
      if(fabs(dist[])<0.9*NB_width){

        coord grad_dist;
        if(dist[]>0.){
        foreach_dimension(){
          double a = max(0.,dist[]    - dist[-1,0]);
          double b = min(0.,dist[1,0] - dist[]    );
          grad_dist.x = max(a,b);
          }
        }
        else{
          foreach_dimension(){
            double a = min(0.,dist[]    - dist[-1,0]);
            double b = max(0.,dist[1,0] - dist[]    );
           grad_dist.x = max(a,b);
          }
        }
        int k1, k2;
        int sig[2] = {1 , 1 };
        if(dist[]>0.){
          if(grad_dist.x>0) sig[0] = -1;
          if(grad_dist.y>0) sig[1] = -1;
        }
        else{
          if(grad_dist.x<0) sig[0] = -1;
          if(grad_dist.y<0) sig[1] = -1;
        }
        
        double ratio = fabs(grad_dist.x)/(max(SEPS,fabs(grad_dist.y))) ;
        if(ratio > sqrt(3.)){
          k1 = 1;
          k2 = 0;
        }
        else{
          if(ratio > 1./sqrt(3.)){
            k1 = 1;
            k2 = 1;   
          }
        else{
          k1 = 0;
          k2 = 1; 
          }
        }
        v_pc2.x[] = (v_pc.x[sig[0]*k1,sig[1]*k2])*
                    exp( -((fabs(dist[]-grad_dist.x*Delta/2.))/(4*Delta)));
        v_pc2.y[] = v_pc.y[sig[0]*k1,sig[1]*k2]*
                    exp( -((fabs(dist[]-grad_dist.y*Delta/2.))/(4*Delta)));
      }
    }
    boundary((scalar *){v_pc2});
    restriction((scalar *){v_pc2});    


    foreach_face(){
      if(v_pc.x[]==0. && v_pc2.x[] !=0.) v_pc.x[] = v_pc2.x[];
    }

    boundary((scalar *){v_pc});
    restriction((scalar *){v_pc});
  }
}


void tracer_fluxes_LS (scalar f,
        face vector uf,
        face vector flux,
        double dt,
        (const) scalar src)
{

  /**
  We first compute the cell-centered gradient of *f* in a locally-allocated
  vector field. */
  
  vector g[];
  gradients ({f}, {g});

  /**
  For each face, the flux is composed of two parts... */

  foreach_face() {

    /**
    A normal component... (Note that we cheat a bit here, `un` should
    strictly be `dt*(uf.x[i] + uf.x[i+1])/((fm.x[] +
    fm.x[i+1])*Delta)` but this causes trouble with boundary
    conditions (when using narrow '1 ghost cell' stencils)). */

    double un = dt*uf.x[]/(Delta + SEPS), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = f[i] + (src[] + src[-1])*dt/4. + s*(1. - s*un)*g.x[i]*Delta/2.;

    /**
    and tangential components... */

  
    double vn = (uf.y[i] + uf.y[i,1])/2.;
    double fyy = vn < 0. ? f[i,1] - f[i] : f[i] - f[i,-1];
    f2 -= dt*vn*fyy/(2.*Delta);
   

    flux.x[] = f2*uf.x[];
  }

  /**
  Boundary conditions ensure the consistency of fluxes across
  variable-resolution boundaries (on adaptive meshes). */

  boundary_flux ({flux});
}


void advection_LS (struct Advection p) //never takes into account the solid
// boundary, the velocity field is defined on both phases.
{

  /**
  If *src* is not provided we set all the source terms to zero. */
  
  scalar * lsrc = p.src;
  if (!lsrc) {
    const scalar zero[] = 0.;
    for (scalar s in p.tracers)
      lsrc = list_append (lsrc, zero);
  }

  assert (list_len(p.tracers) == list_len(lsrc));
  scalar f, src;
  for (f,src in p.tracers,lsrc) {
    face vector flux[];
    tracer_fluxes_LS (f, p.u, flux, p.dt, src);
    foreach()
      foreach_dimension()
        f[] += p.dt*(flux.x[] - flux.x[1])/(Delta); // careful we have removed
        // cm[]
  }
  boundary (p.tracers);

  if (!p.src)
    free (lsrc);
}



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
  // #if TREE
  //   refine (level < MAXLEVEL && plane(x,  y, (H0 - dH_refine)) > 0.
  //           && plane(x, y, (H0 + dH_refine)) < 0.);
  // #endif


  DT = CFL*L0 / (1 << MAXLEVEL);
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
    phase_change_velocity_LS_embed (cs, fs ,TL, TS, v_pc, dist, L_H, 1.05*NB_width);

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

  if(i%2==0){

  foreach_face()
    muv.x[] = fs.x[];
  boundary((scalar *) {muv});

    mgT = diffusion(TL, dt, fs);
  }
  else{

  foreach_face()
    muv.x[] = fs.x[];
  boundary((scalar *) {muv});

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

event movies ( i+=40,last;t<500.)
{
  if(t>2.){
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
set xrange[2:1200]
plot 'log' u 1:7 w l t 'Phase Change Velocity', exp(f(log(x))) t ftitle(a,b)
~~~

*/
 