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


![Animation of cs*u.x + (1-cs)*u2.x.](update/movie.mp4)(loop)
*/

#define DOUBLE_EMBED  1
#define LevelSet      1

#include "embed.h"
#include "../centered_alex.h"
#include "../level_set.h"
#include "diffusion.h"
#include "tracer.h"
#include "view.h"

#define MIN_LEVEL 7
#define MAXLEVEL 7
#define latent_heat 10.
#define T_eq         0.
#define TL_inf       1.
#define TS_inf       -1.


#define H0 0.001*L0
#define DT_MAX  1.

#define T_eq         0.
#define plane(x, y, H) (y - H)



scalar TL[], TS[], dist[];
scalar * tracers = {TL, TS};
scalar * level_set = {dist};
face vector v_pc[];
face vector muv[];
mgstats mgT;
scalar grad1[], grad2[] ;

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
  
  face vector gtr[], gtr2[];
  foreach_face(){
    if(fabs(dist[])<NB_width)
      gtr.x[] = face_gradient_x(tr,0);
    else
      gtr.x[] = 0.;
  }


  boundary((scalar*){gtr});
  foreach(){
      cs[]      = 1.-cs[];
  }
  foreach_face()
    fs.x[]      = 1.-fs.x[];

  boundary({cs,fs});
  restriction({cs,fs});
  
  foreach_face(){
    if(fabs(dist[])<NB_width)
      gtr2.x[] = face_gradient_x(tr2,0);
    else
      gtr2.x[] = 0.;
  }
  boundary((scalar*){gtr2});
  
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

  foreach_face(y){
    grad1[] = gtr.y[];
    grad2[] = gtr2.y[];
  }

  coord n = {0, 1.}; // here the normal must be properly defined.
  foreach_face(){
    if ( ( (fs.x[] !=0. && fs.x[] != 1.) || (fs.y[] !=0. && fs.y[] != 1.)   ) 
              && fabs(dist[])<=0.9*NB_width){

      v_pc.x[] = (4.*gtr.x[]*n.x
                  - gtr2.x[]*n.x)*0.01;
      v_pc2.x[] = v_pc.x[];
    }
  }
  boundary((scalar *){v_pc});

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
      v_pc2.x[] = v_pc.x[sig[0]*k1,sig[1]*k2];
      v_pc2.y[] = v_pc.y[sig[0]*k1,sig[1]*k2];
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
    // tracer_fluxes (f, p.u, flux, p.dt, src);
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
  
  N = 1 << MAXLEVEL;
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
  DT = L0 / (1 << MAXLEVEL);
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
    if(t>0.1){  
  phase_change_velocity_LS_embed (cs, fs ,TL, TS, v_pc, dist, L_H, 1.05*NB_width);
    }
    else{
      foreach_face()
        v_pc.x[] = 0.;

      foreach(){
        grad1[] = 0.;
        grad2[] = 0.;
      }
    }

    advection_LS (level_set, v_pc, dt);
    boundary ({dist});
    restriction({dist});

    scalar cs0[];
    foreach()
      cs0[] = cs[];
    boundary({cs0});
    restriction({cs0});


    fractions (dist, cs, fs);
    boundary({cs,fs});
    restriction({cs,fs});

    event ("properties");
  }
}

event tracer_diffusion(i++){

  foreach_face()
    muv.x[] = fm.x[]*0.125/160.;
  boundary((scalar *) {muv});

  if(i%2==0){
    mgT = diffusion(TL, dt, D = 1.);
  }
  else{
    mgT = diffusion(TS, dt, D = 1.);
  }
}


event LS_reinitialization(i+=2,last){
  if(i>0){
    LS_reinit2(dist,L0/(1 << MAXLEVEL), 1.4*NB_width,
      3.*nb_cell_NB);
  }
}


/**
We produce an animation of the tracer field. */

event movies ( t +=0.2,last)
{
  boundary({TL,TS});
  scalar visu[];
  foreach(){
    visu[] = cs[]*TL[]+(1.-cs[])*TS[] ;
    // visu[] = TL[] ;
  }
  boundary({visu});
  view (fov = 16.642, quat = {0,0,0,1}, tx = -0.0665815, 
    ty = -0.00665815, bg = {1,1,1}, width = 600, height = 600, samples = 1);
  draw_vof("cs");
  squares("visu", min =-1, max = 1);
  save ("visu.mp4");
  view (fov = 5., quat = {0,0,0,1}, tx = 0, ty = 0, 
    bg = {1,1,1}, width = 600, height = 600, samples = 1);
  draw_vof("cs");
  squares("grad1", min =-1, max = 1);
  save ("grad1.mp4");
  view (fov = 5., quat = {0,0,0,1}, tx = 0, ty = 0, 
    bg = {1,1,1}, width = 600, height = 600, samples = 1);
  draw_vof("cs");
  squares("grad2", min =-1, max = 1);
  save ("grad2.mp4");

  view (fov = 5., quat = {0,0,0,1}, tx = 0, ty = 0, 
    bg = {1,1,1}, width = 600, height = 600, samples = 1);
  boundary({dist});
  draw_vof("cs");
  squares("dist", min =-NB_width, max = NB_width);
  save ("dist.mp4");
  

  clear();

  draw_vof("cs");
  cells();
  squares("v_pc.y", min =-0.02, max=0.02);
  save ("v_pc.mp4");
}

/**
We check the number of iterations of the Poisson and viscous
problems. */

event logfile (i++;t<40.){
  stats s;
  if(i%2==0){
    s = statsf (TS);
    fprintf (stderr, "# %g TS %.12f %.9f\n", t, s.min, s.max);
  }
  if(i%2==1){
    s = statsf (TL);
    fprintf (stderr, "# %g TL %.12f %.9f\n", t, s.min, s.max);
  }
  stats s2 = statsf(v_pc.y);
  fprintf (stderr, "##  v_pc.y  %.12f %.12f\n", s2.min, s2.max);
}

event sauv(i+=200,last){
  dump();
}