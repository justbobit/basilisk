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
  
  foreach_face() {

    v_pc.x[] = 0.;
  }



  foreach_face(){
    if (fabs(dist[])<=0.5*Delta){

      coord n = facet_normal (point, cs, fs);
      normalize(&n);

      v_pc.x[] = (1.4*gtr.x[]*n.x*fs.x[] 
                  + gtr2.x[]*n.x*(1.-fs.x[]))*0.1;
    }
    else{
      if(fabs(dist[])<0.9*NB_width){
         v_pc.x[] = (1.4*gtr.x[]*fs.x[] + gtr2.x[]*(1.-fs.x[]))
                  * pow(0.5*Delta/dist[],2.); 
      }
    }
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
  DT = 5.e-4;
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
    TL[] = TL_inf;
    TS[] = TS_inf;
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
  if(t>0.5){  
  phase_change_velocity_LS_embed (cs, fs ,TL, TS, v_pc, dist, L_H, NB_width);
  }
  else{
    foreach_face()
      v_pc.x[] = 0.;
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

  int sum = 0;
  foreach(reduction(+:sum)) {
    
    if(cs0[] == 0.){
      cs0[] -= cs[];
      if(fabs(cs0[])> SEPS){ 
        TL[] = T_eq; // 0-order approx. better to use gradients if possible
                     // will be fixed
        sum++;
      }
    }
    else if(cs0[] == 1.){
      cs0[] -= cs[];
      if(fabs(cs0[])> SEPS){ 
        TS[] = T_eq;
        sum++;
      }
    }
    else if(cs[] == 0)  TL[] = T_eq;  
    else if(cs[] == 1.) TS[] = T_eq;
  }
  
  boundary({TL,TS});
  restriction({TL,TS});
  if(sum>0){
    stats s = statsf (v_pc.y);
    fprintf (stderr, "##%d %.12f %.9f %g\n",sum, s.sum, s.min, s.max);
  }

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
    LS_reinit2(dist,L0/(1 << MAXLEVEL), 0.9*NB_width,
      20.*nb_cell_NB);
  }
}


/**
We produce an animation of the tracer field. */

event movies ( i +=24,last)
{
  boundary({TL,TS});
  scalar visu[];
  foreach(){
    visu[] = (1.-cs[])*TL[]+cs[]*TS[] ;
    // visu[] = TL[] ;
  }
  boundary({visu});
  view (fov = 16.642, quat = {0,0,0,1}, tx = -0.0665815, 
    ty = -0.00665815, bg = {1,1,1}, width = 600, height = 600, samples = 1);
  draw_vof("cs");
  squares("visu", min =-1, max = 1);
  save ("visu.mp4");
  clear();
  view (fov = 16.642, quat = {0,0,0,1}, tx = -0.0665815, 
    ty = -0.00665815, bg = {1,1,1}, width = 600, height = 600, samples = 1);
  draw_vof("cs");
  squares("v_pc.y");
  save ("v_pc.mp4");
  clear();
  view (fov = 16.642, quat = {0,0,0,1}, tx = -0.0665815, 
    ty = -0.00665815, bg = {1,1,1}, width = 600, height = 600, samples = 1);
  draw_vof("cs");
  cells();
  squares("TS", min = -1 , max = 1);
  save ("TS.mp4");
}

/**
We check the number of iterations of the Poisson and viscous
problems. */

event logfile (i++;t<0.8){
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

