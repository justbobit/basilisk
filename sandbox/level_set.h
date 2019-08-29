/**
# Level-Set functions

The level set function is defined and initialized elsewhere (typically by the 
user), the face vector field `uf` and the timestep `dt` are defined by a
solver.

This event integrates advection equations of the form
$$
\partial_tLS+\mathbf{u_f}\cdot\nabla LS=0
$$
where $\mathbf{u_f}$ is the velocity field and $LS$ is the level set function.

 */
extern scalar * level_set;
extern face vector uf;
extern double dt;


/**
The integration is performed using the Bell-Collela-Glaz scheme. */
#include "bcg.h"
#include "alex_functions.h"


/**
This function is to be used only with the embed boundary module. It relies on
the calculation of gradients on both side of a boundary and modifies the value
of v_pc which is the phase change velocity using the Stefan relation.
*/

void phase_change_velocity_LS_embed (scalar cs, face vector fs, scalar tr,
 scalar tr2, face vector v_pc, scalar dist, double L_H, 
 double NB_width, int nb_cell_NB, double lambda[2]) {
  
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
        double grad = dirichlet_gradient(point, tr, cs , n, p, 
          temp, &c);
        foreach_dimension(){
          gtr.x[] += grad*n.x;
        }
      }
    }
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
        double grad = dirichlet_gradient(point, tr2, cs , n, p, 
          temp, &c);
        foreach_dimension(){
          gtr2.x[] += grad*n.x;
        }
      }
    }
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


  foreach_face(){
    if ( (  fs.y[] !=0. && fs.y[] != 1.   ) 
              && fabs(dist[])<=0.9*NB_width){
      v_pc.x[]     =  (lambda[1]*gtr2.x[] - lambda[0]*gtr.x[])/L_H;
      v_pc2.x[]    = v_pc.x[];
    }
  }

  boundary((scalar *){v_pc});
  boundary((scalar *){v_pc2});

  /** We now propagate the velocity obtained in the interfacial cells away from
  the interface. The level set function obtained after propagation is not a
  distance function, it only eases the reinit process of the level set.

  */

  int ii;
  for (ii=1; ii<=nb_cell_NB; ii++){
    foreach(){
      // if(fabs(dist[])<0.9*NB_width){

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
        v_pc2.x[] = (v_pc.x[sig[0]*k1,sig[1]*k2])
                    *exp( -((fabs(dist[]-grad_dist.x*Delta/2.))/
                      (2.*nb_cell_NB*Delta)));
        v_pc2.y[] = v_pc.y[sig[0]*k1,sig[1]*k2]
                    *exp( -((fabs(dist[]-grad_dist.y*Delta/2.))/
                      (2.*nb_cell_NB*Delta)));
      // }
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
If needed (no tracer.h included), on can use this event for the level_set advection
*/

event LS_advection (i++,last) {
}
/**
Reinitialization method can be chosen by overloading this function */
event LS_reinitialization(i++,last) {
}
/**
# LS_reinit function  

V2 of the reinit function with subcell correction.
Based on the work of Russo2000.
$$\phi^n+1 = \phi^n - \Delta t S(\phi) G(\phi)$$
far from the interface.

Near the interface it is modified to :
$$\phi^n+1 = phi^n - \Delta t/\Delta x ( sgn(\phi^0) |\phi^n| - Di)$$

with:
$$Di = \Delta x * \phi_i^0/\Delta \phi_0^i$$.

with:  
$$\Delta \phi_0^i = max((\phi^0_{i-1}-\phi^0_{i+1})/2,\phi^0_{i-1}-\phi^0_{i},
\phi^0_{i}-\phi^0_{i+1})$$  
  

Based on the work of Russo (2000)
 */
 double dirichlet_embed_LS(Point point, scalar cs, face vector fs, 
    face vector v_pc){

  coord n = facet_normal (point, cs, fs);
  normalize(&n);
  double a = 0;
  foreach_dimension()
    a += n.x*v_pc.x[];
  return a;
}

void LS_reinit2(scalar dist, double dt, double NB, int it_max){
  vector gr_LS[];
  int i ;
  double eps = dt/100., eps2 = eps/2.;
  scalar dist0[], dist_eps[];
  foreach()
    dist0[] = dist[] ;
  boundary({dist0});

  for (i = 1; i<=it_max ; i++){
    double res=0.;
    foreach(){
      dist_eps[] = dist[] ;
    }
    boundary({dist_eps});
    
    double xCFL = 1.;
// 1) we make a copy of dist before iterating on it
// 2) we determine xCFL according to the local size
    int sum = 0;
    foreach(reduction(min:xCFL) reduction(+:sum)){
        sum ++;
        //min_neighb : variable for detection if cell is near
        //             the zero of the level set function

        double min_neighb = 1.;
        foreach_dimension(){
          min_neighb = min (min_neighb, dist_eps[-1,0]*dist_eps[]); 
          min_neighb = min (min_neighb, dist_eps[ 1,0]*dist_eps[]);
        }

        if(min_neighb < 0.){
          double dist1= 0., dist2= 0.,dist3= 0.;
          foreach_dimension(){
            dist1 += pow((dist0[1,0]-dist0[-1,0])/2.,2.);
            dist2 += pow((dist0[1,0]-dist0[    ]),2.);
            dist3 += pow((dist0[   ]-dist0[-1,0]),2.);
          }
          double Dij = Delta*dist0[]/
                  max(eps2,sqrt(max(dist1,max(dist2,dist3))));
  // stability condition near the interface is modified
          xCFL = min(xCFL,fabs(Dij)/(Delta));
        }
    }
    foreach(reduction(max:res)){
      double delt =0.;
        //min_neighb : variable for detection if cell is near
        //             the zero of the level set function

        double min_neighb = 1.;
        foreach_dimension(){
          min_neighb = min (min_neighb, dist_eps[-1,0]*dist_eps[]);
          min_neighb = min (min_neighb, dist_eps[ 1,0]*dist_eps[]);
        }

        if(min_neighb < 0.){
          double dist1= 0., dist2= 0.,dist3= 0.;
          foreach_dimension(){
            dist1 += pow((dist0[1,0]-dist0[-1,0])/2.,2.);
            dist2 += pow((dist0[1,0]-dist0[    ]),2.);
            dist3 += pow((dist0[   ]-dist0[-1,0]),2.);
          }
          double Dij = Delta*dist0[]/
                  max(eps2,sqrt(max(dist1,max(dist2,dist3))));
          delt = (sign2(dist0[])*fabs(dist_eps[])-Dij)/Delta;
        }
        else 
          if(dist0[]>0){
          foreach_dimension(){
            double a = max(0.,(dist_eps[]    - dist_eps[-1,0])/Delta);
            double b = min(0.,(dist_eps[1,0] - dist_eps[]    )/Delta);
            delt   += max(pow(a,2.),pow(b,2.));
          }
          delt = sign2(dist0[])*(sqrt(delt) - 1.);
        }
        else{
          foreach_dimension(){
            double a = min(0.,(dist_eps[]    - dist_eps[-1,0])/Delta);
            double b = max(0.,(dist_eps[1,0] - dist_eps[]    )/Delta);
            delt   += max(pow(a,2.),pow(b,2.));
           }
           delt = sign2(dist0[])*(sqrt(delt) - 1.);
        }
        dist[] -= xCFL*dt*delt;
        if(fabs(delt)>=res) res = fabs(delt);
    }
    
    boundary({dist});
    restriction({dist});
  // if((i%10 == 0) & (i<it_max)) 
  //   fprintf(stderr,"# REINIT_LS %d %d %6.2e %6.2e 
  //   %6.2e %f %f\n",i, sum, res,eps, xCFL,dt, NB);

    if(res<eps){
      // fprintf(stderr,"# REINIT_LS %d %d %6.2e %6.2e %6.2e %f %f\n",i, sum, res,
        // eps, xCFL,dt, NB);
      break;
    }
    // if(i==it_max)fprintf(stderr,"# REINIT NOT CONVERGED %d %6.2e %6.2e \n", sum,
    //   res,eps);
  }
}

