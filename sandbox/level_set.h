/**
# Level-Set functions

The level set function is defined and initialized elsewhere (typically by the 
user), the face vector field `uf` and the timestep `dt` are defined by a
solver.

This event integrates advection equations of the form
$$
\partial_t\phi+\mathbf{u_f}\cdot\nabla \phi=0
$$
where $\mathbf{u_f}$ is the velocity field and $\phi$ is the level set 
function.

 */
extern scalar * level_set;
extern face vector uf;
extern double dt;


/**
The integration is performed using the Bell-Collela-Glaz scheme. */
#include "bcg.h"
#include "alex_functions.h"
#include "embed.h"


#include "curvature.h"

/**
# Various functions

## Fluxes for functions that do not take into account embedded boundaries

We redefine the fluxes for the level set advection and reconstruction of the
associated velocity field, they do not take into account the embedded boundaries
*/ 


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

/**
## Advection without taking into account embedded boundaries
*/

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
    foreach(){
      foreach_dimension(){
        f[] += p.dt*(flux.x[] - flux.x[1])/(Delta); 
      }
    }
  }
  boundary (p.tracers);

  if (!p.src)free (lsrc);
}

/**
If needed (no tracer.h included), on can use this event for the level_set advection
*/


/**
## Anisotropy of the Gibbs-Thomson relation

This function is to be used only with the embed boundary module. It relies on
the calculation of gradients on both side of a boundary and modifies the value
of v_pc which is the phase change velocity using the Stefan relation.
*/
double fac1(double x, double y){
#if GT_aniso // prefactor for GT-formulation
  double theta = atan2 (y-L0/2., x-L0/2.);
  return (1.-0.5*cos(4*theta));
#else
  return 1.;
#endif
}

/**
# Velocity on the fluid-solid interface
*/
void phase_change_velocity_LS_embed (scalar cs, face vector fs, scalar tr,
 scalar tr2, vector v_pc, scalar dist, double L_H, 
 double NB_width, int nb_cell_NB, double lambda[2],
 double epsK, double epsV) {

  scalar T_eq[];

  scalar curve[];
  curvature (cs,curve);
  boundary({curve});

/**
The temperature on the interface is defined by :
$$
 T_{\Gamma} = T_{m} - \epsilon_{\kappa} * \kappa - \epsilon_{v} v_{pc}
$$
*/
#if Gibbs_Thomson

  foreach(){
    // here we suppose that Tm = 0
    T_eq[] = (epsK*curve[] 
      -epsV*sqrt(v_pc.x[]*v_pc.x[]+v_pc.y[]*v_pc.y[]))*fac1(x,y);
  }

#else
  foreach(){
    T_eq[] = 0;
  }
#endif

  boundary({T_eq});
  restriction({T_eq});

 /**
  The phase change velocity $\mathbf{v}_{pc}$ is

  $$
\mathbf{v}_{pc} = \frac{1}{L_H}\left(\lambda_{L}\left.\nabla T_{L}\right|_{\Gamma} - \lambda_{S}\left.\nabla T_{S}\right|_{\Gamma}\right) 
  $$
  
  here we use the embed_gradient_face_x defined in embed that gives a proper
  definition of the gradients with embedded boundaries. */
  vector gtr[], gtr2[];
  boundary({tr});

  foreach(){
    foreach_dimension(){
      gtr.x[] = 0.;
    }
    if(fabs(dist[])<NB_width){
      if(interfacial(point, cs)){
        coord n       = facet_normal( point, cs ,fs) , p;
        normalize(&n);
        double alpha  = plane_alpha (cs[], n);
        line_length_center (n, alpha, &p);
        double c    = 0.;
        double temp = T_eq[];
        double grad = dirichlet_gradient(point, tr, cs , n, p, 
          temp, &c);
        foreach_dimension(){
          gtr.x[] += grad*n.x+tr[]*c;
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

  boundary({tr2});

  vector p_sauv[], n_sauv[];

  foreach(){
    foreach_dimension(){
      gtr2.x[] = 0.;
      p_sauv.x[] = nodata;
      n_sauv.x[] = nodata;
    }
    if(fabs(dist[])< NB_width){
      if(interfacial(point, cs)){
        coord n       = facet_normal( point, cs ,fs) , p;
        normalize(&n);
        double alpha  = plane_alpha (cs[], n);
        line_length_center (n, alpha, &p);
        foreach_dimension(){
          p_sauv.x[] = p.x;
          n_sauv.x[] = n.x;
        }
        double c=0.;
        double temp = T_eq[];
        // foreach_dimension(){
        //   temp += v_pc.x[]*v_pc.x[];
        // }
        // temp = -0.2*sqrt(temp);
        double grad = dirichlet_gradient(point, tr2, cs , n, p, 
          temp, &c);
        foreach_dimension(){
          gtr2.x[] += grad*n.x+tr2[]*c;
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
  With the the normal vector and the gradients of the tracers we can now 
  compute the phase change velocity $\mathbf{v}_{pc}$. */

  foreach_face(){
    v_pc.x[]  = 0.; 
  }

  foreach_face(){
    if(interfacial(point, cs)){
      v_pc.x[]     =  (lambda[1]*gtr2.x[] - lambda[0]*gtr.x[])/L_H;
    }
  }
  boundary((scalar *){v_pc});
}



event LS_advection (i++,last) {
}
/**
Reinitialization method can be chosen by overloading this function */
event LS_reinitialization(i++,last) {
}

double dirichlet_embed_LS(Point point, scalar cs, face vector fs, 
  face vector v_pc){

  coord n = facet_normal (point, cs, fs);
  normalize(&n);
  double a = 0;
  foreach_dimension()
  a += n.x*v_pc.x[];
  return a;
}

/**
# Reinitialization of the level-set function

Redistancing function with subcell correction, see the work of [Russo et al.,
1999](#russo_remark_2000). 

Let $\phi$ be a function close to a signed function that has been perturbed by
numerical diffusion. By iterating on this equation :
$$
\dfrac{\partial \phi}{\partial t} = sign(\phi^{0}) \left(1- \nabla \phi\right)
$$
we can correct or redistance $\phi$ to make it a signed function.

We discretize far from the interface into :
$$\phi_i^{n+1} = \phi_i^n - \Delta t * sign(\phi_i^0) G(\phi)_i$$
with:
$$
G(\phi_i) = \left\{ \begin{array}{ll}
max ( |a_+|, |b_-|) &, if \phi_i^0 >0 \\
max ( |a_-|, |b_+|) &, if \phi_i^0 <0 \\
\end{array}
\right.
$$
and:
$$
a= D_x^-\phi_i = (\phi_i -\phi_{i-1})/\Delta x\\
b= D_x^+\phi_i = (\phi_{i+1}- \phi_{i})/\Delta x
$$
which is a classical upwind scheme.

Near the interface, *i.e.* for cells where:
$$
\phi^0_i\phi^0_{i+1} \leq 0 \text{ or } \phi^0_i\phi^0_{i-1} \leq 0
$$
the scheme must stay truly upwind, meaning that the movement 0 level-set of the
function must be as small as possible. Therefore the upwind numerical scheme is
modified to :
$$\phi_i^{n+1} = \phi_i^n - \frac{\Delta t}{\Delta x} ( sign(\phi_i^0)
|\phi_i^n| - D_i)$$

with:
$$D_i = \Delta x * \frac{\phi_i^0}{\Delta \phi_i^0}$$

and:  
$$\Delta \phi_0^i = \max((\phi^0_{i-1}-\phi^0_{i+1})/2,\phi^0_{i-1}-\phi^0_{i},
\phi^0_{i}-\phi^0_{i+1})$$  

 */
void LS_reinit2(scalar dist, double dt, double NB, int it_max){
  vector gr_LS[];
  int i ;
  double eps = dt/100., eps2 = eps/2.;

/**
We create `dist0[]` which will be a copy of the initial level-set function
before the iterations and `dist_n[]` which will be $\phi^{n}$ used for the
iterations.
*/
  scalar dist0[], dist_n[];
  foreach(){
    dist0[] = dist[] ;
  }
  boundary({dist0});


/**
Iteration loop
*/
  for (i = 1; i<=it_max ; i++){
    double res=0.;
    foreach(){
      dist_n[] = dist[] ;
    }
    boundary({dist_n});


/**
First step, calculate $\Delta t$ which depends on the local position of the
0-level set. This time step is only used for the cells near the interface, we
tag them using `min_neighb` ($<0$ for interfacial cells).
*/
    double xCFL = 1.;
    foreach(reduction(min:xCFL)){
      double min_neighb = 1.;
      foreach_dimension(){
        min_neighb = min (min_neighb, dist_n[-1,0]*dist_n[]); 
        min_neighb = min (min_neighb, dist_n[ 1,0]*dist_n[]);
      }

/**
Then we calculate the CFL for interfacial cells :

CFL  = \text{min}(Delta x * \phi^0_i)/\Delta \phi_i
*/
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
      if(dist0[]< NB){
/**
We only iterate in the narrow band
*/
        double delt =0.;
        //min_neighb : variable for detection if cell is near
        //             the zero of the level set function

        double min_neighb = 1.;
        foreach_dimension(){
          min_neighb = min (min_neighb, dist_n[-1,0]*dist_n[]);
          min_neighb = min (min_neighb, dist_n[ 1,0]*dist_n[]);
        }

        if(min_neighb < 0.){ // the cell contains the interface
          double dist1= 0., dist2= 0.,dist3= 0.;
          foreach_dimension(){
            dist1 += pow((dist0[1,0]-dist0[-1,0])/2.,2.);
            dist2 += pow((dist0[1,0]-dist0[    ]),2.);
            dist3 += pow((dist0[   ]-dist0[-1,0]),2.);
          }
          double Dij = Delta*dist0[]/
          max(eps2,sqrt(max(dist1,max(dist2,dist3))));
          delt = (sign2(dist0[])*fabs(dist_n[])-Dij)/Delta;
        }
        else{ 
          if(dist0[]>0){ // no interface in the cell
            foreach_dimension(){
              double a = max(0.,(dist_n[]    - dist_n[-1,0])/Delta);
              double b = min(0.,(dist_n[1,0] - dist_n[]    )/Delta);
              delt   += max(pow(a,2.),pow(b,2.));
            }
            delt = sign2(dist0[])*(sqrt(delt) - 1.);
          }
          else{
            foreach_dimension(){
              double a = min(0.,(dist_n[]    - dist_n[-1,0])/Delta);
              double b = max(0.,(dist_n[1,0] - dist_n[]    )/Delta);
              delt   += max(pow(a,2.),pow(b,2.));
            }
            delt = sign2(dist0[])*(sqrt(delt) - 1.);
          }
        }
        if(fabs(dist0[])>0.3*NB){
          dist[] -= 0.5*dt*delt;
        }
        else{
          dist[] -= xCFL*dt*delt;
        }
        if(fabs(delt)>=res) res = fabs(delt);
      }
    }

    boundary({dist});
    restriction({dist});


/**
Iterations are stopped when $max(|\phi_i^{n+1}-\phi_i^n| < eps$
*/
    if(res<eps){
      break;
    }
    if(i==it_max){
      if(res>10.){
        fprintf(stderr,"BAD LS_REINIT %g %g \n", xCFL, res);
        exit(1);
      } 
    }
  }
}
/**
# Reconstruction of a velocity off an interface
*/  
void recons_speed(scalar dist, double dt, int nb_cell_NB,
  double NB, scalar * LS_speed){
/**
The phase change field is only defined in interfacial cells, it must now be
reconstructed on the cells in the direct vicinity. We use here the method
described by [Peng et al., 1999](#peng_pde-based_1999) by solving the following
equation:

$$
\dfrac{\partial v_{pc}}{\partial t}  + S(\phi) \frac{\nabla \phi}{|\nabla \phi|}. \nabla v_{pc}= 0
$$

First part : calculate 
$$ 
S(\phi) \frac{\nabla \phi}{|\nabla \phi|}
$$

using a centered approximation
*/

  vector n_dist[], grad_dist[];

  foreach(){

    double sum=1e-15;
    foreach_dimension(){
      // grad_dist.x[] = (dist[]-dist[-1,0])/(2.*Delta),
      grad_dist.x[] = min(dist[1,0]-dist[],
        min(dist[]-dist[-1,0],(dist[1,0]-dist[-1,0])/2))/Delta;
      sum += grad_dist.x[]*grad_dist.x[];

    }
    foreach_dimension(){
      n_dist.x[] = grad_dist.x[]/sqrt(sum)*sign2(dist[]);
    }
  }
  boundary((scalar *) {n_dist});


/**
Second part do the advection a few times to extend the velocity from the 
surface along the normal to the surface. 
*/ 

  int ii;
  scalar delt_err[];
      
  for (scalar f in LS_speed){
    fprintf(stderr, "COMPOSANTE\n" );
    scalar f_i[];

    foreach()
      f_i[] = f[];
    boundary({f_i});

    for (ii=1; ii<=10*nb_cell_NB; ii++){
      scalar f2[];
      foreach(){
        f2[] = f[];
      }
      boundary({f2});
      restriction({f2});

      foreach(){
        if(fabs(dist[])<0.9*NB){
          foreach_dimension(){
            if(cs[] ==0. || cs[] == 1.){
              f[] -= dt *(max(n_dist.x[],0.)*(f2[   ]-f2[-1,0])/Delta
               +min(n_dist.x[],0.)*(f2[1,0]-f2[    ])/Delta);
            }
          }
        }
      }
      boundary({f});
      restriction({f});

      foreach()
        delt_err[] = nodata;

      foreach(){
        if(interfacial(point, cs)){
          coord n       = facet_normal( point, cs ,fs) , p;
          normalize(&n);
          double alpha  = plane_alpha (cs[], n);
          line_length_center (n, alpha, &p);
          
          int Stencil[2];
          InterpStencil(p, Stencil);

          coord p_interp = {x+p.x*Delta, y + p.y*Delta};

          double f_temp = mybilin(point, f, Stencil, p_interp);
          // double f_temp = bilinear_noembed(point,f);
          double error  = f_i[] - f_temp;
          delt_err[] = error;
          if(ii>(5*nb_cell_NB) && ii%(2*nb_cell_NB) == nb_cell_NB-1){
            correct_values(point, f, Stencil, error);
          }
        }
      }

      boundary({f});
      restriction({f});
    }

    stats s2  = statsf(delt_err);
    fprintf(stderr, "%g %g %g\n", s2.min, s2.max, s2.stddev);

  }
}

/**
## References

~~~bib

@article{russo_remark_2000,
  title = {A remark on computing distance functions},
  volume = {163},
  number = {1},
  journal = {Journal of Computational Physics},
  author = {Russo, Giovanni and Smereka, Peter},
  year = {2000},
  pages = {51--67}
}

@article{peng_pde-based_1999,
  title = {A {PDE}-Based Fast Local Level Set Method},
  volume = {155},
  issn = {0021-9991},
  url = {http://www.sciencedirect.com/science/article/pii/S0021999199963453},
  doi = {10.1006/jcph.1999.6345},
  pages = {410--438},
  number = {2},
  journaltitle = {Journal of Computational Physics},
  author = {Peng, Danping and Merriman, Barry and Osher, Stanley and Zhao, Hongkai and Kang, Myungjoo},
  urldate = {2019-09-09},
  date = {1999-11}
}
~~~

*/