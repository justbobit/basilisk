#ifndef BASILISK_HEADER_18
#define BASILISK_HEADER_18
#line 1 "./../level_set.h"
/**
# Level-Set advection event

This event integrates advection equations of the form
$$
\partial_tLS+\mathbf{u_f}\cdot\nabla LS=0
$$
where $\mathbf{u_f}$ is the velocity field and $LS$ is the level set function.

The level set function is defined and initialized elsewhere (typically by the 
user), the face vector field `uf` and the timestep `dt` are defined by a
solver. */
extern scalar * level_set;
extern face vector uf;
extern double dt;


/**
The integration is performed using the Bell-Collela-Glaz scheme. */
#include "bcg.h"
#include "alex_functions.h"
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
phi^n+1 = phi^n - Delta t S(phi) G(phi) far from the interface
Near the interface it is modified to :
phi^n+1 = phi^n - Delta t/Delta x ( sgn(phi^0) |phi^n| - Di)

with:
Di = Delta x * phi_i^0/Delta phi_0^i

with:  
$$Delta phi_0^i = max((phi^0_{i-1}-phi^0_{i+1})/2,phi^0_{i-1}-phi^0_{i},
phi^0_{i}-phi^0_{i+1})$$  
  

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
      if(fabs(dist_eps[])<NB){
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
    }
    foreach(reduction(max:res)){
      double delt =0.;
      if(fabs(dist_eps[])<NB){
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
    }
    
    boundary({dist});
    restriction({dist});
  // if((i%10 == 0) & (i<it_max)) 
  //   fprintf(stderr,"# REINIT_LS %d %d %6.2e %6.2e 
  //   %6.2e %f %f\n",i, sum, res,eps, xCFL,dt, NB);

    if(res<eps){
      fprintf(stderr,"# REINIT_LS %d %d %6.2e %6.2e %6.2e %f %f\n",i, sum, res,
        eps, xCFL,dt, NB);
      break;
    }
    // if(i==it_max)fprintf(stderr,"# REINIT NOT CONVERGED %d %6.2e %6.2e \n", sum,
    //   res,eps);
  }
}
#endif
