/**
# Reconstruction of a field off an interface


This is a test case for the reconstruction of a field only defined on an
interface off of it. The following method is employed in the level_set.h library
of function in my sandbox.
*/

#include "run.h"
#include "timestep.h"
#include "bcg.h"

#include "embed.h"
#include "basic_geom.h"
#include "view.h"
#define MAXLEVEL 7


scalar dist[];
vector v_pc[];
face vector v_pc_f[];
scalar * level_set = {dist};

scalar * velocity = {v_pc.x,v_pc.y};

int nb_cell_NB;
double  NB_width ;              // length of the NB

/**
The initial interface defined by the level set function is a circle.
*/

double geometry(double x, double y, double Radius) {

  coord center_circle, center_rectangle, size_rectangle;
  center_circle.x = 0.5;
  center_circle.y = 0.5;

  center_rectangle.x = -0.0;
  center_rectangle.y =  -0.5;

  size_rectangle.x = 0.51;
  size_rectangle.y = 1.21; 

  double s = circle (x, y, center_circle, Radius);

  double zalesak = -s;
  // double zalesak = y-0.45;

  return zalesak;
}

double sign2 (double x)
{
  return(x > 0. ? 1. : x<0 ? -1. : 0.);
}

void LS_reinit2(scalar dist, double dt, double NB, int it_max){
  vector gr_LS[];
  int i ;
  double eps = dt/100., eps2 = eps/2.;
  scalar dist0[], d2[], dist_eps[];
  foreach(){
    dist0[] = dist[] ;
    foreach_dimension(){
      d2[] = min(fabs(dist[]),min(fabs(dist[-1,0]),fabs(dist[1,0])));
    }
    d2[] = min(d2[], min(fabs(dist[1,1]),dist[-1,-1]));
  }
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
      if(d2[]< NB){
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
        if(fabs(dist0[])>0.3*NB_width){
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
  // if((i%10 == 0) & (i<it_max)) 
  //   fprintf(stderr,"# REINIT_LS %d %d %6.2e %6.2e 
  //   %6.2e %f %f\n",i, sum, res,eps, xCFL,dt, NB);

    if(res<eps){
      fprintf(stderr,"%d %6.2e %6.2e \n", 
    	xCFL,sum, res,eps);
      if(res>10.) exit(1);
      break;
    }
    if(i==it_max)fprintf(stderr,"%d %6.2e %6.2e \n", 
    	xCFL,sum, res,eps);
      if(res>10.) exit(1);
  }
}



int main(){
	int N = 1 << MAXLEVEL;
	init_grid(N);
	nb_cell_NB = 1 << 3 ;
	NB_width = nb_cell_NB * L0 / (1 << MAXLEVEL);


  foreach_vertex() {
    dist[] = clamp(-geometry(x,y,L0/5.),-1.5*NB_width,1.5*NB_width);
  }
  int iloop;


  for(iloop=1; iloop<=1000*MAXLEVEL;iloop++){

  fprintf(stderr, "# %d\n", iloop);

  if(iloop%100 ==0){
	 	foreach_vertex() {
	    dist[] = clamp(dist[],-(nb_cell_NB+2) * L0 / (1 <<MAXLEVEL),
	    	(nb_cell_NB+2) * L0 / (1 << MAXLEVEL));
	  }	
  }

  boundary ({dist});
  restriction({dist});
  fractions (dist, cs, fs);

  boundary({cs,fs});
  restriction({cs,fs});

  

  /**
  We set the initial velocity phase change field : 
  $$
  v_{pc} =  { cos(2\theta), sin(2\theta)}
  $$
   */

  foreach(){
    if(cs[] != 0. && cs[] != 1.){
      double theta = atan2 (y-L0/2., x-L0/2.);
      // v_pc.x[] = 0.01;
      // v_pc.x[] = 0.7*L0/(1<< MAXLEVEL)+0.7*L0/(1<< MAXLEVEL)*cos(2.*theta);
      v_pc.x[] = 0.4*L0/(1<< MAXLEVEL);
      v_pc.y[] = 0.;
      // v_pc.y[] = 0.7*L0/(1<< MAXLEVEL)*sin(4.*theta);
     }
     else{
      v_pc.x[] = 0.;
      v_pc.y[] = 0.;
    }
  }
  
  boundary((scalar*){v_pc});

  vector n_dist[], grad_dist[];

/**
The phase change field is only defined in interfacial cells, it must now be
reconstructed on the cells in the direct vicinity. We use here the method
described by [Peng et al., 1999](#peng_pde-based_1999) by solving the following
equation:

$$
\partial_t v_{pc}  + S(\phi) \frac{\nabla \phi}{|\nabla \phi|}. \nabla v_{pc}= 0
$$

First part : calculate 
$$ 
S(\phi) \frac{\nabla \phi}{|\nabla \phi|}
$$

using a centered approximation
*/
  
  foreach(){

    double sum=1e-15;
    foreach_dimension(){
      grad_dist.x[] = min(dist[1,0]-dist[],
      			min(dist[]-dist[-1,0],(dist[1,0]-dist[-1,0])/2))/Delta;
      sum += grad_dist.x[]*grad_dist.x[];

    }
    if(sum==0.){
        fprintf(stderr, "%g %g %g %g %g %g\n", dist[1,0], dist[-1,0], 
          dist[0,1], dist[0,-1], x , y);
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

  dt = L0/(1 << MAXLEVEL);

  int ii;
  for (ii=1; ii<=10*nb_cell_NB; ii++){
    for (scalar f in velocity){
    	scalar f2[];
    	foreach(){
    		f2[] = f[];
    	}
    	boundary({f2});
    	restriction({f2});

      foreach(){
      	if(fabs(dist[])<0.9*NB_width){
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
    }
  }

  foreach(){
  	for(scalar f in velocity){
  		f[] = f[]*exp(2.*((2.*NB_width-fabs(dist[]))/(2.*NB_width)-1.));
  	}
  }
  boundary((scalar *){v_pc});
  restriction((scalar *){v_pc});
/**
Advection with the reconstructed speed
*/   
	foreach_dimension(){
    scalar dvpc[];
    foreach(){
      dvpc[] = dist[] >0. ? (v_pc.x[]-v_pc.x[-1,0])/Delta : 
                (v_pc.x[1,0]-v_pc.x[])/Delta;
    }
    boundary({dvpc});
    restriction({dvpc});


    foreach(){
    		v_pc_f.x[] = 0.5*(v_pc.x[-1,0] + v_pc.x[] +(dvpc[-1,0]-dvpc[])*Delta/2.);
    }
  }
  boundary((scalar * ){v_pc_f});
  restriction((scalar *){v_pc_f});


    advection(level_set, v_pc_f, L0 / (1 << MAXLEVEL));

  LS_reinit2(dist,0.49*L0/(1 << MAXLEVEL), 
  	(nb_cell_NB+4) * L0 / (1 <<MAXLEVEL), // n'a pas d'influence
      15);
  if( iloop%60==0){
  	view (fov = 14.6629, quat = {0,0,0,1}, tx = -0.628336, ty = -0.470744, 
  		bg = {1,1,1}, width = 600, height = 600, samples = 1);
	  draw_vof("cs");
	  squares("dist", min =-NB_width, max = NB_width);
	  save ("dist.mp4");
  clear();  

	  draw_vof("cs");
	  squares("v_pc.x", min =-0.2*L0/(1 << MAXLEVEL), 
	  	max = 0.2*L0/(1 << MAXLEVEL));
	  save ("vpcx.mp4");
  }
	
	}
	
		dump();
}



/**
![Animation of the level set function](recons_speed/dist.mp4)(loop)

## References

~~~bib

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
