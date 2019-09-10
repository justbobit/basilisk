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
#define MAXLEVEL 5


double (* gradient) (double, double, double) = NULL;

scalar dist[];
face vector v_pc[];
scalar * level_set = {dist};

scalar * velocity = {v_pc.x,v_pc.y};

int nb_cell_NB;

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

  return zalesak;
}

double sign2 (double x)
{
  return(x > 0. ? 1. : x<0 ? -1. : 0.);
}


int main(){
	int N = 1 << MAXLEVEL;
	init_grid(N);
	nb_cell_NB = 1 << 2 ;


  foreach_vertex() {
    dist[] = -geometry(x,y,L0/4.);
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
      double theta = atan2 (y, x);
      v_pc.x[] = cos(2.*theta);
      v_pc.y[] = sin(2.*theta);
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
*/
  
  foreach(){

    double sum=1e-15;
    foreach_dimension(){
      grad_dist.x[] = (dist[1,0]-dist[-1,0])/(2*Delta);
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
  scalar * components;

  components = {v_pc.x,v_pc.y};

  dt = L0/(1 << MAXLEVEL);
  int ii;
  for (ii=1; ii<=1; ii++){
    for (scalar f in components){
      foreach(){
        foreach_dimension(){
          f[] -= dt *(max(n_dist.x[],0.)*(dist[   ]-dist[-1,0])/Delta
                     +min(n_dist.x[],0.)*(dist[1,0]-dist[    ])/Delta);
        }
      }
      boundary({f});
    }
  }

  view (fov = 20.0645, quat = {0,0,0,1}, tx = -0.501847, ty = -0.497771, 
    bg = {1,1,1}, width = 600, height = 600, samples = 1);
  cells();
  draw_vof("cs");
  squares("v_pc.x", min =-1., max = 1.);
  save ("vpcx.png");
  clear();  
  view (fov = 20.0645, quat = {0,0,0,1}, tx = -0.501847, ty = -0.497771, 
    bg = {1,1,1}, width = 600, height = 600, samples = 1);
  cells();
  draw_vof("cs");
  squares("v_pc.y", min =-1., max = 1.);
  save ("vpcy.png");

/**
	Here we tag the cells where we will reconstruct the speed away from the
	interface
*/
  scalar x_tag[],x_tag1[];

  foreach(){
  	if(cs[] != 0. && cs[] != 1.){
  			x_tag[]  = 1.;
  			x_tag1[] = 1.;
  	}
  	else{
  			x_tag[]  = 0.;
  			x_tag1[] = 0.;
  	}
  }
  boundary({x_tag,x_tag1});

  int j;
  for (j=1;j<=nb_cell_NB;j++){
  		foreach(){
  		if(x_tag[] == 1.){
  			foreach_neighbor(1)
  				x_tag1[] = 1.;
  		}
    }
    boundary({x_tag1});

    foreach(){
    	x_tag[] = x_tag1[];
    }
    boundary({x_tag});

  }
  view (fov = 20.0645, quat = {0,0,0,1}, tx = -0.501847, ty = -0.497771, 
  	bg = {1,1,1}, width = 600, height = 600, samples = 1);
  cells();
	draw_vof("cs");
  squares("x_tag", min =0., max = 1.);
  save ("x_tag.png");
  dump();  


}

/**
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
