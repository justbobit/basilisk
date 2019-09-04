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

double geometry(double x, double y, double Radius) {

  coord center_circle, center_rectangle, size_rectangle;
  center_circle.x = 0.5;
  center_circle.y = 0.5;

  center_rectangle.x = -0.0;
  center_rectangle.y =  -0.5;

  size_rectangle.x = 0.51;
  size_rectangle.y = 1.21; 

  double s = circle (x, y, center_circle, Radius);
  double r = -rectangle (x, y, center_rectangle, size_rectangle);

  // double zalesak = -difference (s , r);
  double zalesak = -s;
  // double zalesak = -r;

  return zalesak;
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
  We set the initial velocity field and set a rotating field inside the notched
  disk. */
  foreach_face(){
    if(cs[] != 0. && cs[] != 1.){
      v_pc.x[] = 1.;
     }
     else{
      v_pc.x[] = 0.;
    }
  }
  
  boundary((scalar*){v_pc});


  face vector grad_dist[];
  foreach(){

    if(dist[]>0.){
    foreach_dimension(){
      double a = max(0.,dist[]    - dist[-1,0]);
      double b = min(0.,dist[1,0] - dist[]    );
      grad_dist.x[] = max(a,b);
      }
    }
    else{
      foreach_dimension(){
        double a = min(0.,dist[]    - dist[-1,0]);
        double b = max(0.,dist[1,0] - dist[]    );
       grad_dist.x[] = max(a,b);
      }
    }
  }

  for (scalar f in velocity)
    f.gradient = gradient;

	
  
  advection (velocity, grad_dist, L0/(1 << MAXLEVEL));
  view (fov = 20.0645, quat = {0,0,0,1}, tx = -0.501847, ty = -0.497771, 
  	bg = {1,1,1}, width = 600, height = 600, samples = 1);
  cells();
	draw_vof("cs");
  squares("v_pc.y", min =-1., max = 1.);
  save ("vpcx.png");
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


