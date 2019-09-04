#include "embed.h"
#include "basic_geom.h"
#include "view.h"


scalar dist[];
face vector v_pc[];
scalar * level_set = {dist};

double  NB_width = 1 ;              // length of the NB
#define MAXLEVEL 7
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
	nb_cell_NB = 1 << 3 ;


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
      v_pc.x[] = 0.;
    }
  }
  
  boundary((scalar*){v_pc});

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

  int j;
  for (j=1;j<=nb_cell_NB;j++){
  		foreach(){
  		if(x_tag[] == 1.){
  			foreach_neighbor(1)
  				x_tag1[] = 1.;
  		}
    }

    foreach(){
    	x_tag[] = x_tag1[];
    }
  }
  cells();
	draw_vof("cs");
  squares("x_tag", min =0., max = 1.);
  save ("x_tag.png");


// /**
// Here we try define method to reconstruct a speed calculated in interfacial cells
// away from the interface
// */


//     foreach(){
//       // if(fabs(dist[])<0.9*NB_width){

//         coord grad_dist;
//         if(dist[]>0.){
//         foreach_dimension(){
//           double a = max(0.,dist[]    - dist[-1,0]);
//           double b = min(0.,dist[1,0] - dist[]    );
//           grad_dist.x = max(a,b);
//           }
//         }
//         else{
//           foreach_dimension(){
//             double a = min(0.,dist[]    - dist[-1,0]);
//             double b = max(0.,dist[1,0] - dist[]    );
//            grad_dist.x = max(a,b);
//           }
//         }
//         int k1, k2;
//         int sig[2] = {1 , 1 };
//         if(dist[]>0.){
//           if(grad_dist.x>0) sig[0] = -1;
//           if(grad_dist.y>0) sig[1] = -1;
//         }
//         else{
//           if(grad_dist.x<0) sig[0] = -1;
//           if(grad_dist.y<0) sig[1] = -1;
//         }
        
//         double ratio = fabs(grad_dist.x)/(max(SEPS,fabs(grad_dist.y))) ;
//         if(ratio > sqrt(3.)){
//           k1 = 1;
//           k2 = 0;
//         }
//         else{
//           if(ratio > 1./sqrt(3.)){
//             k1 = 1;
//             k2 = 1;   
//           }
//         else{
//           k1 = 0;
//           k2 = 1; 
//           }
//         }
//         v_pc2.x[] = (v_pc.x[sig[0]*k1,sig[1]*k2])
//                     *exp( -((fabs(dist[]-grad_dist.x*Delta/2.))/
//                       (2.*nb_cell_NB*Delta)));
//         v_pc2.y[] = v_pc.y[sig[0]*k1,sig[1]*k2]
//                     *exp( -((fabs(dist[]-grad_dist.y*Delta/2.))/
//                       (2.*nb_cell_NB*Delta)));
//       // }
//     }
//     boundary((scalar *){v_pc2});
//     restriction((scalar *){v_pc2});    


//     foreach_face(){
//       if(v_pc.x[]==0. && v_pc2.x[] !=0.) v_pc.x[] = v_pc2.x[];
//     }

//     boundary((scalar *){v_pc});
//     restriction((scalar *){v_pc});
  



}