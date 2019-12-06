/**
## My set of functions 


A function to rescale normals so that they are unit vectors w.r.t. the
2-norm (by default, the 1-norm is adopted for efficiency purposes). */

coord normal (Point point, scalar c) {
  coord n = mycs (point, c);
  double nn = 0.;
  foreach_dimension()
    nn += sq(n.x);
  nn = sqrt(nn);
  foreach_dimension()
    n.x /= nn;
  return n;
}

/**
A function to compute 2-norm normal in every cell. */

void compute_normal (scalar f, vector normal_vector) {
  foreach() {
    coord n = normal (point, f);
    foreach_dimension() 
      normal_vector.x[] = n.x;
  }
  boundary((scalar*){normal_vector});
}

/** 
#sign function
*/
double sign2 (double x)
{
  return(x > 0. ? 1. : x<0 ? -1. : 0.);
}

void InterpStencil (coord p, int Stencil[]){
  if(p.x<0){
    Stencil[0] = -1;
  } 
  else{
    Stencil[0] = 1;
  } 

  if(p.y<0) {Stencil[1] = -1;}
  else {Stencil[1] = 1;}
}


double capteur(Point point, scalar s){
  return s[];
}

double linearInterpolation(double x1, double f_x1, double x2, double f_x2, double x)
{
  double result = (x - x1)/(x2-x1)*f_x2  + (x2-x)/(x2-x1)*f_x1;
  return result;
}

/**
Bilinear interpolation of a scalar s
*/
double mybilin (Point point, scalar s, int Stencil[], coord p_interp){
  #if dimension == 2

      double x1 = x + Stencil[0]*Delta;
      double x2 = x + (Stencil[0]+1)*Delta;  
      double y1 = y + Stencil[1]*Delta;
      double y2 = y + (Stencil[1]+1)*Delta;

      double f11 = s[Stencil[0],Stencil[1]];
      double f21 = s[Stencil[0]+1,Stencil[1]];
      double f12 = s[Stencil[0],Stencil[1]+1];
      double f22 = s[Stencil[0]+1,Stencil[1]+1];

      double R1 =  linearInterpolation(x1, f11, x2, f21, p_interp.x);
      double R2 =  linearInterpolation(x1, f12, x2, f22, p_interp.x);
      double  P =  linearInterpolation(y1,  R1, y2,  R2, p_interp.y);
      return P;
  #endif

}

void correct_values(Point point, scalar f, int Stencil[], double error){
  f[Stencil[0]  ,Stencil[1]  ] = f[Stencil[0]  ,Stencil[1]  ] + error/32.;
  f[Stencil[0]+1,Stencil[1]  ] = f[Stencil[0]+1,Stencil[1]  ] + error/32.;
  f[Stencil[0]  ,Stencil[1]+1] = f[Stencil[0]  ,Stencil[1]+1] + error/32.;
  f[Stencil[0]+1,Stencil[1]+1] = f[Stencil[0]+1,Stencil[1]+1] + error/32.;
}
