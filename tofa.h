#include <cmath>

double tofa ( double a, double om0 ) 
{  

  double ai3 = a*a*a;

  return   2./3./sqrt(1.-om0)*log( ( sqrt(om0+(1-om0)*ai3) + sqrt((1-om0)*ai3) )/sqrt(om0) ) / 100.;
}
