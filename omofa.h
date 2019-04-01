#include <cmath>

double omofa ( double a, double om0 ) 
{  
  return om0/(om0 + (1.-om0)*a*a*a) ;
}
