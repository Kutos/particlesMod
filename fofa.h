#include <cmath>

double fofa ( double a, double om0 ) 
{  

  double omai = omofa(a, om0) ;
double alpha = 6./11. - 15./2057.*log(omai) + 205./540421.*log(omai)*log(omai) ;

  return  pow(omai,alpha) ;
}
