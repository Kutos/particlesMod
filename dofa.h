#include <cmath>

double dofa ( double a, double om0 ) 
{  
  	double omai = omofa(a,om0) ;
  	double alpha0 = 6./11. - 15./2057.*log(om0) + 205./540421.*log(om0)*log(om0) ;
  	double fai = fofa(a,om0) ;
	return  a/(om0 + (1.-om0)*a*a*a)*( pow(om0,alpha0) + 3.*om0/2. ) / ( fai + 3.*omai/2.) ;
}
