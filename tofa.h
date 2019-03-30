#include <cmath>
double tofa(double a, double om0){ return 2./3./sqrt(1.-om0)*log((sqrt(om0+(1.-om0)*a*a*a) + sqrt((1.-om0)*a*a*a))/sqrt(om0))/100.;}
