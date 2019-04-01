#include <random>
using namespace std ;

void ic_generator (int num, double lbox, double hf, int flag, double * pk, double * positions, double * velocities) 
{ 
  int numk = num/2 + 1 ;
  double pi = 3.1415927 ;
  double kf = 2.*pi/lbox ;
  double s = lbox/num ;
  double psi, qj, qjkf ; 
  double en[numk], an[numk] ;
  
  //INITIALIZE RANDOMS:
  //int t1 = mm_systime() ;
  random_device r1, r2;
  default_random_engine generator1(r1()), generator2(r2());
  uniform_real_distribution<double> distribution1(0., 1.);
  uniform_real_distribution<double> distribution2(0., 1.);
  //normal_distribution<double> distribution2(0.,sigmav/10.);
 

  //GENERATE THE RANDOM PHASES
  for (int n=1 ; n < numk ; n++) 
        en[n] =  2.*pi*distribution1(generator1) ;
  //distribution1.reset();

  //PUT AT ZERO FIRST ELEMENT
  en[0] = 0 ;
  an[0] = 0 ;
  //cout << flag << endl ;
//CHECK IF ONLY RANDOM PHASES IS REQUIRED:
switch (flag)
   {
   case 1:
   for (int n=1 ; n < numk ; n++) 
        an[n] =  sqrt( *(pk+n)/kf ) / double(n) ;   
   break;
   default :
   for (int n=1 ; n < numk ; n++) 
        an[n] =  sqrt( - *(pk+n)/kf*log( 1 - distribution2(generator2) ) ) / double(n) ;   
   }

//distribution2.reset();

//COMPUTE THE DISPLACEMENT:
for (int j=0 ; j < num ; j++) 
    { 
      qj = *(positions+j);
      qjkf = qj*kf ;
      psi = 0;

      for (int n=1 ; n < numk ; n++) 
        {
	  psi+= an[n] * sin( en[n] + qjkf*double(n) ) ;
        }

      *(positions+j)-= 2.*psi ;
      *(velocities+j)= - 2.*psi*hf ;
    } 

 
}
