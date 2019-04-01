void estimate_pn (double * positions, int num, double lbox, double * pn) 
{ 
  int numk = num/2 + 1 ;
  double real, imag, x , kn;
  double pi = 3.1415927 ;
  double kf = 2.*pi/lbox ;
  double oneovernp2 = 1./double(num)/double(num) ;

 for (int n=1 ; n < numk ; n++) 
    { 
      kn = n*kf;
      real = 0 ;
      imag = 0 ;
      for (int i=0 ; i < num ; i++) 
	{
          x = *(positions+i);
	  real+= cos(kn*x);
          imag+= sin(kn*x);  
        }
      *(pn+n)=(real*real + imag*imag)*oneovernp2 ;
    }
  *pn = 0 ; 
}

