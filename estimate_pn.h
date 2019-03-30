#include <cmath>
void estimate_pn(double *x, int num, double lbox, double *kn, double *pn){
	int numk = num/2 + 1;
	double real, imag, kf = 2*M_PI/lbox;

	for(int i=0;i<numk;i++){
		real = 0;
		imag = 0;
		*(kn+i) = i*kf;
		for(int j=0;j<num;j++){
			real += cos(kn[i]*x[j]);
			imag += sin(kn[i]*x[j]);
		}
		pn[i] = (real*real + imag*imag)/(double(num)*double(num));
	}
}
