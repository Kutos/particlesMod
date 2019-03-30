#include <random>
#include <math.h>
using namespace std;
void ic_generator(int num, double lbox, double hd, int flag, double *pk, double *positions, double *velocities){
	int numk = num/2+1;
	double kf = 2.*M_PI/lbox;
	double s = lbox/num;
	double psi, qj, qjkf;
	double en[numk], an[numk];

	random_device r1, r2;
	default_random_engine generator1(r1()), generator2(r2());
	uniform_real_distribution<double> distribution1(0.,1.);
	uniform_real_distribution<double> distribution2(0.,1.);

	for(int n=1;n<numk;n++){
		en[n] = 2*M_PI*distribution1(generator1);
	}
	en[0] = 0;
	an[0] = 0;

	switch(flag){
		case 1:
			for(int n=1;n<num;n++){
				an[n] = sqrt( *(pk+n)/kf ) / double(n);
			}
			break;

		default:
			for(int n=1;n<num;n++){
				an[n] = sqrt( - *(pk+n) / kf * log( 1 - distribution2(generator2))) /double(n);
			}
		}

	for(int j=0;j<num;j++){
		qj = *(positions+j);
		qjkf = qj*kf;
		psi = 0;

		for(int n=1; n<numk;n++){
			psi+=an[n] * sin(en[n]+qjkf*double(n));
		}

		*(positions+j) -= 2.*psi;
		*(velocities+j) = -2.-psi*kf;
	}
}
