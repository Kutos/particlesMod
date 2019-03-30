#include <cstring>
#include <math.h>
#include "particles.h"
#include "tofa.h"
#include "omofa.h"
#include "eofa.h"
#include "fofa.h"
#include "dofa.h"

// argc : nombre de argument donnés lors de l'appelle du programme
// argv[0] : final (nom de l'executable)
// argv[1] : nombre de particules
// argv[2] : nombre d'étape à effectuer dans la simulation
// argv[3] : taille de la boite

int main ( int argc, char *argv[]){
// num : nombre de particules
// nstep : le nombre de pas à effectuer dans la simulation
// lbox : la taille de la boite
	int num, nstep;
	double lbox;

	switch(argc){
		case 2:
			num = atof(argv[1]);
			nstep = 1000;
			lbox = 1000.;
			break;
		case 3:
			num = atof(argv[1]);
			nstep = atof(argv[2]);
			lbox = 1000.;
			break;
		case 4:
			num = atof(argv[1]);
			nstep = atof(argv[2]);
			lbox = atof(argv[3]);
			break;
		default:
			num = 2048;
			nstep = 1000;
			lbox = 1000.;
	}

// Initialisation des constantes cosmologiques
	double H=100., a=0.1, om0=0.32, G=4.30035e-9, pc0=2.7757e11;

// Réglage du redshit initial et final de la simulation
	double smin=sqrt(a),smax=sqrt(1./(1.+0.5)),ds;
	ds = (smax-smin)/nstep;

// On instancie quelques variable utile dans nos calculs
	double coeff, eoa, foa, t1, t2, s;
// On calculs leurs valeurs initiales
	eoa = eofa(a,om0);
	foa = fofa(a,om0);
	s=smin;
	t1 = tofa(a,om0);
	coeff = 2*M_PI*G*om0*pc0*lbox/(num*a*a*a);

// On créer la variable de classe particles qui contient la majorité de l'algorithme
	particles sysParts(num);

// On fait calculé les valeurs initiales dans l'espace des phases à la classe particles
	sysParts.setInitialValuesWithPk(lbox, H, eoa, foa, coeff);
	for(int i=1;i<nstep+1;i++){
		s += ds;
		a = s*s;

		eoa = eofa(a,om0);
		t2 = tofa(a,om0);
		coeff = 2*M_PI*G*om0*pc0*lbox/(num*a*a*a);

		sysParts.move(t2-t1,lbox);
		sysParts.calGama(coeff, lbox,H, eoa);

		t1 = t2;
	}

	sysParts.savePhaseSpace();
	sysParts.savePn(lbox,1);

	return 0;
}
