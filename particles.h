#include <math.h>
#include <fstream>
#include <iostream>
#include <cstring>
#include "estimate_pn.h"
#include "ic_generator.h"

using namespace std;

class particles{
	private :

// nombre de pas effectué
	long step;
// retourne le signe
	int Sg(double a){return 1-(a<0)*2;}
// Enregistre un fichier en Ascii //
// N la taille du tableau //
// enabledRef = 1 : on a donné un tableau ref que l'on veut en première colonne du fichier //
// enabledRef = 2 : on veut avoir le numero de la ligne en première colone du fichier //
// enabledRef = ? : on veut avoir qu'une colone remplite avec le tableau data //
	void saveAsciiFile(string filename, double *data, int N, int enabledRef, double *ref){
		filename += ".txt";
		ofstream outFile;
		outFile.open(filename.c_str());
		switch(enabledRef){
			case 1:
				for(int i=0;i<N;i++){outFile << *(ref+i) << " " << *(data+i) << endl;}
				break;
			case 2:
				for(int i=0;i<N;i++){outFile << i << " " << *(data+i) << endl;}
				break;
			default:
				for(int i=0;i<N;i++){outFile << *(data+i) << endl;}
		}
		outFile.close();
	}
	void saveBinaryFile(string filename, double *data, int N){
		filename += ".data";
		std::ofstream out(filename.c_str(), std::ios::out | std::ios::binary);
		out.write((char *) data, N*sizeof(double));
		out.close();
	}

	public :
	long num;
	double *x;
	double *v;
	double *gama;
// Constructeurs
	particles(long numo){
		step = 0;
		num = numo;
		x = new double [num];
		v = new double [num];
		gama = new double [num];

		std::memset(x, 0, num*sizeof(double));
		std::memset(v, 0, num*sizeof(double));
		std::memset(gama, 0, num*sizeof(double));
	}

	particles(long numo, double *xp){
		step = 0;
		num = numo;
		x = new double [num];
		v = new double [num];
		gama = new double [num];

		for(int i=0;i<num;i++){
			x[i] = *(xp++);
		}
		std::memset(v,0, num*sizeof(double));
		std::memset(gama, 0, num*sizeof(double));
	}

// Permet de répartir les num particules uniformément dans un boite de longueur lbox
	void uniformDistrib(double lbox){
		for(int i=0;i<num;i++){
			*(x+i) = i*lbox/num;
		}
	}
// Permet d'assigner des valeurs à x et v à partir d'une distribution pk linéaire théorique
	void setInitialValuesWithPk(double lbox, double H, double eofa, double fofa, double coeff){
		int numk = num/2 + 1;
		double pk[numk];

		string ref = "./PARTICLES/initial_pk_linear_theo.data";

		std::ifstream in(ref.c_str(), std::ios::in | std::ios::binary);
		in.read((char *) pk, (numk)*sizeof(double));

		in.close();

		ic_generator(num, lbox, H*eofa*fofa, 1, pk, x, v);
		calGama(coeff,lbox,H,eofa);
	}
// Permet de retrouver les vitesses initiales à partir de positions random des particules
	void findVInit(double lbox, double H, double eoa, double foa){
		for(int i=0;i<num;i++){*(v+i)=H*eoa*foa*( *(x+i) - i*lbox/num );}
	}
// Permet de sauvegarder en asci, flag =1 ou en binary les pn
	void savePk(double lbox, double dofai2, double dofa02, int flag){
		int numk = num/2 + 1;
		double pk[numk], kn[numk];
		estimate_pn(x,num,lbox,kn,pk);
		for(int i=0;i<numk;i++){ 2*M_PI/lbox * dofai2 /dofa02 *(*pk+i);}
		string filename = "./DATA/pk_step_" + to_string(step);
		switch(flag){
			case 1:
				saveAsciiFile(filename, pk, numk, 1, kn);
				break;
			default:
				saveBinaryFile(filename,pk, numk);
		}
	}
// Permet de sauvegarder en asci les particules de l'espace des phases
	string savePhaseSpace(){
		string filename = "./DATA/phase_space_step_" + to_string(step);
		saveAsciiFile(filename, v, num, 1, x);
		return filename;
	}
// Permet de plot l'espace des phase
	void plotPhaseSpace(){
		string filename;
		filename = savePhaseSpace();

		FILE *gnuplot = popen("gnuplot -persistent","w");

		fprintf(gnuplot,"set terminal postscript color\n");
		fprintf(gnuplot,"set output \"| ps2pdf - ./DATA/phase_space_%i.pdf'\"\n",step);

		fprintf(gnuplot,"set title \"Phase space at %d step\n",step);
		fprintf(gnuplot,"set xlabel \"x\"\n");
		fprintf(gnuplot,"set ylabel \"v\"\n");
    	fprintf(gnuplot,"plot '%s' u 1:2 \n",filename);
		pclose(gnuplot);
	}
// Permet de faire bouger les particules d'un pas de temps
	void move(double dt, double lbox){
		for(int i=0;i<num;i++){
			*(x+i) += *(v+i)*dt;
			*(v+i) += *(gama+i)*dt;

			*(x+i) += (( *(x+i) < 0 ) - ( *(x+i) > lbox ))*lbox;
		}
		step += 1;
	}

// Calcul le champs d'accélération des particules
	void calGama(double coeff, double lbox, double H, double eofa){
		double dx;
		std::memset(gama, 0, num*sizeof(double));
		for(int i=0;i<num;i++){
			for(int j=0;j<num;j++){
				dx = *(x+j) - (*(x+i));
				dx -= (dx>lbox/2)*lbox;
				*(gama+i) += Sg(dx)-2*dx/lbox;
			}
			*(gama+i) = (*(gama+i))*coeff - 2*H*eofa*(*v+i);
		}
	}
};
