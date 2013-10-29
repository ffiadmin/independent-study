#include <iostream>
#include <complex>
#include <string>
#include <fstream>

#define _USE_MATH_DEFINES
#include <math.h>

using std::ifstream;
using std::ofstream;
using std::cout;
using std::cin;
using std::string;

typedef std::complex<double> Complex;

enum {MAXBAS=500, MAXSTP=500, MAXBS3=300, MAXNUC=20};

//many global variables discarded. will look to declare them as they are used

int main(){

//=============================================================================
// Initilization (main.f lines 55 - 98)
//=============================================================================

	ifstream fin("input.txt");
	ofstream fout("output.txt");
	string title; 
	int nbasis;
	fin >> title >> nbasis; 
	fout << title << "\nThere are " << nbasis << " basis functions."; //line 59

	//do loop at line 63:
	int *l,*m,*n;
	l = new int[nbasis]; 
	m = new int[nbasis]; 
	n = new int[nbasis];
	
	double *x,*y,*z,*alpha;
	x = new double[nbasis]; 
	y = new double[nbasis]; 
	z = new double[nbasis]; 
	alpha = new double[nbasis];

	for (int i=1; i <= nbasis; i ++){
		fin >> l[i] >> m[i] >> n[i] >> x[i] >> y[i] >> z[i] >> alpha[i];
		x[i] /= 0.529177249;
		y[i] /= 0.529177249;
		z[i] /= 0.529177249;
		fout << l[i] << ' ' << m[i] << ' ' << n[i] << ' ' << x[i] << ' ' << y[i] << ' ' << z[i] << ' ' << alpha[i];
	}

	int nnuc; fin >> nnuc;
	fout << "\nThere are " << nnuc << " nuclei."; // line 72

	//do loop at line 75:
	double *nx, *ny, *nz, *nchg;
	nx = new double[nnuc];
	ny = new double[nnuc];
	nz = new double[nnuc];
	nchg = new double[nnuc];

	for (int i=1; i<=nnuc; i++){
		fin >> nx[i] >> ny[i] >> nz[i] >> nchg[i];
		nx[i] /= 0.529177249;
		ny[i] /= 0.529177249;
		nz[i] /= 0.529177249;
		fout << nx[i] << ' ' << ny[i] << ' ' << nz[i] << ' ' << nchg[i];
	}

	int nthet, nalp;
	double alpstrt, alpstp, thstart, thstep, min, max;
	fin >> nthet >> nalp >> alpstrt >> alpstp >> thstart >> thstep >> min >> max;

	fout << "\nThe initial values for r and theta are " << alpstrt << " and " << thstart;
	fout << "\nThe program will take " << nthet << " steps.";
	fout << "\nEach step in theta will be " << thstep << " radians.";
	fout << "\nThe program will take " << nalp << " steps.";
	fout << "\nEach step in r will be " << alpstp << " au";
	fout << "\nThe range for the real part of the energy is from " << min << " to " << max;

//=============================================================================
// Calculate normalization constants (main.f lines 98-183
//=============================================================================
	
	//do loop at line 100:
	double t1, t2, t3, t4;
	for (int i=0; i<=nbasis; i++){
		t1 = pow(2, (l[i]+ m[i]+ n[i]));
		t2 = pow(alpha[i], (2*l[i]+2*m[i]+2*n[i]+3));
		t2 = pow(t2, 0.25);
		t3 = pow((2.0/M_PI), 0.75);
		t4 = fact2(2*l[i]-1)*fact2(2*m[i]-1)*fact2(2*n[i]-1);
		t4 = sqrt(t4);
	}

	//Start at line 111




	return 0;
}