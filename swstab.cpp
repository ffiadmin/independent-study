/**********************************************************************
 * This program uses the one-electron integrals from:
 * cints.c  C implementation of simple math functions in pyutil.
 *
 * The equations herein are based upon
 * 'Gaussian Expansion Methods for Molecular Orbitals.' H. Taketa,
 * S. Huzinaga, and K. O-ohata. H. Phys. Soc. Japan, 21, 2313, 1966.
 * [THO paper].
 *
 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 

 **********************************************************************/

#include "_swFunc.h"

using namespace std;

#if defined(_WIN32)
double lgamma(double x);
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define ITMAX 100
#define EPS 3.0e-12
#define FPMIN 1.0e-30
#define SMALL 0.0000000000001

extern "C" void FORTRANNAME(dsyev)(char *job, char *uplo, int *n, double *a, int *lda,
                                   double *w, double *aux, int *naux, int *info);
extern "C" void FORTRANNAME(dgemm)(char *transa, char *transb, int *m, int *n, int *k,
                                   double *alpha, double *a, int *lda, double *b, int *ldb, double *beta,
                                   double *c, int *ldc);

extern "C" void FORTRANNAME(qromb)(int *l, double *alp,  double *answer);                     

int main(){
	
	ifstream fin("model_input");
   
    cout.setf(ios::fixed);
    cout.setf(ios::showpoint);
    cout.precision(16);
		
	double minalp, maxalp, step;
	fin >> minalp >> maxalp >> step;

    cout << "minalp = "  << minalp << endl;
    cout << "maxalp = "  << maxalp << endl;
    cout << "step   = "  << step   << endl;

	int nbasis,nnuc;
    fin >> nbasis >> nnuc;

    cout << "nbasis =  " << nbasis << endl;
    cout << "nnuc =  " << nnuc << endl;
	
	int		l[nbasis], m[nbasis], n[nbasis],
			naux = 3*nbasis-1,
			info;

	double	alpha[nbasis],
			x[nbasis], y[nbasis], z[nbasis],
			norm[nbasis],
			nchg[nnuc],
			nx[nnuc], ny[nnuc], nz[nnuc],
			Tint[nbasis][nbasis],
			Smat[nbasis][nbasis],
			Vint[nbasis][nbasis],
			core [nbasis][nbasis],
			Mint[nbasis][nbasis],
			Cmat[nbasis][nbasis],
			Xmat[nbasis][nbasis],
			Fmat[nbasis][nbasis],
			Fmat1[nbasis][nbasis],
			XmatTrans[nbasis][nbasis],
			alp_standard[nbasis],
			fvec[nbasis],
			svec[nbasis],
			aux[naux],
			tempalpha = 1.0,
			tempbeta = 0.0,
			t1,t2,t3,t4,
			min, max;

  	for(int i=0; i<nbasis; i++) fin >> l[i] >> m[i] >> n[i] >> alp_standard[i] >> x[i] >> y[i] >> z[i];        

  	for(int i=0; i<nnuc; i++){
            fin >> nx[i] >> ny[i] >> nz[i] >> nchg[i];
            nx[i]=nx[i]/0.529177249;	//What is this number?? Could be a constant
            ny[i]=ny[i]/0.529177249;
            nz[i]=nz[i]/0.529177249;
        }
	
	fin >> min >> max;
	cout << "min =  " << min  << endl;
	cout << "max =  " << max  << endl;
	fin.close();

	for (double scale=minalp; scale <=maxalp; scale=scale+step){ 
		
		for (int i=0; i<nbasis; i++) alpha[i]=scale*alp_standard[i];	

		//calculate normalization constants
		for(int i=0; i<nbasis; i++){
			t1 = pow(2.0,(double)(l[i]+m[i]+n[i]));
			t2 = pow(alpha[i],(double)(2*(l[i]+m[i]+n[i])+3));
			t2 = pow(t2,.25);
			t3 = pow(2/M_PI,0.75);
			t4 = sqrt((double)(fact2(2*l[i]-1)*fact2(2*m[i]-1)*fact2(2*n[i]-1)));
			norm[i]=t1*t2*t3/t4;
		}
		
		//calculate the T, S, and V ints and form core

		for (int row=0; row<nbasis; row++){
			for (int col=0; col<nbasis; col++){
				Tint[row][col] = 
					norm[row] *
					norm[col] *
					kinetic(alpha[row],l[row],m[row],n[row],x[row],y[row],z[row],
							alpha[col],l[col],m[col],n[col],x[col],y[col],z[col]);
	
				Smat[row][col] = 
					norm[row] * 
					norm[col] *
					overlap(alpha[row],l[row],m[row],n[row],x[row],y[row],z[row],
							alpha[col],l[col],m[col],n[col],x[col],y[col],z[col]);

 			// This is a call to a generic model potential numerical integration routine.
				Mint[row][col]= 
					4.0 * M_PI * 
					norm[row] * 
					norm[col] * 
					Model_pot(l[row],l[col],m[row],m[col],n[row],n[col],
							alpha[row],alpha[col]);

				Vint[row][col]=0.0;

		     	core[row][col] = Tint[row][col] + Mint[row][col];
			}
		}

//              find eigenvals and vecs for the S matrix, form the X and Xtranspose and transfor core to orthonormal
//              basis set before diagonalizing.

        FORTRANNAME(dsyev)("V", "U", &nbasis, *Smat, &nbasis,  svec, aux, &naux, &info);	

        for (int row=0; row<nbasis; row++){
       		for (int col=0; col<nbasis; col++){               	              
           		 Xmat[row][col] = Smat[col][row]/sqrt(svec[col]);
				 XmatTrans[col][row]=Xmat[row][col];
               	}
        }

//		transform core to orthonormal basis and store in Fmat1

		FORTRANNAME(dgemm)("N","N",&nbasis,&nbasis,&nbasis,&tempalpha,*Xmat,&nbasis,*core,&nbasis,&tempbeta,*Fmat,&nbasis);
		FORTRANNAME(dgemm)("N","N",&nbasis,&nbasis,&nbasis,&tempalpha,*Fmat,&nbasis,*XmatTrans,&nbasis,&tempbeta,*Fmat1,&nbasis);
		FORTRANNAME(dsyev)("V", "U", &nbasis, *Fmat1, &nbasis,  fvec, aux, &naux, &info);
		FORTRANNAME(dgemm)("N","N",&nbasis,&nbasis,&nbasis,&tempalpha,*XmatTrans,&nbasis,*Fmat1,&nbasis,&tempbeta,*Cmat,&nbasis);
		
		for (int i=0; i<nbasis; i++){
			fvec[i]=fvec[i];
			if ((fvec[i]>min)&&(fvec[i]<max)){
				cout << scale << "    " << fvec[i]  << endl;
			}
		}

	}
}       
     
