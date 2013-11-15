/*
 * main.cpp
 * 
 * This file is a conversion of Dr. Falcetta's FORTRAN main.f file
 * into C++ and all of its glory.  
 *
 * Look! Objects and stuff!
 *
 */

#include <iostream>
#include <complex>
#include <string>
#include <fstream>
#include <libmints\mints.h>
using namespace psi;

#include "InputParser.h"

#define _USE_MATH_DEFINES
#include <math.h>

typedef std::complex<double> Complex;

enum {MAXBAS=500, MAXSTP=500, MAXBS3=300, MAXNUC=20};
const char FILENAME[] = "";

//many global variables discarded. will look to declare them as they are used

#define PTR boost::shared_ptr

int main(){

//=============================================================================
// Initilization (main.f lines 55 - 98)
//=============================================================================
	
	PTR<BasisSetParser> parser(new Gaussian94BasisSetParser());	
	parser->load_file(FILENAME);

	//These two need actual information; what that is, I don't know.
	string type;
	PTR<Molecule> mol;	

	PTR<BasisSet> bs = BasisSet::construct(parser, mol, type);
	
	
	
	
	InputParser data("input.txt");

	while(auto temp = data.next1()){
		temp->x /= 0.529177249;
		temp->y /= 0.529177249;
		temp->z /= 0.529177249;
		data.print1();
	}	
	data.reset1(); // reset iterator to beginning

	//do loop at line 75:
	while(auto temp = data.next2()){		
		temp->nx /= 0.529177249;
		temp->ny /= 0.529177249;
		temp->nz /= 0.529177249;
		data.print2();
	}
	data.reset2();

	data.printBounds();

//=============================================================================
// Calculate normalization constants and integrate (main.f lines 98-183
//=============================================================================
	
	//do loop at line 100:

	double t1, t2, t3, t4;
	while(auto temp = data.next1()){
		t1 = pow(2, (temp->l + temp->m + temp->n));
		t2 = pow(temp->alpha, (2*temp->l + 2*temp->m + 2*temp->n + 3));
		t2 = pow(t2, 0.25);
		t3 = pow((2.0/M_PI), 0.75);
		//t4 = fact2(2*temp->l-1)*fact2(2*temp->m-1)*fact2(2*temp->n-1);
		t4 = sqrt(t4);
	}

	//do loop at line 111
	//can probably use psi4 integrals -- must...learn...how...
	data.reset1();
	while(auto temp = data.next1()){}
	
	//dsyev? line 134

	//skip loop from 136-138

	Complex vint[MAXBAS][MAXBAS];
	Complex core[MAXBAS][MAXBAS];
	
	//loop at 142
	for (int i=0; i < data.getBounds()->nalp; i++){
		for (int j=0; j<data.getBounds()->nthet; j++){
			
			double r = data.getBounds()->alpstart;
			double theta = data.getBounds()->thstart;
			Complex scale(r*cos(theta), r*sin(theta));
			int nbasis = data.getNBasis();

			for (int k=0; k < nbasis; k++){
				for (int l=0; l < nbasis; l++){
					
					//inner loop at line 160
					while(auto temp = data.next2()){
						//psi4 attraction function

						Complex c6(temp->nchg, 0);
						vint[k][l] = c6 * venergy; // <- output parameter in attract function
					}

					core[k][l] = scale*scale*tint[k][l]+vint[k][l];
				}
			}
		}
	}

//=============================================================================
// Find eigenvalues and eigenvectors (main.f lines 186-253
//=============================================================================

	return 0;
}