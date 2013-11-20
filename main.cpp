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
#include <Complex>
#include <string>
#include <fstream>
#include <libmints\mints.h>
using namespace psi;

#include "InputParser.h"

#define _USE_MATH_DEFINES
#include <math.h>

typedef std::Complex<double> Complex;

enum {MAXBAS=500, MAXSTP=500, MAXBS3=300, MAXNUC=20};
const char FILENAME[] = "input.txt";

//many global variables discarded. will look to declare them as they are used

#define PTR boost::shared_ptr

int main(){

//=============================================================================
// Initilization (main.f lines 55 - 98)
//=============================================================================
	
	//Psi 4 attempt
					//PTR<BasisSetParser> parser(new Gaussian94BasisSetParser());	
					//parser->load_file(FILENAME);

					////These two need actual information; what that is, I don't know.
					//string type;
					//PTR<Molecule> mol;	

					//PTR<BasisSet> bs = BasisSet::construct(parser, mol, type);
	
	//Homebrew attempt

	InputParser data(FILENAME);

	while(auto temp = data.nextBasis()){
		temp->x /= 0.529177249;
		temp->y /= 0.529177249;
		temp->z /= 0.529177249;
		data.printBasis();
	}	

	//do loop at line 75:
	while(auto temp = data.nextCharge()){		
		temp->nx /= 0.529177249;
		temp->ny /= 0.529177249;
		temp->nz /= 0.529177249;
		data.printNuclei();
	}
	data.reset(); //reset iterators

	data.printBounds();

//=============================================================================
// Calculate normalization constants and integrate (main.f lines 98-183
//=============================================================================
	
	//do loop at line 100:

	double t1, t2, t3, t4;
	while(auto temp = data.nextBasis()){
		t1 = pow(2, (temp->l + temp->m + temp->n));
		t2 = pow(temp->alpha, (2*temp->l + 2*temp->m + 2*temp->n + 3));
		t2 = pow(t2, 0.25);
		t3 = pow((2.0/M_PI), 0.75);
		t4 = fact2(2*temp->l-1)*fact2(2*temp->m-1)*fact2(2*temp->n-1);
		t4 = sqrt(t4);
		temp->norm = t1*t2*t3/t4;
	}

	//do loop at line 111
	//can probably use psi4 integrals -- must...learn...how...
	
	double r1; //What is this?
	int i=0, j=0;

	Complex tint[500][500];
	double squareMatrix[500][500];
	auto it1 = data.basisBegin();
	do{		
		auto it2 = data.basisBegin();
		do{
			double overlap = calculateOverlap(*it1, *it2);
			double kineticEnergy = calculateKineticEnergy(*it1, *it2);
			squareMatrix[i][j] = (*it1)->norm * (*it2)->norm * overlap;
			tint[i][j] = Complex((*it1)->norm * (*it2)->norm * overlap, kineticEnergy);
		}while(++it2 != data.basisEnd());
		i++; j=0;
	}while(++it1 != data.basisEnd());	
	
	//dsyev? line 134

	//skip loop from 136-138

	data.reset();
	Complex vint[MAXBAS][MAXBAS];
	Complex core[MAXBAS][MAXBAS];	
	Complex xmat[MAXBAS][MAXBAS];
	Complex scale;
	
	//loop at 142
	i=0, j=0;
	double alpReal = data.alpStart;	
	while(data.moreAlp())
	{
		double theta = data.thetaStart;
		while (data.moreTheta())
		{			
			double rr1 = alpReal * cos(theta);
			double rr2 = -alpReal * sin(theta);
			scale = Complex(rr1, rr2);

			auto it1 = data.basisBegin();
			do{		
				auto it2 = data.basisBegin();
				do{	
					auto it3 = data.nucleiBegin();
					do{
						Complex venergy = calculatAttraction(*it1, *it2, *it3);
						vint[i][j] = vint[i][j] + Complex((*it3)->nchg,0) * venergy;

					}while(++it3 != data.nucleiEnd());
				
					core[i][j] = scale * scale * tint[i][j] + vint[i][j];
					
					//extra loop from 186-193 moved inside this loop
					xmat[i][j] = Complex(squareMatrix[i][j]/sqrt(svec[j])); //svec must be set in a call to dsyev
					xmat[j][i] = xmat[i][j];
				
				}while(++it2 != data.basisEnd());
				i++; j=0;
			}while(++it1 != data.basisEnd());

//=============================================================================
// Find eigenvalues and eigenvectors (main.f lines 186-253
//=============================================================================
		
			//function calls at 195
			zgemm('N','N',nbasis,nbasis,nbasis,XX1,core,maxbas,xmat,maxbas,XX2,xtmp,maxbas)
			zgemm('N','N',nbasis,nbasis,nbasis,XX1,xmatx,maxbas,xtmp,maxbas,XX2,core,maxbas)
			zgeev('V','V',nbasis,core,maxbas,eigval,VL,maxbas,VR,maxbas,work,maxbs3,rwork,info)
			//the above functions must populate each Basis set with an eigenvalue

			//do loop at 204
			auto it1 = data.basisBegin();
			do{
				Complex temp = ((*it1)->eigenValue *= 27.2114);
				if (temp.real < data.min || temp.real > data.max) continue;
				eta[nscan] = scale;
				val[nscan] = temp;
				ang[nscan] = theta;	

			}while (it1++ != data.basisEnd());

			theta += data.thetaStep;
		}

		double dermin = 1;
		for (int i=0; i<data.thetaCount-1; i++)
		{
			deriv[i] = (val[i+1] - val[i])/(eta[i+1] - eta[i]);
			double derMagnitude = (deriv[i].real * deriv[i].imag);
			if (derMagnitude < dermin)
			{
				dermin = derMagnitude;
				etamin = eta[i];
				valmin = val[i];
				angmin = ang[i];
			}
		}

		data.print(alpReal);
		data.print(angmin);
		data.print(etamin);
		data.print(valmin);
		data.print(ermin);
		data.print('\n');

		alpReal += data.alpStep;
	}

	//END line 253
	//Function library begins

	return 0;
}