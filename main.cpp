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

//many global variables discarded. will look to declare them as they are used

int main(){

//=============================================================================
// Initilization (main.f lines 55 - 98)
//=============================================================================
	IntegralFactory integrals(BasisSet(;
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
// Calculate normalization constants (main.f lines 98-183
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
	data.reset1();
	
	vector<SphericalTransform> v;	
	while(auto temp = data.next1())	v.push_back(SphericalTransform(temp->l));
	}
	OverlapInt(SphericalTransform())
			

	
	
	

	return 0;
}