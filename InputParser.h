#ifndef INPUT_PARSER_H
#define INPUT_PARSER_H

#include <fstream>
#include <string>
#include <vector>

#include "Manager.h"

using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;

//needs some overloads -- check old code from COMP 220
struct complex
{
	complex():real(0), imag(0){}
	complex(double d) : real(d), imag(0){}
	complex(double r, double i) : real(r), imag(i){}
	double real, imag;
}


struct Basis {
public : 
	int l;
	int m;
	int n;
	double x;
	double y;
	double z;
	double alpha;
	double norm;
	complex eigenValue;
};

struct Nuclei {
public:
	double nx;
	double ny;
	double nz;
	double nchg;
};


class InputParser {
public : 
	InputParser(string inputFile, string outputFile = "output.txt");
	~InputParser();

	//Iterate
	Basis* nextBasis()		{return *(it1++);}
	Nuclei* nextCharge()	{return *(it2++);}
	bool moreAlp()			{return currentAlp++ < alpCount;}
	bool moreTheta()		{return currentTheta++ < thetaCount;}
	bool moreBasis()		{return currentBasis++ < basisCount;}
	bool moreNuclei()		{return currentNuclei++ < nucleiCount;}

	//Reset iterators and incremetors
	void reset();	

	//Print
	void printBasis();
	void printNuclei();
	void printBounds();

	//Public local variables
	int basisCount;
	int nucleiCount;
	int thetaCount;
	int alpCount;
	double alpStart;
	double alpStep;
	double thetaStart;
	double thetaStep;
	double min;
	double max;

	//Access start and end iterators
	vector<Basis*>::iterator basisBegin() {return basisSets.begin();}
	vector<Basis*>::iterator basisEnd() {return basisSets.end();}
	vector<Nuclei*>::iterator nucleiBegin() {return nucleiSet.begin();}
	vector<Nuclei*>::iterator nucleiEnd() {return nucleiSet.begin();}
	

private : 
	//File streams
	ifstream fin;
	ofstream fout;
	
	//Parsing functions
	void parseBasisSets();
	void parseNucleiSets();
	void parseBounds();
	
	//Data storage
	vector<Basis*> basisSets;	
	vector<Nuclei*> nucleiSet;	

	//Iteration
	int currentAlp, currentTheta, currentBasis, currentNuclei; 	//used to iterate over alp and theta
	vector<Basis*>::iterator it1;
	vector<Nuclei*>::iterator it2;
};

#endif