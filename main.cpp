#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <ctype.h>

using std::cin;
using std::cout;
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;

//Exception class
class FileNotFound{};

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
	InputParser(string inputFile);
	~InputParser() {fin.close();}

	//Iterate
	Basis* nextBasis()		{return *(it1++);}
	Nuclei* nextCharge()	{return *(it2++);}

	//Reset iterators and incremetors
	void reset();	

	//Print
	void printBasis();
	void printNuclei();
	void printBounds();
	
	template<class T>
	void print(T output) {fout << output << ' ';}

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

InputParser::InputParser(string inputFile) {
	fin.open(inputFile.c_str());
	if (fin.fail()) throw FileNotFound();

// The first line of the input file will contain two numbers which indicate
// the number of rows which make up the next two sets of data
	fin >> basisCount >> nucleiCount;
	parseBasisSets();
	parseNucleiSets();
	parseBounds();
	reset();
}

void InputParser::parseBasisSets() {
	Basis* p;
	basisSets.resize(basisCount);

//Read in each of the rows from the text file
	for (int i = 0; i < basisCount; ++i) {
		p = new Basis();

		fin >>	p->l
				  >> p->m
				  >> p->n
				  >> p->x
				  >> p->y
				  >> p->z
				  >> p->alpha;

		basisSets[i] = p;
	}
}

void InputParser::parseNucleiSets(){
	Nuclei* p;
	nucleiSet.resize(nucleiCount);

	for (int i=0; i < nucleiCount; i++){
		p = new Nuclei();

		fin >> p->nx
				  >> p->ny
				  >> p->nz
				  >> p->nchg;

		nucleiSet[i] = p;
	}
}

void InputParser::parseBounds(){
	fin >> alpCount
		>> alpStart
		>> alpStep
		>> thetaCount
		>> thetaStart
		>> thetaStep
		>> min
		>> max;
}

void InputParser::printNuclei(){
	fout << (*it2)->nx << ' '
		       << (*it2)->ny << ' '
			   << (*it2)->nz << ' '
			   << (*it2)->nchg << '\n';			
}

void InputParser::printBasis(){
	fout << (*it1)->l << ' '
		       << (*it1)->m << ' '
			   << (*it1)->n << ' '
			   << (*it1)->x << ' '
			   << (*it1)->y << ' '
			   << (*it1)->z << ' '
			   << (*it1)->alpha << '\n';
}

void InputParser::printBounds(){
	fout << "\nThe initial values for r and theta are " << alpStart << " and " << thetaStart;
	fout << "\nThe program will take " << thetaCount << " steps.";
	fout << "\nEach step in theta will be " << thetaStep << " radians.";
	fout << "\nThe program will take " << alpCount << " steps.";
	fout << "\nEach step in r will be " << alpStep << " au";
	fout << "\nThe range for the real part of the energy is from " << min << " to " << max;
}

void InputParser::reset(){
	it1 = basisSets.begin(); 
	it2 = nucleiSet.begin(); 
	currentAlp = currentTheta = currentBasis = 0;
}

void runFortran(InputParser& data, int thetaCount, double thetaStart);

#ifndef FORTRANNAME
	#ifdef IBM
		#define FORTRANNAME(name) name
	#else
		#define FORTRANNAME(name) name##_
	#endif
#endif

extern "C"
{
	void FORTRANNAME(setnbasis)(int*);
	void FORTRANNAME(setnnuc)(int*);
	void FORTRANNAME(setbasiscoords)(int*, int*, int*, int*, double*, double*, double*, double*);
	void FORTRANNAME(setnucleuscoords)(int*, double*, double*, double*, double*);	
	void FORTRANNAME(setfinalparams)(int* nThetVal, int* nAlpVal, double* alpStrtVal, double* alpStpVal, double* thStartVal, double* thStepVal, double* minVal, double* maxVal);
	void FORTRANNAME(go)();
}

int main(int argc, char *argv[])
{	
	if (argc < 3)
	{
		cout << "Invalid number of arguments\n"
				<< "<executable> <input filename> <\"real\" OR \"complex\">\n"; 
		return -1;
	}	

	InputParser data(argv[1]);	

	FORTRANNAME(setnbasis)(&data.basisCount);
	FORTRANNAME(setnnuc)(&data.nucleiCount);	

	string input = argv[2];
	for (int i=0; i<input.size(); i++) input[i] = tolower(input[i]);
	
	if (input == "real") runFortran(data, 1, 0);
	else if (input == "complex") runFortran(data, data.thetaCount, data.thetaStart);
	else
	{
		cout << "Invalid second parameter.  Choose either \"real\" or \"complex\".\n";
		return -1;
	}

	return 0;
}

void runFortran(InputParser& data, int thetaCount, double thetaStart)
{
	int index = 1;
	
	for (vector<Basis*>::iterator it = data.basisBegin(); it < data.basisEnd(); ++it)
	{
		FORTRANNAME(setbasiscoords)(&index,	&(*it)->l, &(*it)->m, &(*it)->n, 
								&(*it)->x, &(*it)->y, &(*it)->z, 
								&(*it)->alpha);
		index++;
	}


	index = 1;
	for (vector<Nuclei*>::iterator it = data.nucleiBegin(); it < data.nucleiEnd(); ++it)
	{
		FORTRANNAME(setnucleuscoords)(&index, &(*it)->nx, &(*it)->ny, &(*it)->nz, &(*it)->nchg);  
		index++;
	}

	FORTRANNAME(setfinalparams)(&thetaCount, &data.alpCount, &data.alpStart, &data.alpStep, 
					&thetaStart, &data.thetaStep, &data.min, &data.max);

	FORTRANNAME(go)();
}
