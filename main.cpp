#include "InputParser.h"
#include <iostream>
using std::cin;
using std::cout;

void runFortran(InputParser& data, int thetaCount, double thetaStart);

enum {REAL = 1, COMPLEX = 2};

extern "C"
{
	void _stdcall setNBasis(int);
	void _stdcall setNNuc(int);
	void _stdcall setBasisCoords(int, int, int, int, double, double, double, double);
	void _stdcall setNucleusCoords(int, double, double, double, double);	
	void _stdcall setFinalParams(int nThetVal, int nAlpVal, double alpStrtVal, double alpStpVal, double thStartVal, double thStepVal, int minVal, int maxVal);
	void _stdcall go();
}

int main()
{
	cout << "Input file: ";
	string inputFile;
	cin >> inputFile;	
	
	InputParser data(inputFile);

	setNBasis(data.basisCount);
	setNNuc(data.nucleiCount);	

	cout << "1 Real\n2 Complex\n\n";	
	int ID; cin >> ID;
	
	switch(ID){
	case REAL: runFortran(data, 1, 0); break;
	case COMPLEX: runFortran(data, data.thetaCount, data.thetaStart); break;
	}

	return 0;
}

void runFortran(InputParser& data, int thetaCount, double thetaStart)
{
	int index = 1;
	for (auto it = data.basisBegin(); it < data.basisEnd(); ++it)
	{
		setBasisCoords(index,	(*it)->l, (*it)->m, (*it)->n, 
								(*it)->x, (*it)->y, (*it)->z, 
								(*it)->alpha);
		index++;
	}

	index = 1;
	for (auto it = data.nucleiBegin(); it < data.nucleiEnd(); ++it)
	{
		setNucleusCoords(index, (*it)->nx, (*it)->ny, (*it)->nz, (*it)->nchg);  
		index++;
	}

	setFinalParams(	thetaCount, data.alpCount, data.alpStart, data.alpStep, 
					thetaStart, data.thetaStep, data.min, data.max);

	go();
}