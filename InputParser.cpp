#include "InputParser.h"

InputParser::InputParser(string inputFile, string outputFile) {
	fin.open(inputFile);
	fout.open(outputFile);	

// The first line of the input file will contain two numbers which indicate
// the number of rows which make up the next two sets of data
	fin >> basisCount >> nucleiCount;

	parseBasisSets();
	parseNucleiSets();
	parseBounds();
	reset();
}

InputParser::~InputParser() {
	fin.close();
}

void InputParser::parseBasisSets() {
	Basis* p;
	basisSets.resize(basisCount);

//Read in each of the rows from the text file
	for (int i = 0; i < basisCount; ++i) {
		p = new Basis();

		fin >> p->l
				  >> p->m
				  >> p->n
				  >> p->x
				  >> p->y
				  >> p->z
				  >> p->alpha;

		basisSets[i] = p;
	}

	basisSets.push_back(0); // add a null terminator to ease iteration
	it1 = basisSets.begin(); //set class iterator to beginning of data set

	//Manager::msg("So much boom sauce!", OUTPUT_LEVEL_VERBOSE);
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

	nucleiSet.push_back(0); // add a null terminator to ease iteration
	it2 = nucleiSet.begin(); //set class iterator to beginning of data set

	//Manager::msg("Boom sauce level 2", OUTPUT_LEVEL_VERBOSE);
}

void InputParser::parseBounds(){
	fin >> thetaCount
		>> alpCount
		>> alpStart
		>> alpStep
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
