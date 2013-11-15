#include "InputParser.h"

InputParser::InputParser(string inputFile, string outputFile) {
	this->fin.open(inputFile);
	this->fout.open(outputFile);

// The first line of the input file will contain two numbers which indicate
// the number of rows which make up the next two sets of data
	this->fin >> this->dataSet1RowCount >> this->dataSet2RowCount;

	this->parseDataSet1();
	this->parseDataSet2();
	this->parseBounds();
}

InputParser::~InputParser() {
	this->fin.close();
}

void InputParser::parseDataSet1() {
	DataSet1* ds;
	this->dataSet1.resize(this->dataSet1RowCount);

//Read in each of the rows from the text file
	for (int i = 0; i < this->dataSet1RowCount; ++i) {
		ds = new DataSet1();

		this->fin >> ds->l
				  >> ds->m
				  >> ds->n
				  >> ds->x
				  >> ds->y
				  >> ds->z
				  >> ds->alpha;

		this->dataSet1[i] = ds;
	}

	dataSet1.push_back(0); // add a null terminator to ease iteration
	it1 = dataSet1.begin(); //set class iterator to beginning of data set

	Manager::msg("So much boom sauce!", OUTPUT_LEVEL_VERBOSE);
}

void InputParser::parseDataSet2(){
	DataSet2* ds;
	dataSet2.resize(dataSet2RowCount);

	for (int i=0; i < this->dataSet1RowCount; i++){
		ds = new DataSet2();

		this->fin >> ds->nx
				  >> ds->ny
				  >> ds->nz
				  >> ds->nchg;

		dataSet2[i] = ds;
	}

	dataSet2.push_back(0); // add a null terminator to ease iteration
	it2 = dataSet2.begin(); //set class iterator to beginning of data set

	Manager::msg("Boom sauce level 2", OUTPUT_LEVEL_VERBOSE);
}

void InputParser::parseBounds(){
	bounds = new Bounds();

	this->fin >> bounds->nthet
		      >> bounds->nalp
			  >> bounds->alpstart
			  >> bounds->alpstep
			  >> bounds->thstart
			  >> bounds->thstep
			  >> bounds->min
			  >> bounds->max;
}

void InputParser::print2(){
	this->fout << (*it2)->nx << ' '
		       << (*it2)->ny << ' '
			   << (*it2)->nz << ' '
			   << (*it2)->nchg << '\n';			
}

void InputParser::print1(){
	this->fout << (*it1)->l << ' '
		       << (*it1)->m << ' '
			   << (*it1)->n << ' '
			   << (*it1)->x << ' '
			   << (*it1)->y << ' '
			   << (*it1)->z << ' '
			   << (*it1)->alpha << '\n';
}

void InputParser::printBounds(){
	fout << "\nThe initial values for r and theta are " << bounds->alpstart << " and " << bounds->thstart;
	fout << "\nThe program will take " << bounds->nthet << " steps.";
	fout << "\nEach step in theta will be " << bounds->thstep << " radians.";
	fout << "\nThe program will take " << bounds->nalp << " steps.";
	fout << "\nEach step in r will be " << bounds->alpstep << " au";
	fout << "\nThe range for the real part of the energy is from " << bounds->min << " to " << bounds->max;
}
