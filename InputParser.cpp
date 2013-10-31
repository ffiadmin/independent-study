#include "InputParser.h"

InputParser::InputParser(string inputFile) {
	this->fin.open(inputFile);

// The first line of the input file will contain two numbers which indicate
// the number of rows which make up the next two sets of data
	this->fin >> this->dataSet1RowCount >> this->dataSet2RowCount;

	this->parseDataSet1();
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

		this->fin >> ds->one
				  >> ds->two
				  >> ds->three
				  >> ds->four
				  >> ds->five
				  >> ds->six
				  >> ds->seven;

		this->dataSet1[i] = ds;
	}

	Manager::msg("So much boom sauce!", OUTPUT_LEVEL_VERBOSE);
}