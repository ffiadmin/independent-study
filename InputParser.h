#ifndef INPUT_PARSER_H
#define INPUT_PARSER_H

#include <fstream>
#include <string>
#include <vector>

#include "Manager.h"

using std::ifstream;
using std::string;
using std::vector;

struct DataSet1 {
public : 
	int one;
	int two;
	int three;
	double four;
	double five;
	double six;
	double seven;
};

class InputParser {
public : 
	InputParser(string inputFile);
	~InputParser();

private : 
	vector<DataSet1*> dataSet1;
	int dataSet1RowCount;
	int dataSet2RowCount;
	ifstream fin;

	void parseDataSet1();
};

#endif