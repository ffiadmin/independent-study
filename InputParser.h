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

struct DataSet1 {
public : 
	int l;
	int m;
	int n;
	double x;
	double y;
	double z;
	double alpha;
};

struct DataSet2 {
public:
	double nx;
	double ny;
	double nz;
	double nchg;
};

struct Bounds {
public:	
	int nthet;
	int nalp;
	double alpstart;
	double alpstep;
	double thstart;
	double thstep;

	double min;
	double max;
};

class InputParser {
public : 
	InputParser(string inputFile, string outputFile = "output.txt");
	~InputParser();

// returns the next data set in the vector and increments the iterator
	//if there are no more data sets, returns 0 
	//(desired use: while (temp = dataset.next()){} )
	DataSet1* next1(){ return *(it1++); }
	void reset1() { it1 = dataSet1.begin(); }	

	DataSet2* next2(){ return *(it2++); }
	void reset2() { it2 = dataSet2.begin(); }	

	void print1();
	void print2();
	void printBounds();

	//Accessors
	Bounds* getBounds(){return bounds;}
	int getNBasis(){return dataSet1RowCount;}
	

private : 
	vector<DataSet1*> dataSet1;
	vector<DataSet1*>::iterator it1;
	vector<DataSet2*> dataSet2;
	vector<DataSet2*>::iterator it2;
	Bounds* bounds;
	
	int dataSet1RowCount;
	int dataSet2RowCount;
	ifstream fin;
	ofstream fout;
	
	void parseDataSet1();
	void parseDataSet2();
	void parseBounds();

	
};

#endif