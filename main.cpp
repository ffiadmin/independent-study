#include "InputParser.h"
#include <iostream>
using std::cin;
using std::cout;

const char* INPUT_FILE = "input.txt";

void runRealCalc(InputParser& data);
void runComplexCalc(InputParser& data);

enum {REAL = 1, COMPLEX = 2};

int main()
{
	InputParser data(INPUT_FILE);

	cout << "1 Real\n2 Complex\n\n";	
	int ID; cin >> ID;
	
	switch(ID){
	case REAL: runRealCalc(data); break;
	case COMPLEX: runComplexCalc(data); break;
	}

	return 0;
}

void runRealCalc(InputParser& data)
{	
	cout << "Real\n";
}

void runComplexCalc(InputParser& data)
{
	cout << "Complex\n";
}