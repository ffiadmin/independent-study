#include <iostream>

#include "config.h"
#include "InputParser.h"
#include "Manager.h"

int main(int argc, char** argv) {
	Manager::setFloatPrecision();
	InputParser parser(Config::INPUT_FILE);

	return 0;
}