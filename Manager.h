#ifndef MANAGER_H
#define MANAGER_H

#include <iostream>

#include "config.h"
#include "Display.h"

using std::cout;

class Manager {
	
public : 
	static void msg(char* text, int level);
	static void setFloatPrecision(int precision = Config::FLOAT_PRECISION);
};

#endif