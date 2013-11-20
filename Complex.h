/*
Author: Chris Prosser
Course: COMP 220, Computer Programming II
Date: 10 October 2012
Description: This file declares the class Complex
*/

#ifndef __COMPLEX_
#define __COMPLEX_

#include <iostream>
using std::ostream;
using std::istream;

class Complex{
public:
	Complex(){real = 0; imag = 0;}
	Complex(double r, double i){real = r; imag = i;}
	Complex(double real_part){real = real_part; imag = 0;}

	friend ostream& operator<<(ostream& out, const Complex& num);
	friend istream& operator>>(istream& in, Complex& num);

	friend Complex operator+(const Complex& first, const Complex& second);
	friend Complex operator-(const Complex& first, const Complex& second);
	friend Complex operator*(const Complex& first, const Complex&  second);	
	
	friend bool operator==(const Complex& first, const Complex& second);	

	double real, imag;
};

#endif