/*
Author: Chris Prosser
Course: COMP 220, Computer Programming II
Date: 10 October 2012
Description: This file defines the class Complex
*/

#include "Complex.h"

ostream& operator<<(ostream& out, const Complex& num){
	out << num.real;
	if (num.imag < 0) out << " - ";
	else out << " + ";
	out << abs(num.imag) << 'i';
	return out;
}

istream& operator>>(istream& in, Complex& num){
	in >> num.real;
	in.ignore();
	char sign; in.get(sign);
	in.ignore();
	in >> num.imag;
	if (sign == '-') num.imag *= -1;
	in.ignore();
	return in;
}

Complex operator+(const Complex& first, const Complex& second){
	return Complex(first.real + second.real, first.imag + second.imag);
}

Complex operator-(const Complex& first, const Complex& second){
	return Complex(first.real - second.real, first.imag - second.imag);
}

Complex operator*(const Complex& first, const Complex& second){
	double r = (first.real*second.real - first.imag*second.imag);
	double im = (first.real*second.imag + second.real*first.imag);
	return Complex(r, im);
}

bool operator==(const Complex& first, const Complex& second){
	return (first.real == second.real && first.imag == second.imag);
}

Complex& Complex::operator=(const Complex& rComp){
	imag = rComp.imag;
	real = rComp.real;
	return *this;
}

Complex& Complex::operator+=(const Complex& rComp){
	imag += rComp.imag;
	real += rComp.real;
	return *this;
}