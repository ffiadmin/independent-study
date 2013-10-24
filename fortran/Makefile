#
#  This is my Makefile
#


#
#  compiler and path settings on radon:
#
CXX = g++
CXXFLAGS = -g -Wall  -Wno-deprecated 

CC = gcc
CFLAGS = -g -Wall

FC =gfortran
FFLAGS = -g -Wall -O -c 

#LDFLAGS = -L. -L

LIBS =  -llapack -lblas


#
#  object lists
#

OBJ = 	main.o  

all	: ptchg

ptchg	: $(OBJ) 
	$(FC) -o $@ $(OBJ)  $(LDFLAGS) $(LIBS)


DXY = 	dxyptq.o    

dxy   : $(DXY) 
	$(FC) -o $@ $(DXY)  $(LDFLAGS) $(LIBS)

PTQ = 	ptq.o    

ptq   : $(PTQ) 
	$(FC) -o $@ $(PTQ)  $(LDFLAGS) $(LIBS)

PLT = 	plot.o    

plot   : $(PLT) 
	$(FC) -o $@ $(PLT)  $(LDFLAGS) $(LIBS)

XCL = 	exclude.o    

exclude  : $(XCL) 
	$(FC) -o $@ $(XCL)  $(LDFLAGS) $(LIBS)


SRR = 	search.o  

search : $(SRR) 
	$(FC) -o $@ $(SRR)  $(LDFLAGS) $(LIBS)

clean	:
	rm *.o *~





