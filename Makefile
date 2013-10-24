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

FC = gfortran
FFLAGS = -g -Wall -O -c 

LDFLAGS = -L. 
LIBS = -llapack -lblas 

DDIR  = $(HOME)/grid/drudemodel

#
#  object lists
#

OBJ = 	test1.o      \

all	: ptchg

ptchg	: $(OBJ) 
	$(CXX) -o $@ $(OBJ)  $(LDFLAGS) $(LIBS)

SWC = swcap.o

swcap   : $(SWC)
	$(CXX) -o $@ $(SWC) smyeig.o wqromb.o sqromb.o  $(LDFLAGS) $(LIBS)


SWV = swave.o

swave   : $(SWV)
	$(CXX) -o $@ $(SWV) smyeig.o sqromb.o  $(LDFLAGS) $(LIBS)

SWT = swstab.o

swstab   : $(SWT)
	$(CXX) -o $@ $(SWT) sqromb.o  $(LDFLAGS) $(LIBS)

MOD = model.o

model   : $(MOD)
	$(CXX) -o $@ $(MOD) myeig.o qromb.o  $(LDFLAGS) $(LIBS)

CAP = cap.o

cap     : $(CAP)
	$(CXX) -o $@ $(CAP) myeig.o  $(LDFLAGS) $(LIBS)

MYEIG = myeig.o

myeig   : $(MYEIG)
#	$(FC)  $(MYEIG) 


clean	:
	rm *.o *~





#
#  dependencies on header files
#

dvr3d.o	: fulldiag.h larnoldi.h maxdim.h setupdvr.h tsin.h 

fulldiag.o	: constants.h maxdim.h fortran.h

ho_dvr.o	: fortran.h

larnoldi.o	: constants.h fortran.h mtx.h tsin.h writewfcuts.h

mtx.o	: fortran.h maxdim.h

setupdvr.o : cm_dvr.h ho_dvr.h potentials.h sine_dvr.h tsin.h 

