COMPILER 	= g++ 
CFLAGS 		= -Wno-deprecated -ftemplate-depth-30 -Wno-write-strings 
LIBS 		= -L/sw/lib -lgsl  -lgslcblas -w
INCLUDE 	= -I/sw/include
EXECUTABLE  	= ../flux
OBJECT  	= main.o Flux.o Stage1.o Stage2.o EQDSK_Read.o CalcQGP.o CalcrP.o Calcst.o Profile.o\
	Metric.o FTCos.o FTSin.o InterpolatePsi.o Interpolate.o Interpolate1.o Extrapolate1.o\
	Extrapolate2.o Extrapolate3.o Misc.o Boundary.o ToroidalP.o ToroidalQ.o

# Linking everything:

all:		$(EXECUTABLE)

$(EXECUTABLE):	$(OBJECT)
		$(COMPILER) -o $(EXECUTABLE) $(OBJECT) $(LIBS)

%.o:		%.cpp 
		$(COMPILER) $(CFLAGS) -o $*.o -c $(INCLUDE) $*.cpp

# Cleanup

clean: 
		rm -f $(OBJECT)

clear:
		cd ../Stage2; rm -f *.out *.eps

