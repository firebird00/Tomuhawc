COMPILER 	= g++ -O3
CFLAGS 		= -Wno-deprecated -ftemplate-depth-30 -Wno-write-strings 
LIBS 		= -L/sw/lib -lgsl  -lgslcblas -w
INCLUDE 	= -I/sw/include
EXECUTABLE 	= ../thawc
OBJECT 		= main.o Tomuhawc.o Solv.o Boundary.o CalcFee.o CalcFoo.o CalcFeo.o\
		CalcFoe.o CalcEigEe.o CalcEigFee.o Equilibrium.o GGJCalc.o Resonant.o\
		ZeroFind.o Fixup.o Launch.o Launch1.o Segment_Fixup.o Segment.o Segment1.o Rhs.o\
		Couple.o Jump.o Jump1.o L1P1T1.o L1kP1k.o Log.o CalcEe.o CalcEo.o CalcEigEo.o\
		CalcEigFoo.o CalcPsiMax.o CalcEigFeea.o CalcEigEea.o CalcEigFooa.o CalcEigEoa.o\
		Vacuum.o

# Linking everything:

all:		$(EXECUTABLE) 

$(EXECUTABLE):  $(OBJECT)
	       	$(COMPILER) -o $(EXECUTABLE) $(OBJECT) $(LIBS)

%.o:		%.cpp
		$(COMPILER) $(CFLAGS) -o $*.o -c $(INCLUDE) $*.cpp

# Cleanup

clean: 
	  	rm -f $(OBJECT)

clear:
	  	cd ../Stage3; rm -f *.out *.eps


