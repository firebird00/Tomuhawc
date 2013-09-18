COMPILER 	= g++ -O3
CFLAGS 		= -Wno-deprecated -ftemplate-depth-30 -Wno-write-strings 
LIBS 		= -L/sw/lib -lgsl  -lgslcblas -w
INCLUDE 	= -I/sw/include
EXECUTABLE  	= ../ggj
OBJECT  	= main.o GGJ.o 

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
		cd ../Stage4; rm -f *.out *.eps


