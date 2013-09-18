COMPILER 	= g++ 
CFLAGS 		= -Wno-deprecated -ftemplate-depth-30 -Wno-write-strings 
LIBS 		= -L/sw/lib -lgsl -lgslcblas -w
INCLUDE 	= -I/sw/include 
EXECUTABLE  	= error
OBJECT  	= Error.o

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
		rm -f *~


