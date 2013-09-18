COMPILER 	= g++ 
CFLAGS 		= -Wno-deprecated -ftemplate-depth-30 -Wno-write-strings 
LIBS 		= 
INCLUDE 	= 
EXECUTABLE  	= ../profile
OBJECT  	= main.o WriteNameList.o WriteProfile.o

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
		cd ../Stage1; rm -f *.out *.eps

