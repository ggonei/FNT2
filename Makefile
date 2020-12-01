#	Makefile for FNT2
#	Code is C++98 and later compliant but needs a 64-bit system
#	George O'Neill, 21/11/2020

#	Use compiler
CC = g++ -m64

#	Point to source files
SOURCES = $(shell ls FNT2/src/*.cpp)

#	Use source names to get binary object files
OBJECTS = $(SOURCES:.cpp=.o)

#	Include ROOT libraries without warnings
ROOTINC = -isystem /home/gon/root_install/include

#	Get ROOT flags
ROOTCXXFLAGS := $(shell root-config --cflags)
#	LDFLAGS := $(shell root-config --libs)
#LDFLAGS = -L/home/gon/root_install/lib -pthread -rdynamic 
#LDFIX = -Wl,--no-as-needed -lCore -lRIO -lTree -lHist -lSpectrum
LDFIX = -pthread -rdynamic -L/home/gon/root_install/lib/ -Wl,-rpath -Wl,/home/gon/root_install/lib/ -Wl,--no-as-needed -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -ltbb -lm -ldl -lSpectrum

#	All useful warnings are turned on as errors (-Werror), along with debug symbols and optimisation
CXXFLAGS = -c -g -O3 $(ROOTINC) $(ROOTCXXFLAGS) -Wall -pedantic -pedantic-errors -Wextra -Werror -Wcast-align -Wcast-qual -Wchar-subscripts -Wcomment -Wconversion -Wdisabled-optimization -Wfloat-equal -Wformat -Wformat=2 -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wimport -Winit-self -Winline -Winvalid-pch -Wlong-long -Wmissing-braces -Wmissing-field-initializers -Wmissing-format-attribute -Wmissing-include-dirs -Wmissing-noreturn -Wpacked -Wparentheses -Wpointer-arith -Wredundant-decls -Wreturn-type -Wsequence-point -Wshadow -Wsign-compare -Wstack-protector -Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch -Wswitch-default -Wswitch-enum -Wtrigraphs -Wuninitialized -Wunknown-pragmas -Wunreachable-code -Wunused -Wunused-function -Wunused-label -Wunused-parameter -Wunused-value -Wunused-variable -Wvariadic-macros -Wvolatile-register-var -Wwrite-strings -I./include/

#	We have to turn off some errors
CXXFLAGS += -Wno-long-long

#	Make the executable
fnt2: $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(LDFIX) -o $@.exe

.cpp.o:
	$(CC) $(CXXFLAGS) $< -o $@

#	Clean up built program
clean:
	rm -f src/*_cpp.* src/*.o fnt2.exe
