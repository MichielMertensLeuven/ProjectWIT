#MAKEFILE FOR WIT PROJECT

.PHONY: clean all info

CXX = g++

TARGETS := print main
SOURCES := pprint.cpp main.cpp
OBJS    := $(TARGETS:=.o)

OFLAGS := -O2 -O3 -ffast-math
PARALLELFLAGS:= -D_GLIBCXX_PARALLEL -fopenmp -pthread -DUSEOMP
CXXFLAGS := -std=c++14
LDFLAGS  := -lm -ldl
LIBS :=  -lstdc++ -I./Eigen
LIPOPT := -L/localhost/packages/math_eng_soft/CoinIpopt/build/lib -I/localhost/packages/math_eng_soft/CoinIpopt/build/include/coin -lipopt
OPENMPLIB :=-L/usr/lib -fopenmp -pthread

EXAMPLE_DEPS = Makefile

all: main

clean:
	rm -f $(OBJS) $(TARGETS)

info:
	@echo Compiler:	 CXX	  = $(CXX)
	@echo Compile command: COMPILE.cc  = $(COMPILE.cc)
	@echo Link command:    LINK.cc     = $(LINK.cc)

print.o: print.cpp $(EXAMPLE_DEPS)
	@$(CXX) -c $(CXXFLAGS) $(OFLAGS) $(LIBS) -o print.o print.cpp
main.o: main.cpp $(EXAMPLE_DEPS)
	@$(CXX) -c $(CXXFLAGS) $(OFLAGS) $(LIBS) $(LIPOPT) -o main.o main.cpp
main: main.o print.o
	@$(CXX) $(LDFLAGS) $(LIPOPT) -o main main.o print.o $(LIBS)
