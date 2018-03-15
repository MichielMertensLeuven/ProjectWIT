#MAKEFILE FOR WIT PROJECT

.PHONY: clean all info

CXX = g++

TARGETS := plate heat solver print opti main
SOURCES := plate.cpp heat.cpp solver.cpp print.cpp opti.cpp main.cpp
OBJS    := $(TARGETS:=.o)

OFLAGS := -O2 -O3 -ffast-math
PARALLELFLAGS:= -D_GLIBCXX_PARALLEL -fopenmp -pthread -DUSEOMP
CXXFLAGS := -Wall -std=c++14
LDFLAGS  :=
LIBS :=  -lstdc++ -I../Eigen
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

plate.o: plate.cpp $(EXAMPLE_DEPS)
	@$(CXX) -c $(CXXFLAGS) $(OFATSFLAGS) $(LIBS) -o plate.o plate.cpp
heat.o: heat.cpp $(EXAMPLE_DEPS)
	@$(CXX) -c $(CXXFLAGS) $(OFATSFLAGS) $(LIBS) -o heat.o heat.cpp
solver.o: solver.cpp $(EXAMPLE_DEPS)
	@$(CXX) -c $(CXXFLAGS) $(OFLAGS) $(LIBS) $(LIPOPT) -o solver.o solver.cpp
print.o: print.cpp $(EXAMPLE_DEPS)
	@$(CXX) -c $(CXXFLAGS) $(OFLAGS) $(LIBS) -o print.o print.cpp
opti.o: opti.cpp $(EXAMPLE_DEPS)
	@$(CXX) -c $(CXXFLAGS) $(OFLAGS) $(LIBS) $(LIPOPT) -o opti.o opti.cpp
main.o: main.cpp $(EXAMPLE_DEPS)
	@$(CXX) -c $(CXXFLAGS) $(OFLAGS) $(LIBS) -o main.o main.cpp
main: main.o plate.o heat.o solver.o print.o opti.o
	@$(CXX) $(LDFLAGS) -o main main.o plate.o heat.o solver.o print.o opti.o $(LIBS)
