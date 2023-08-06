# make file to provide an easier way to
# build matters than with CMake

EXE=WarfareIBM

CPP=ibm_warfare.cpp main.cpp individual.cpp patch.cpp
CXX=g++
CXXFLAGS=-Wall -O3 -std=c++17 -ggdb

$(EXE) : $(CPP)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(CPP)


clean :
	rm -rf $(EXE)
