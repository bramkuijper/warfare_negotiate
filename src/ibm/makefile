# make file to provide an easier way to
# build matters than with CMake

EXE=MutualismIBM

CPP=ibm_mutualism.cpp ibm_mutualism_main.cpp individual.cpp patch.cpp
CXX=g++
CXXFLAGS=-Wall -O3

$(EXE) : $(CPP)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(CPP)


clean :
	rm -rf $(EXE)
