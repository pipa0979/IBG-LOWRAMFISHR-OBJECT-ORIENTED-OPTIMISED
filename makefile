all:
	g++  -O2  -o IBG-FISHR-LOWRAM-TEST ./src/ErrorFinderMain.cpp ./src/ErrorFinderManager.cpp ./src/Consolidator.cpp ./src/ErrorCalculator.cpp
