all:
	g++  -O2  -o IBG-FISHR-LOWRAM-RIPOFF ./src/ErrorFinderMain.cpp ./src/ErrorFinderManager.cpp ./src/Consolidator.cpp ./src/ErrorCalculator.cpp