all:
	g++  -Wall -WExtra -Wno-sign-compare -g -o IBG-FISHR-LOWRAM-RIPOFF ./src/ErrorFinderMain.cpp ./src/ErrorFinderManager.cpp ./src/Consolidator.cpp ./src/ErrorCalculator.cpp