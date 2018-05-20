#!/bin/sh
#Get the absolute path of the working shell file.
cd `dirname $0`

#compile the cpp file 
g++ -c Boltzmann.cpp 
ar rcs libBoltzmann.a Boltzmann.o
g++ -o test main.cpp -L. -lBoltzmann -lnrutil

#Determine whether the directory exists.
if [ ! -d "./data" ];then
	mkdir ./data
	echo "The new directory <data> has been set up."
else
	echo "The directory <data> exists."
fi

#run the executable file.
cp test ./data 
rm test
cd ./data
./test 



