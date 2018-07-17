#!/bin/sh
#Get the absolute path of the working shell file.
cd `dirname $0`

#compile the cpp file 
make

#Determine whether the directory exists.
if [ ! -d "./data" ];then
	mkdir ./data
	echo "A new directory <data> has been set up."
else
	echo "The directory <data> exists."
fi

#run the executable file.
cp test ./data 
rm test
cd ./data
nohup ./test > test.log 2>&1 &



