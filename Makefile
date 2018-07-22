all:./bin/test

./bin/test:./main/*.o ./lib/lib*.a 
	@g++ -o ./bin/test ./main/*.o ./lib/lib*.a

./main/*.o:./main/*.cpp 
	@cd main;g++ -c *.cpp -I../inc 

./lib/lib*.a:./src/*.cpp
	@cd ./src;g++ -c *.cpp -I../inc
	@cd ./src;for name in `ls *.o`;do libname=lib$${name%.*}.a;ar rcs ../lib/$${libname} $${name};done
	@rm -f ./src/*.o

clean:
	@rm -f */*.o */*.a bin/*

