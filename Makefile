parKVFinder: dictionaryprocessing.o matrixprocessing.o pdbprocessing.o argparser.o tomlprocessing.o resultsprocessing.o move src/parKVFinder.c
	gcc -fopenmp -Isrc -g -o parKVFinder lib/dictionaryprocessing.o lib/matrixprocessing.o lib/pdbprocessing.o lib/argparser.o lib/tomlprocessing.o lib/resultsprocessing.o src/parKVFinder.c -lm
# 	if [ ! "$(KVFinder_PATH)" ]; then echo "\n\nKVFinder_PATH variable not found. Export KVFinder_PATH and save it to your ~/.bashrc."; echo "export KVFinder_PATH=`pwd`\n"; fi
# if [ -f ~/.bash_profile ]; then echo "export "; fi
# if [ -f ~/.bashrc ]; then echo "export "; fi

matrixprocessing.o: src/matrixprocessing.c src/matrixprocessing.h
	gcc -fopenmp -O3 -Isrc -c src/matrixprocessing.c -lm -static

dictionaryprocessing.o: src/dictionaryprocessing.c src/dictionaryprocessing.h
	gcc -Isrc -c src/dictionaryprocessing.c

pdbprocessing.o: src/pdbprocessing.c src/pdbprocessing.h
	gcc -Isrc -c src/pdbprocessing.c

argparser.o: src/argparser.c src/argparser.h
	gcc -Isrc -c src/argparser.c

tomlprocessing.o: src/tomlprocessing.c src/tomlprocessing.h
	gcc -Isrc -c src/tomlprocessing.c

resultsprocessing.o: src/resultsprocessing.c src/resultsprocessing.h
	gcc -Isrc -c src/resultsprocessing.c

move: dictionaryprocessing.o matrixprocessing.o pdbprocessing.o argparser.o tomlprocessing.o resultsprocessing.o
	mv dictionaryprocessing.o matrixprocessing.o pdbprocessing.o argparser.o tomlprocessing.o resultsprocessing.o lib/

link:
	cd ${BASH_SOURCE[0]%/*}
	sudo ln -s `pwd`/parKVFinder /usr/local/bin/parKVFinder

KVFinder_PATH:
	cd ${BASH_SOURCE[0]%/*}
	@echo >> ~/.bashrc
	@echo "# KVFinder environment variable" >> ~/.bashrc
	@echo export KVFinder_PATH=`pwd` >> ~/.bashrc

clean:
	rm -r lib/matrixprocessing.o lib/dictionaryprocessing.o lib/pdbprocessing.o lib/argparser.o lib/tomlprocessing.o lib/resultsprocessing.o parKVFinder
