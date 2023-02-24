parKVFinder: utils.o fileprocessing.o gridprocessing.o argparser.o move src/parKVFinder.c requirements
	gcc -fopenmp -Isrc -o parKVFinder lib/utils.o lib/fileprocessing.o lib/gridprocessing.o lib/argparser.o src/parKVFinder.c -lm -fcommon
	@if [ ! "${KVFinder_PATH}" ]; then \
		printf "\n\nKVFinder_PATH system variable not found. Export KVFinder_PATH to your system variables.\n"; \
		if [ -f ${HOME}/.bashrc ]; then \
			printf "Run the following command:\n"; \
			printf "> echo export KVFinder_PATH=`pwd` >> ~/.bashrc\n\n"; \
		elif [ -f ${HOME}/.bash_profile ]; then \
			printf "Run the following command:\n"; \
			printf "> echo export KVFinder_PATH=`pwd` >> ~/.bash_profile\n\n"; \
		else \
			printf "Set the following path to the KVFinder_PATH system variable in your configuration file:\n"; \
			printf "> KVFinder_PATH=`pwd`\n\n"; \
		fi \
	fi

utils.o: src/utils.c src/utils.h
	gcc -Isrc -c src/utils.c -fcommon

fileprocessing.o: src/fileprocessing.c src/fileprocessing.h utils.o
	gcc -Isrc -c src/fileprocessing.c -fcommon

gridprocessing.o: src/gridprocessing.c src/gridprocessing.h
	gcc -fopenmp -O3 -Isrc -c src/gridprocessing.c -lm -fcommon

argparser.o: src/argparser.c src/argparser.h
	gcc -Isrc -c src/argparser.c -fcommon

move: utils.o fileprocessing.o gridprocessing.o argparser.o
	if [ ! -d "lib" ]; then mkdir lib/; fi
	mv utils.o fileprocessing.o gridprocessing.o argparser.o lib/

requirements: pip pip3

PIP := $(shell command -v pip 2> /dev/null)
PIP3 := $(shell command -v pip3 2> /dev/null)

pip:
ifndef PIP
	printf "\n> python-pip is not available. To use parKVFinder with PyMOL v1, please install python-pip and run make pip\n\n"
else
	pip install -r tools/tk/requirements.txt
endif

pip3:
ifndef PIP3
	@printf "\n> python3-pip is not available. To use parKVFinder with PyMOL v2, please install python3-pip and run make pip3\n\n"
else
	pip3 install -r tools/PyMOL2-parKVFinder-Tools/requirements.txt
endif

link:
	@if [ -f /usr/local/bin/parKVFinder ]; then \
  		printf "[==> parKVFinder symbolic link already exist ...\n"; \
	else \
		sudo ln -s `pwd`/parKVFinder /usr/local/bin/parKVFinder; \
	fi

clean:
	if [ -d "lib" ]; then rm -r lib/; fi
	if [ -f parKVFinder ]; then rm parKVFinder; fi
