# Detect OS and set compiler
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S), Linux)
	CC := gcc
	CFLAGS := -fopenmp -O3
	LDFLAGS := -lm -fcommon
else ifeq ($(UNAME_S), Darwin)
	CC := clang
	CFLAGS := -Xpreprocessor -fopenmp=libomp -O3 -ffast-math
	LDFLAGS := -L/opt/homebrew/opt/libomp/lib -lomp
endif

# General flags
INCLUDES := -Isrc

build: utils.o fileprocessing.o gridprocessing.o argparser.o move src/parKVFinder.c requirements
	$(CC) $(CFLAGS) $(INCLUDES) -o parKVFinder lib/utils.o lib/fileprocessing.o lib/gridprocessing.o lib/argparser.o src/parKVFinder.c $(LDFLAGS)
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
	$(CC) $(CFLAGS) $(INCLUDES) -c src/utils.c $(LDFLAGS)

fileprocessing.o: src/fileprocessing.c src/fileprocessing.h utils.o
	$(CC) $(INCLUDES) -c src/fileprocessing.c $(LDFLAGS)

gridprocessing.o: src/gridprocessing.c src/gridprocessing.h
	$(CC) $(CFLAGS) $(INCLUDES) -c src/gridprocessing.c $(LDFLAGS)

argparser.o: src/argparser.c src/argparser.h
	$(CC) $(INCLUDES) -c src/argparser.c $(LDFLAGS)

move: utils.o fileprocessing.o gridprocessing.o argparser.o
	if [ ! -d "lib" ]; then mkdir lib/; fi
	mv utils.o fileprocessing.o gridprocessing.o argparser.o lib/

requirements: pip3

PIP3 := $(shell command -v pip3 2> /dev/null)

pip3:
ifndef PIP3
	@printf "\n> python3-pip is not available. To use parKVFinder with PyMOL v2, please install python3-pip and run make pip3\n\n"
else
	pip3 install -r tools/PyMOL2-parKVFinder-Tools/requirements.txt
endif

link:
	@if [ ! -d $(HOME)/.local/bin ]; then \
		mkdir -p $(HOME)/.local/bin; \
	fi
	@if [ -f $(HOME)/.local/bin/parKVFinder ]; then \
		printf "[==> parKVFinder symbolic link already exists ...\n"; \
	else \
		ln -s `pwd`/parKVFinder $(HOME)/.local/bin/parKVFinder; \
	fi

clean:
	if [ -d "lib" ]; then rm -r lib/; fi
	if [ -f parKVFinder ]; then rm parKVFinder; fi
