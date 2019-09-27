#------------------------------------------------------------------------------
CC=g++
OPTS=-std=c++11 -O3 -pthread
#------------------------------------------------------------------------------

all: main.o SegAligner.o SAHelper.o
	$(CC) $(OPTS) -lm -o SegAligner main.o SAHelper.o SegAligner.o
	rm -f *.o

main.o: main.cpp SAHelper.h SegAligner.h
	$(CC) $(OPTS) -c main.cpp

SegAligner.o: SegAligner.cpp SegAligner.h SAHelper.h
	$(CC) $(OPTS) -c SegAligner.cpp

SAHelper.o: SAHelper.cpp SAHelper.h
	$(CC) $(OPTS) -c SAHelper.cpp

clean:
	rm -f *.o SegAligner
