#------------------------------------------------------------------------------

SOURCE=main.cpp
PROGRAM=SegAligner
CC=g++

#------------------------------------------------------------------------------

all: $(PROGRAM)

$(PROGRAM): $(SOURCE)

	$(CC) -Wall -std=c++11 $(SOURCE) -o $(PROGRAM) -O3

clean:

	rm -f $(PROGRAM).o
	rm -f $(PROGRAM)

