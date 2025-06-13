CC      = g++
CFLAGS  = -static -Wall -O3 -I. -I./eigen-3.3.9 -I./xdrfile
CXXFLAGS = $(CFLAGS) -std=c++17  

CPP_SRC = $(filter-out PyModule.cpp, $(wildcard *.cpp))
C_SRC   = $(wildcard xdrfile/*.c)
OBJ     = $(CPP_SRC:.cpp=.o) $(C_SRC:.c=.o)
EXE     = MS_dump

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $(CXXFLAGS) $^ -o $@

%.o: %.cpp
	$(CC) $(CXXFLAGS) -c $< -o $@

%.o: %.c
	gcc $(CFLAGS) -c $< -o $@

clean:
	rm -f $(EXE) *.o xdrfile/*.o

.PHONY: clean all
