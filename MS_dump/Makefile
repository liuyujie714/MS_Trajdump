CC  = g++ -static -fPIC -std=c++17
CFLAGS = -O3 -I. -I./eigen-3.3.9

SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)
EXE = MS_dump

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $(CFLAGS) $^ -o $@

%.o: %.cpp
	$(CC) -c $< $(CFLAGS) -o $@ 

clean:
	rm -f $(EXE) *.o

.PHONY: clean all
