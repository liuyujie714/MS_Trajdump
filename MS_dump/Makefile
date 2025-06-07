CC  = g++ -static -fPIC -std=c++17
src = $(wildcard *.c *.cpp)
OBJ = $(patsubst %.c,%.o,$(patsubst %.cpp,%.o,$(src)))

CLFAGS = -O3 -I.
EXE = MS_dump

all:$(EXE)

$(EXE):$(OBJ)
	$(CC) $(CLFAGS) $(OBJ) -o $(EXE) 

%.o:%.c
	$(CC) -c $< -o $@ $(CLFAGS) 

%.o:%.cpp
	$(CC) -c $< -o $@ $(CLFAGS) 

clean:
	rm -f $(EXE) *.o

.PHONY: clean all

