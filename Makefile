CPPFLAGS = -O3
CC = g++

BIN_DIR = bin

GSLFLAGS_C = -I/opt/local/include 
GSLFLAGS_L = -L/opt/local/lib -lgsl -lgslcblas -lm 


TARGETS = $(BIN_DIR)/ftumch2 

all: clean $(TARGETS)

clean: 
	-rm -rf $(TARGETS)
	-rm -rf $(BIN_DIR)/*.o


utils.o: src/utils.cpp
	$(CC) $(CPPFLAGS) $(DFLAGS) $(GSLFLAGS_C) -c src/utils.cpp -o $(BIN_DIR)/utils.o

$(BIN_DIR)/ftumch2: src/ftumch2.cpp utils.o 
	$(CC) $(CPPFLAGS) $(DFLAGS) $(GSLFLAGS_C) $(GSLFLAGS_L)  src/ftumch2.cpp $(BIN_DIR)/utils.o  -o $(BIN_DIR)/ftumch2

