CC       := gcc 
LIBS     := -lm -lblas -llapack -Wall -lgsl -lgslcblas -fopenmp -O0 -ggdb -g
DEBUGGER := -O0 -ggdb
FLAGS    := -std=gnu99
TARGET   := ./../CG
OBJS     := cg.c.o functions.c.o draw.c.o io.c.o verlet.c.o aux.c.o

.SUFFIXES: .c .o  

.PHONY: clean 

%.c.o: %.c 
	$(CC) -c $< -o $@ $(LIBS) $(FLAGS)

$(TARGET): $(OBJS) 
	$(CC) -o $(TARGET) $(OBJS) $(LIBS) $(FLAGS)

all: $(TARGET) 

clean:
	rm *.o $(TARGET)
