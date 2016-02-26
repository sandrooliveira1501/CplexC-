CC=g++
IDIR=./include
ODIR=./obj
VPATH=./src

CONCERTDIR = /opt/ibm/ILOG/CPLEX_Studio1261/concert
CPLEXDIR = /opt/ibm/ILOG/CPLEX_Studio1261/cplex

TARGET=ilp_dist

#_DEPS = ilp.h prob.h
_DEPS = 
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

#_OBJ = main.o ilp.o prob.o
_OBJ = example.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

CFLAGS=-DIL_STD -O -DNDEBUG -I$(IDIR) -I$(CPLEXDIR)/include	\
-I$(CONCERTDIR)/include -fstrict-aliasing -pedantic -g -Wall		\
-fexceptions -ffloat-store -DILOUSEMT -D_REENTRANT -DILM_REENTRANT

LDFLAGS = -L$(CPLEXDIR)/lib/x86-64_linux/static_pic \
-lilocplex -lcplex -L$(CONCERTDIR)/lib/x86-64_linux/static_pic \
-lconcert -lm -lpthread

.PHONY: all clean

all: $(TARGET)

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LDFLAGS)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(ODIR)/*.o *~ $(IDIR)/*~ $(TARGET)
