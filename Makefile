EXEC = dplus.x

OBJS = main.o tools.o io.o numrec.o cosmology.o 

INCL = tools.h io.h numrec.h cosmology.h Makefile

OPTIONS = $(OPT)

CC       = gcc
OPTIMIZE = -O3 -march=native 
LIBS     = -lm  
CFLAGS   = $(OPTIONS) $(OPTIMIZE) 

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS) -o  $(EXEC)  

$(OBJS): $(INCL) 


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)



