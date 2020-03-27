EXEC = dplus.x

OBJS = main.o numrec.o allvars.o io.o cosmology.o 

INCL = allvars.h proto.h Makefile

OPTIONS = $(OPT)

CC       = gcc
OPTIMIZE = -O3 
LIBS     = -lm  
CFLAGS   = $(OPTIONS) $(OPTIMIZE) 

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS) -o  $(EXEC)  

$(OBJS): $(INCL) 


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)



