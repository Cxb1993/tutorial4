CC = mpicc
CFLAGS = -Wall -pedantic

SRC_SIM 	= 	helper 		\
  		    	init 		\
      			boundary_val 	\
		      	uvp 		\
			visual 		\
			parallel 	\
			sor 		\
		      	main

SRC_JM		=	helper		\
			join_matrix

all:
	@make sim
	@make join_matrix

sim: $(SRC_SIM:%=%.o)
	$(CC) $(CFLAGS) -o sim $(SRC_SIM:%=%.o)  -lm

join_matrix: $(SRC_JM:%=%.o)
	$(CC) $(CFLAGS) -o join_matrix $(SRC_JM:%=%.o)  -lm

%.o : %.c
	$(CC) -c $(CFLAGS) $*.c -o $*.o

clean:
	/bin/rm -f $(SRC_SIM:%=%.o) sim $(SRC_JM:%=%.o) join_matrix

helper.o     	: helper.h 
init.o       	: helper.h init.h 
boundary_val.o  : helper.h boundary_val.h 
uvp.o        	: helper.h uvp.h
visual.o     	: helper.h
parallel.o 	: parallel.h
join_matrix.o	: helper.h

$main.o       	: helper.h init.h boundary_val.h uvp.h visual.h parallel.h sor.h
$join_matrix.o	: helper.h
