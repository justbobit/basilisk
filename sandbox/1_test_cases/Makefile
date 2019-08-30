CFLAGS +=  -O2

include $(BASILISK)/Makefile.defs

moving_embed.tst  	: CFLAGS += -fopenmp

mini_cell.tst  		: CFLAGS +=

test_advection.tst	: CFLAGS += -fopenmp

update.tst        	: CFLAGS += -fopenmp

reversed.tst      	: CFLAGS += -fopenmp

solidification_Gibbs_Thomson.tst : CFLAGS += -events

dirichlet.tst: CFLAGS += -DDIRICHLET=1 -DMTRACE=3
dirichlet.c: neumann_double.c
	ln -s neumann_double.c dirichlet.c
