CFLAGS +=  -O2

include $(BASILISK)/Makefile.defs

anisotropy.tst  					: CFLAGS += -fopenmp

crystal.tst  						: CFLAGS += -fopenmp

dirichlet.tst 						: CFLAGS += -DDIRICHLET=1 -DMTRACE=3
	
dirichlet.c                         : neumann_double.c
					ln -s neumann_double.c dirichlet.c

gaussian_solidification.tst  		: CFLAGS += -fopenmp

moving_embed.tst  					: CFLAGS += -fopenmp

recons_speed.tst					: CFLAGS += -fopenmp

reversed.tst      					: CFLAGS += -fopenmp

solidification_Gibbs_Thomson.tst 	: CFLAGS += -events

solidification_mwe.tst 				: CFLAGS += -fopenmp

solidification_sine.tst				: CFLAGS += -fopenmp

test_advection.tst					: CFLAGS += -fopenmp

update.tst     					   	: CFLAGS += -fopenmp
