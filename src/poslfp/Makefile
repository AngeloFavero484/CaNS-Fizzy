PROGRAM = distrlfp_shell

FC = gfortran 

FFLAGS := -g -O3 -fopenmp #-Wall -fcheck=all # -ffixed-line-length-none -fdefault-real-8

all: $(PROGRAM)

OBJS = lfp.o loadd.o param.o

SRCS = $(OBJS:.o=.f90)

$(PROGRAM): $(OBJS)
	mpif90 -o $@ $(OBJS) $(FFLAGS)

clean:
	rm -f core $(OBJS) *.mod $(PROGRAM) .depend

ifneq ($(MAKECMDGOALS),clean)
-include .depend
endif

.depend:
	./makedepo $(SRCS) > .depend

%.o: %.f90
	$(FC) -c -o $@ $(FFLAGS) $<

include .depend
