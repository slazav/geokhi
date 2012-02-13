############################################################
ANIHOME = /home/sla/ani2D
include $(ANIHOME)/src/Rules.make

LIBS	= -L $(ANILIB) $(LIBAFT) $(LIBFEM) $(LIBLU) $(LIBLMR) $(LIBMBA) \
	              $(LIBVIEW) $(LIBLAPACK) $(LIBBLAS) -lm 


############################################################
all: main

run: main
	./main

runi: run
	make -C png

OBJ  =  main.o solve_u.o solve_t.o solve_q.o solve_n.o solve_abcd.o\
         mesh.o read_cfg.o


############################################################
main: $(OBJ) 
	$(FLINKER) $(LDFLAGS) -o $@ $(OBJ) $(LIBS)

############################################################
$(OBJ): th.fh


############################################################
.f.o:
	$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	$(CC)  $(CFLAGS) -c $*.c -o $*.o

clean:
	make -C ps clean
	rm -f *.o main