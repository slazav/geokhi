############################################################
ANIHOME = /home/sla/ani2D
include $(ANIHOME)/src/Rules.make

LIBS	= -L $(ANILIB) $(LIBAFT) $(LIBFEM) $(LIBLU) $(LIBLMR) $(LIBMBA) \
	              $(LIBVIEW) $(LIBLAPACK) $(LIBBLAS) -lm 


############################################################
all: main

OBJ  =  main.o solve_u.o solve_t.o solve_q.o mesh.o read_cfg.o


############################################################
main: $(OBJ) 
	$(FLINKER) $(LDFLAGS) -o $@ $(OBJ) $(LIBS)

############################################################
$(OBJ): th.fh


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f 
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c 
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o

clean:
	make -C res clean
	rm -f *.o main