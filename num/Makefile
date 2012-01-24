############################################################
ANIHOME = /home/sla/ani2D
include $(ANIHOME)/src/Rules.make

LIBS	= -L$(ANILIB) $(LIBAFT) $(LIBFEM) $(LIBLU) $(LIBLMR) $(LIBMBA) \
	              $(LIBVIEW) $(LIBLAPACK) $(LIBBLAS) -lm 


############################################################
all: main

OBJ  =  main.o forlibfem.o mkcell.o


############################################################
main: $(OBJ) 
	@echo '$(FLINKER) -o $@' 
	@$(FLINKER) $(LDFLAGS) -o $@ $(OBJ) $(LIBS)


############################################################
#main.o      : main.f
#forlibaft.o : forlibaft.f
#forlibfem.o : forlibfem.f


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f 
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c 
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o

