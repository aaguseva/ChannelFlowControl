###################################################################
#           Makefile for  DNS CODE
###################################################################
#
########################################################## compiler
#F77 = f90
#F77= mpif77
F77 = mpiifort #mpif90
CMPI_LIB =
OMP_FLAGS =
#F77 = mpxlf
#
#################################################### compiler flags
#F77FLAGS = -O3 -c -qmaxmem=-1 -qcache=auto -q32\
	   -qarch=ppc970 -qtune=ppc970 -qfixed=72


#LFLAGS  = -O3 -qmaxmem=-1 -qcache=auto -q32\
	  -qarch=ppc970 -qtune=ppc970 -qfixed=72
LFLAGS = -limf -O3 -lm -heap-arrays -g -traceback #-check all -fpe0 -warn -debug extended

F77FLAGS = -c -O3 -heap-arrays -g -traceback #-check all -fpe0 -warn -debug extended

#CCFLAGS = -O3 -c -qmaxmem=-1 -qcache=auto -qaltivec \
#          -q64 -qarch=ppc970 -qtune=ppc970 -pg
#
#
#
#
########################################################### machine
#SYSTEM = ibm6000
##################################################### objects alpha
INIT = .
SRCDIR = $(INIT)
OBJ=$(INIT)
OBJDIR = $(OBJ)
CALCDIR =$(INIT)

OBJECTS = $(OBJ)/main.o\
$(OBJ)/cftsingle.o\
$(OBJ)/rftsingle.o\
$(OBJ)/fou3D.o\
$(OBJ)/cross.v2.o\
$(OBJ)/var.o\
$(OBJ)/stats.o\
$(OBJ)/io.o\
$(OBJ)/laps.v7.o\
$(OBJ)/change.allfac.o\
$(OBJ)/cfdiff8.v7.o\
$(OBJ)/modules.v7.o\
$(OBJ)/genexp.o

#
#################################################### compile
#
$(OBJDIR)/modules.v7.o : $(SRCDIR)/modules.v7.f90  $(SRCDIR)/ctes3D $(SRCDIR)/makefile
	@echo compiling modules.v7.f90
	@cd $(OBJDIR); $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/modules.v7.f90
#
$(OBJDIR)/main.o : $(SRCDIR)/main.f90  $(SRCDIR)/ctes3D $(SRCDIR)/makefile $(SRCDIR)/modules.v7.o
	@echo compiling main.f90
	@cd $(OBJDIR); $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/main.f90

$(OBJDIR)/cftsingle.o : $(SRCDIR)/cftsingle.f $(SRCDIR)/makefile
	@echo compiling cftsingle.f
	@cd $(OBJDIR); $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/cftsingle.f
#
$(OBJDIR)/rftsingle.o : $(SRCDIR)/rftsingle.f $(SRCDIR)/makefile
	@echo compiling rftsingle.f
	@cd $(OBJDIR); $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/rftsingle.f
#
$(OBJDIR)/fou3D.o : $(SRCDIR)/fou3D.f $(SRCDIR)/ctes3D $(SRCDIR)/makefile
	@echo compiling fou3D.f
	@cd $(OBJDIR); $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/fou3D.f
#
$(OBJDIR)/cross.v2.o : $(SRCDIR)/cross.v2.f90  $(SRCDIR)/ctes3D $(SRCDIR)/makefile
	@echo compiling cross.v2.f90
	@cd $(OBJDIR); $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/cross.v2.f90
#
$(OBJDIR)/var.o : $(SRCDIR)/var.f90  $(SRCDIR)/ctes3D $(SRCDIR)/makefile
	@echo compiling var.f90
	@cd $(OBJDIR); $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/var.f90
#
$(OBJDIR)/stats.o : $(SRCDIR)/stats.f90  $(SRCDIR)/ctes3D $(SRCDIR)/makefile
	@echo compiling stats.f90
	@cd $(OBJDIR); $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/stats.f90
#
$(OBJDIR)/io.o : $(SRCDIR)/io.f90  $(SRCDIR)/ctes3D $(SRCDIR)/makefile
	@echo compiling io.f90
	@cd $(OBJDIR); $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/io.f90
#
$(OBJDIR)/laps.v7.o : $(SRCDIR)/laps.v7.f  $(SRCDIR)/ctes3D $(SRCDIR)/makefile $(SRCDIR)/modules.v7.o
	@echo compiling laps.v7.f
	@cd $(OBJDIR); $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/laps.v7.f
#
$(OBJDIR)/change.allfac.o : $(SRCDIR)/change.allfac.f  $(SRCDIR)/ctes3D $(SRCDIR)/makefile
	@echo compiling change.allfac.f
	@cd $(OBJDIR); $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/change.allfac.f
#
$(OBJDIR)/cfdiff8.v7.o : $(SRCDIR)/cfdiff8.v7.f90  $(SRCDIR)/ctes3D $(SRCDIR)/makefile $(SRCDIR)/modules.v7.o
	@echo compiling cfdiff8.v7.f90
	@cd $(OBJDIR); $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/cfdiff8.v7.f90
#
$(OBJDIR)/genexp.o : $(SRCDIR)/genexp.f  $(SRCDIR)/ctes3D $(SRCDIR)/makefile
	@echo compiling genexp.f
	@cd $(OBJDIR); $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/genexp.f
#
############################################################ build
CTRL    :printmsgA $(OBJECTS)
	@echo Linking...
	@$(F77) -o $@ $(OBJECTS) $(LFLAGS)
	@echo CTRL built, congratulations.
#
########################################################## message
printmsgA :
	@echo
	@echo Building CTRL for $(SYSTEM)
	@echo Compiler flags : $(F77FLAGS)
	@echo Linker   flags : $(LFLAGS)
#
###################################################### end of file
clean	:
	rm *.o *.mod
