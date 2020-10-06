
MOD=mod
OBJ=obj

## PGI Fortran:
FC=pgf90
FFLAGS= -module $(MOD) -O2 -acc -ta=nvidia,time -Minfo=acc
#FFLAGS= -module $(MOD) -O2
LIBS=

## GNU Fortran
#FC=gfortran
#FFLAGS=-J$(MOD) -O2
#LIBS=

.SUFFIXES: .F90 .o 

# These are the object files (compiled from the F90 files) that are needed to
# create the execs.
SHEL_OBJS=$(addprefix $(OBJ)/, constants.o w3servmd.o w3gsrumd.o w3arrymd.o \
	w3gdatmd.o w3odatmd.o w3idatmd.o w3timemd.o w3dispmd.o w3adatmd.o \
	wmmdatmd.o w3src4md.o w3triamd.o w3snl1md.o w3parall.o \
	w3wdatmd.o w3iogrmd.o w3iorsmd.o w3nmlshelmd.o w3iopomd.o \
	w3iogomd.o w3fldsmd.o w3partmd.o w3cspcmd.o w3sln1md.o \
	w3iotrmd.o w3wdasmd.o w3updtmd.o w3srcemd.o w3iobcmd.o \
	w3iosfmd.o w3wavemd.o w3initmd.o ww3_shel.o)

GRID_OBJS=$(addprefix $(OBJ)/, constants.o w3servmd.o w3gsrumd.o w3arrymd.o \
	w3gdatmd.o w3odatmd.o w3idatmd.o w3timemd.o w3nmlgridmd.o w3dispmd.o \
	w3adatmd.o w3snl1md.o w3src4md.o w3iogrmd.o w3triamd.o ww3_grid.o)

OUTF_OBJS=$(addprefix $(OBJ)/, constants.o w3servmd.o w3timemd.o w3arrymd.o \
	w3dispmd.o w3gsrumd.o w3parall.o wmmdatmd.o w3gdatmd.o w3wdatmd.o \
	w3adatmd.o w3idatmd.o w3odatmd.o w3src4md.o w3snl1md.o w3iogrmd.o \
	w3iogomd.o ww3_outf.o)


.PHONY: all objdir moddir

# By default, build both the ww3_shel, ww3_grid and ww3_outf targets below:
all: ww3_shel ww3_grid ww3_outf

## Build and link the WW3_SHEL executable:
ww3_shel: objdir moddir $(SHEL_OBJS)
	$(FC) $(FFLAGS) $(LIBS) -o ww3_shel $(SHEL_OBJS)
	
## Build and link the WW3_GRID executable:
ww3_grid: objdir moddir $(GRID_OBJS)
	$(FC) $(FFLAGS) $(LIBS) -o ww3_grid $(GRID_OBJS)

## Build and link the WW3_GRID executable:
ww3_outf: objdir moddir $(OUTF_OBJS)
	$(FC) $(FFLAGS) $(LIBS) -o ww3_outf $(OUTF_OBJS)

objdir:
	@mkdir -p obj

moddir:
	@mkdir -p mod

	
## General rule for building objects from f90 source files:
## This takes an F90 file ($<) and compiles it into an object file ($@)
$(OBJ)/%.o:%.F90
	$(FC) $(FFLAGS) -c $< -o $@

## Rule to prevent make from itentifying Fortran .mod files as Modulas 2
%.o : %.mod

## Clean:
clean: 
	rm -f ww3_grid ww3_shel ww3_outf $(MOD)/*.mod $(OBJ)/* 


## Specific dependencies of each module:
$(OBJ)/constants.o : constants.F90 w3macros.h

$(OBJ)/w3adatmd.o: w3adatmd.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3gdatmd.o $(OBJ)/w3idatmd.o $(OBJ)/w3odatmd.o $(OBJ)/w3servmd.o

$(OBJ)/w3arrymd.o : w3arrymd.F90 w3macros.h $(OBJ)/w3servmd.o

$(OBJ)/w3cspcmd.o : w3cspcmd.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3servmd.o

$(OBJ)/w3dispmd.o : w3dispmd.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3gdatmd.o $(OBJ)/w3odatmd.o $(OBJ)/w3servmd.o

$(OBJ)/w3cspcmd.o : w3cspcmd.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3servmd.o

$(OBJ)/w3dispmd.o : w3dispmd.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3gdatmd.o $(OBJ)/w3odatmd.o $(OBJ)/w3servmd.o

$(OBJ)/w3fldsmd.o : w3fldsmd.F90 w3macros.h $(OBJ)/w3gsrumd.o $(OBJ)/w3idatmd.o $(OBJ)/w3timemd.o

$(OBJ)/w3gdatmd.o : w3gdatmd.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3gsrumd.o $(OBJ)/w3servmd.o

$(OBJ)/w3gsrumd.o : w3gsrumd.F90 w3macros.h $(OBJ)/w3servmd.o

$(OBJ)/w3idatmd.o : w3idatmd.F90 w3macros.h $(OBJ)/w3gdatmd.o $(OBJ)/w3servmd.o

$(OBJ)/w3initmd.o : w3initmd.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3adatmd.o $(OBJ)/w3arrymd.o $(OBJ)/w3dispmd.o $(OBJ)/w3gdatmd.o $(OBJ)/w3idatmd.o $(OBJ)/w3iogomd.o $(OBJ)/w3iogrmd.o $(OBJ)/w3iopomd.o $(OBJ)/w3iorsmd.o $(OBJ)/w3odatmd.o $(OBJ)/w3parall.o $(OBJ)/w3servmd.o $(OBJ)/w3timemd.o $(OBJ)/w3wdatmd.o

$(OBJ)/w3iobcmd.o : w3iobcmd.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3adatmd.o $(OBJ)/w3cspcmd.o $(OBJ)/w3gdatmd.o $(OBJ)/w3gsrumd.o $(OBJ)/w3odatmd.o $(OBJ)/w3servmd.o $(OBJ)/w3triamd.o $(OBJ)/w3wdatmd.o

$(OBJ)/w3iogomd.o : w3iogomd.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3adatmd.o $(OBJ)/w3dispmd.o $(OBJ)/w3gdatmd.o $(OBJ)/w3odatmd.o $(OBJ)/w3parall.o $(OBJ)/w3servmd.o $(OBJ)/w3wdatmd.o

$(OBJ)/w3iogrmd.o : w3iogrmd.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3dispmd.o $(OBJ)/w3gdatmd.o $(OBJ)/w3odatmd.o $(OBJ)/w3servmd.o $(OBJ)/w3snl1md.o $(OBJ)/w3src4md.o $(OBJ)/w3timemd.o

$(OBJ)/w3iopomd.o : w3iopomd.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3adatmd.o $(OBJ)/w3gdatmd.o $(OBJ)/w3gsrumd.o $(OBJ)/w3odatmd.o $(OBJ)/w3servmd.o $(OBJ)/w3triamd.o $(OBJ)/w3wdatmd.o

$(OBJ)/w3iorsmd.o : w3iorsmd.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3gdatmd.o $(OBJ)/w3odatmd.o $(OBJ)/w3parall.o $(OBJ)/w3servmd.o $(OBJ)/w3triamd.o $(OBJ)/w3wdatmd.o

$(OBJ)/w3iosfmd.o : w3iosfmd.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3adatmd.o $(OBJ)/w3gdatmd.o $(OBJ)/w3odatmd.o $(OBJ)/w3parall.o $(OBJ)/w3partmd.o $(OBJ)/w3servmd.o $(OBJ)/w3wdatmd.o

$(OBJ)/w3iotrmd.o : w3iotrmd.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3adatmd.o $(OBJ)/w3gdatmd.o $(OBJ)/w3gsrumd.o $(OBJ)/w3odatmd.o $(OBJ)/w3parall.o $(OBJ)/w3servmd.o $(OBJ)/w3timemd.o $(OBJ)/w3wdatmd.o

$(OBJ)/w3nmlgridmd.o : w3nmlgridmd.F90 w3macros.h  $(OBJ)/w3odatmd.o $(OBJ)/w3servmd.o

$(OBJ)/w3nmlshelmd.o : w3nmlshelmd.F90 w3macros.h  $(OBJ)/w3servmd.o $(OBJ)/wmmdatmd.o

$(OBJ)/w3odatmd.o : w3odatmd.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3gdatmd.o $(OBJ)/w3servmd.o

$(OBJ)/w3parall.o : w3parall.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3adatmd.o $(OBJ)/w3gdatmd.o $(OBJ)/w3idatmd.o $(OBJ)/w3odatmd.o

$(OBJ)/w3partmd.o : w3partmd.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3dispmd.o $(OBJ)/w3gdatmd.o $(OBJ)/w3odatmd.o

$(OBJ)/w3servmd.o : w3servmd.F90 w3macros.h $(OBJ)/constants.o

$(OBJ)/w3sln1md.o : w3sln1md.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3gdatmd.o $(OBJ)/w3odatmd.o $(OBJ)/w3servmd.o

$(OBJ)/w3snl1md.o : w3snl1md.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3adatmd.o $(OBJ)/w3gdatmd.o $(OBJ)/w3odatmd.o

$(OBJ)/w3src4md.o : w3src4md.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3dispmd.o $(OBJ)/w3gdatmd.o $(OBJ)/w3odatmd.o $(OBJ)/w3servmd.o

$(OBJ)/w3srcemd.o : w3srcemd.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3dispmd.o $(OBJ)/w3gdatmd.o $(OBJ)/w3idatmd.o $(OBJ)/w3odatmd.o $(OBJ)/w3sln1md.o $(OBJ)/w3snl1md.o $(OBJ)/w3src4md.o $(OBJ)/w3wdatmd.o

$(OBJ)/w3timemd.o : w3timemd.F90 w3macros.h $(OBJ)/w3odatmd.o $(OBJ)/w3servmd.o

$(OBJ)/w3triamd.o : w3triamd.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3adatmd.o $(OBJ)/w3gdatmd.o $(OBJ)/w3idatmd.o $(OBJ)/w3odatmd.o $(OBJ)/w3servmd.o

$(OBJ)/w3updtmd.o : w3updtmd.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3adatmd.o $(OBJ)/w3dispmd.o $(OBJ)/w3gdatmd.o $(OBJ)/w3idatmd.o $(OBJ)/w3odatmd.o $(OBJ)/w3parall.o $(OBJ)/w3servmd.o $(OBJ)/w3timemd.o $(OBJ)/w3triamd.o $(OBJ)/w3wdatmd.o

$(OBJ)/w3wavemd.o : w3wavemd.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3adatmd.o $(OBJ)/w3gdatmd.o $(OBJ)/w3idatmd.o $(OBJ)/w3iobcmd.o $(OBJ)/w3iogomd.o $(OBJ)/w3iogrmd.o $(OBJ)/w3iopomd.o $(OBJ)/w3iorsmd.o $(OBJ)/w3iosfmd.o $(OBJ)/w3iotrmd.o $(OBJ)/w3odatmd.o $(OBJ)/w3parall.o $(OBJ)/w3servmd.o $(OBJ)/w3srcemd.o $(OBJ)/w3timemd.o $(OBJ)/w3triamd.o $(OBJ)/w3updtmd.o $(OBJ)/w3wdatmd.o

$(OBJ)/w3wdasmd.o : w3wdasmd.F90 w3macros.h $(OBJ)/w3adatmd.o $(OBJ)/w3gdatmd.o $(OBJ)/w3odatmd.o $(OBJ)/w3wdatmd.o

$(OBJ)/w3wdatmd.o : w3wdatmd.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3gdatmd.o $(OBJ)/w3odatmd.o $(OBJ)/w3parall.o $(OBJ)/w3servmd.o

$(OBJ)/wmmdatmd.o : wmmdatmd.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3gdatmd.o $(OBJ)/w3odatmd.o $(OBJ)/w3servmd.o

$(OBJ)/ww3_grid.o : ww3_grid.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3arrymd.o $(OBJ)/w3dispmd.o $(OBJ)/w3gdatmd.o $(OBJ)/w3gsrumd.o $(OBJ)/w3iogrmd.o $(OBJ)/w3nmlgridmd.o $(OBJ)/w3odatmd.o $(OBJ)/w3servmd.o $(OBJ)/w3timemd.o $(OBJ)/w3triamd.o 

$(OBJ)/ww3_shel.o : ww3_shel.F90 w3macros.h $(OBJ)/w3adatmd.o $(OBJ)/w3fldsmd.o $(OBJ)/w3gdatmd.o $(OBJ)/w3idatmd.o $(OBJ)/w3initmd.o $(OBJ)/w3iogomd.o $(OBJ)/w3iogrmd.o $(OBJ)/w3iopomd.o $(OBJ)/w3nmlshelmd.o $(OBJ)/w3odatmd.o $(OBJ)/w3servmd.o $(OBJ)/w3timemd.o $(OBJ)/w3wavemd.o $(OBJ)/w3wdasmd.o $(OBJ)/w3wdatmd.o 

$(OBJ)/ww3_outf.o: ww3_outf.F90 w3macros.h $(OBJ)/constants.o $(OBJ)/w3adatmd.o $(OBJ)/w3arrymd.o $(OBJ)/w3gdatmd.o $(OBJ)/w3iogomd.o $(OBJ)/w3iogrmd.o $(OBJ)/w3odatmd.o $(OBJ)/w3servmd.o $(OBJ)/w3timemd.o $(OBJ)/w3wdatmd.o
