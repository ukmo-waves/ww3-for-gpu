
FC=pgf90
FFLAGS= -module obj/ -O2 -g -mp 
LIBS=

.SUFFIXES: .F90 .o 

# These are the object files (compiled from the F90 files) that are needed to
# create the WW3_SHEL executable. The order may be important here, as some modules
# rely on symbols from other modules. Although this should be handled by the
# dependencies.

SHEL_OBJS=$(addprefix obj/, constants.o w3servmd.o w3gsrumd.o w3arrymd.o w3gdatmd.o \
	w3odatmd.o w3idatmd.o w3timemd.o w3dispmd.o w3adatmd.o \
	wmmdatmd.o w3src4md.o w3triamd.o w3snl1md.o w3parall.o \
	w3wdatmd.o w3iogrmd.o w3iorsmd.o w3nmlshelmd.o w3iopomd.o \
	w3iogomd.o w3fldsmd.o w3partmd.o w3cspcmd.o w3sln1md.o \
	w3iotrmd.o w3wdasmd.o w3updtmd.o w3srcemd.o w3iobcmd.o \
	w3iosfmd.o w3wavemd.o w3initmd.o \
	ww3_shel.o)

GRID_OBJS=$(addprefix obj/, constants.o w3servmd.o w3gsrumd.o w3arrymd.o w3gdatmd.o \
	w3odatmd.o w3idatmd.o w3timemd.o w3nmlgridmd.o w3dispmd.o \
	w3adatmd.o w3snl1md.o w3src4md.o w3iogrmd.o \
	w3triamd.o ww3_grid.o)

# By default, build both the ww3_shel and ww3_grid targets below:
all: ww3_shel ww3_grid

## Build and link the WW3_SHEL executable:
ww3_shel: $(SHEL_OBJS)
	$(FC) $(FFLAGS) $(LIBS) -o ww3_shel $(SHEL_OBJS)
	
## Build and link the WW3_GRID executable:
ww3_grid: $(GRID_OBJS)
	$(FC) $(FFLAGS) $(LIBS) -o ww3_grid $(GRID_OBJS)
	
## General rule for building objects from f90 source files:
## This takes an F90 file ($<) and compiles it into an object file ($@)

obj/%.o:%.F90
	$(FC) $(FFLAGS) -c $< -o $@

## Rule to prevent make from itentifying Fortran .mod files as Modulas 2
%.o : %.mod

## Clean:
clean: 
	rm -f ww3_grid ww3_shel *.mod *.o obj/* 

## Specific dependencies of each module (to make building more robust 
## or enable parallel building), we could specify each dependency
## below, e.g.:

constants.o : constants.F90 w3macros.h

w3adatmd.o : w3adatmd.F90 w3macros.h constants.o w3gdatmd.o w3idatmd.o w3odatmd.o w3servmd.o

w3arrymd.o : w3arrymd.F90 w3macros.h w3servmd.o

w3cspcmd.o : w3cspcmd.F90 w3macros.h constants.o w3servmd.o

w3dispmd.o : w3dispmd.F90 w3macros.h constants.o w3gdatmd.o w3odatmd.o w3servmd.o

w3cspcmd.o : w3cspcmd.F90 w3macros.h constants.o w3servmd.o

w3dispmd.o : w3dispmd.F90 w3macros.h constants.o w3gdatmd.o w3odatmd.o w3servmd.o

w3fldsmd.o : w3fldsmd.F90 w3macros.h w3gsrumd.o w3idatmd.o w3timemd.o

w3gdatmd.o : w3gdatmd.F90 w3macros.h constants.o w3gsrumd.o w3servmd.o

w3gsrumd.o : w3gsrumd.F90 w3macros.h w3servmd.o

w3idatmd.o : w3idatmd.F90 w3macros.h w3gdatmd.o w3servmd.o

w3initmd.o : w3initmd.F90 w3macros.h constants.o w3adatmd.o w3arrymd.o w3dispmd.o w3gdatmd.o w3idatmd.o w3iogomd.o w3iogrmd.o w3iopomd.o w3iorsmd.o w3odatmd.o w3parall.o w3servmd.o w3timemd.o w3wdatmd.o

w3iobcmd.o : w3iobcmd.F90 w3macros.h constants.o w3adatmd.o w3cspcmd.o w3gdatmd.o w3gsrumd.o w3odatmd.o w3servmd.o w3triamd.o w3wdatmd.o

w3iogomd.o : w3iogomd.F90 w3macros.h constants.o w3adatmd.o w3dispmd.o w3gdatmd.o w3odatmd.o w3parall.o w3servmd.o w3wdatmd.o

w3iogrmd.o : w3iogrmd.F90 w3macros.h constants.o w3dispmd.o w3gdatmd.o w3odatmd.o w3servmd.o w3snl1md.o w3src4md.o w3timemd.o

w3iopomd.o : w3iopomd.F90 w3macros.h constants.o w3adatmd.o w3gdatmd.o w3gsrumd.o w3odatmd.o w3servmd.o w3triamd.o w3wdatmd.o

w3iorsmd.o : w3iorsmd.F90 w3macros.h constants.o w3gdatmd.o w3odatmd.o w3parall.o w3servmd.o w3triamd.o w3wdatmd.o

w3iosfmd.o : w3iosfmd.F90 w3macros.h constants.o w3adatmd.o w3gdatmd.o w3odatmd.o w3parall.o w3partmd.o w3servmd.o w3wdatmd.o

w3iotrmd.o : w3iotrmd.F90 w3macros.h constants.o w3adatmd.o w3gdatmd.o w3gsrumd.o w3odatmd.o w3parall.o w3servmd.o w3timemd.o w3wdatmd.o

w3nmlgridmd.o : w3nmlgridmd.F90 w3macros.h  w3odatmd.o w3servmd.o

w3nmlshelmd.o : w3nmlshelmd.F90 w3macros.h  w3servmd.o wmmdatmd.o

w3odatmd.o : w3odatmd.F90 w3macros.h constants.o w3gdatmd.o w3servmd.o

w3parall.o : w3parall.F90 w3macros.h constants.o w3adatmd.o w3gdatmd.o w3idatmd.o w3odatmd.o

w3partmd.o : w3partmd.F90 w3macros.h constants.o w3dispmd.o w3gdatmd.o w3odatmd.o

w3servmd.o : w3servmd.F90 w3macros.h constants.o

w3sln1md.o : w3sln1md.F90 w3macros.h constants.o w3gdatmd.o w3odatmd.o w3servmd.o

w3snl1md.o : w3snl1md.F90 w3macros.h constants.o w3adatmd.o w3gdatmd.o w3odatmd.o

w3src4md.o : w3src4md.F90 w3macros.h constants.o w3dispmd.o w3gdatmd.o w3odatmd.o w3servmd.o

w3srcemd.o : w3srcemd.F90 w3macros.h constants.o w3dispmd.o w3gdatmd.o w3idatmd.o w3odatmd.o w3sln1md.o w3snl1md.o w3src4md.o w3wdatmd.o

w3timemd.o : w3timemd.F90 w3macros.h w3odatmd.F90 w3servmd.F90

w3triamd.o : w3triamd.F90 w3macros.h constants.o w3adatmd.o w3gdatmd.o w3idatmd.o w3odatmd.o w3servmd.o

w3updtmd.o : w3updtmd.F90 w3macros.h constants.o w3adatmd.o w3dispmd.o w3gdatmd.o w3idatmd.o w3odatmd.o w3parall.o w3servmd.o w3timemd.o w3triamd.o w3wdatmd.o

w3wavemd.o : w3wavemd.F90 w3macros.h constants.o w3adatmd.o w3gdatmd.o w3idatmd.o w3iobcmd.o w3iogomd.o w3iogrmd.o w3iopomd.o w3iorsmd.o w3iosfmd.o w3iotrmd.o w3odatmd.o w3parall.o w3servmd.o w3srcemd.o w3timemd.o w3triamd.o w3updtmd.o w3wdatmd.o

w3wdasmd.o : w3wdasmd.F90 w3macros.h w3adatmd.o w3gdatmd.o w3odatmd.o w3wdatmd.o

w3wdatmd.o : w3wdatmd.F90 w3macros.h constants.o w3gdatmd.o w3odatmd.o w3parall.o w3servmd.o

wmmdatmd.o : wmmdatmd.F90 w3macros.h constants.o w3gdatmd.o w3odatmd.o w3servmd.o

ww3_grid.o : ww3_grid.F90 w3macros.h constants.o w3arrymd.o w3dispmd.o w3gdatmd.o w3gsrumd.o w3iogrmd.o w3nmlgridmd.o w3odatmd.o w3servmd.o w3timemd.o w3triamd.o 

ww3_shel.o : ww3_shel.F90 w3macros.h w3adatmd.o w3fldsmd.o w3gdatmd.o w3idatmd.o w3initmd.o w3iogomd.o w3iogrmd.o w3iopomd.o w3nmlshelmd.o w3odatmd.o w3servmd.o w3timemd.o w3wavemd.o w3wdasmd.o w3wdatmd.o 
