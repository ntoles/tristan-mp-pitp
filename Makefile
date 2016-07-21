#
# Tristan makefile
#

SHELL=/bin/sh

cc= h5pcc
FC= h5pfc
LD= h5pfc 

#this is the user file used to initialize the problem and provide user boundary conditions
#change the name of the user file for your problem; it will look fore $USER.F90 file. 
# 
USER= user_weibel
USER= user_twostream
USER= user_shock

# Executable name
EXECUTABLE= tristan-mp2d

#Performance options
PERFORMANCE= -O3  
#-xhost -qopt-report -g -ipo

#Precompiler options for enabling different algorithms
CUSTOM= -DMPI -DHDF5 -DserIO -DtwoD 

# -DMPI use MPI 
# -DHDF5 use HDF5
# -DserIO use serial IO, all processors send chunks to rank 0 which saves into file
#         omitting this option will use parallel IO; performance may vary
# -DtwoD 2D code (ignorable z coordinate), omitting this defaults to 3D
# -DDEBUG when enabled produces a lot of debugging information

# undocumented options: do not use
# -Dfilter2  
# -DABSORB 
# -DLOCALSCRATCH 
# -DLOCALRESTART

# flags normaly used for C and Fortran compilers
CFLAGS= $(PERFORMANCE)
FFLAGS= $(CUSTOM) $(PERFORMANCE)

# Executable name
EXECUTABLE= tristan-mp2d

# Objects
OBJS= system.o systemf.o par.o inputparser.o inputparserf.o fparser.o globaldata.o aux.o communications.o fields.o fieldboundaries.o particles.o domain.o dynamic_domain.o $(USER).o output.o restart.o particles_movedeposit.o overload.o  initialize.o tristanmainloop.o tristan.o

###########################

all: $(EXECUTABLE)

%.o:%.c
	$(cc) $(CFLAGS) $(INCPATH) -c $^

%.o:%.F90
	$(FC) $(FFLAGS) $(INCPATH) -c $^

$(EXECUTABLE): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) 
	tar -cf code.tar *.c *.h *.F90 Makefile* input*

clean: 
	rm -f $(OBJS)
	rm -f *.mod
	rm -f $(EXECUTABLE)
	rm -f code.tar
