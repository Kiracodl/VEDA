PROGS = az_ddaskr force_viewer amp_red_averaging  attard_tester

# intel fortran compiler is about 20% - 25% faster than g95.  
# however, no longer supported by nanohub (at least as of Aug 2020)
# options are g95 or gfortran.  have not done any speed testing between them.
# abitrarily picking gfortran and may reconsider later
F77 = gfortran 
#F77 = g95 
#F77 = ifort #intel fortran

ifeq ($(wildcard /etc/environ.sh),)
    #when this file does not exist, we are not compiling on nanohub
    LIB_DIR = -L/apps/rappture/lib   -L/opt/intel/mkl/lib/intel64/ 
    MATH = lapack
    NANOHUB_LIBS = -lz -lexpat
else
    #here we are on nanohub.  for whatever reason, when I moved from the default gfortran 4.7 to the newer gfortran 7.2, the expat and z libraries are no longer necessary (and don't work) but compiler can't find blas on its own
    RP_DIR = $(shell . /etc/environ.sh; use -e -r rappture ;  echo $${RAPPTURE_PATH})
    MKL_DIR = $(shell . /etc/environ.sh; use -e -r intel-mkl-14.0.0; echo $${MKL_LIB_DIR})
    LIB_DIR = -L$(RP_DIR)/lib -L$(MKL_DIR)
    MATH = lapack
    NANOHUB_LIBS = -L /usr/lib/libblas -lblas
endif

ifeq  '$(strip $(MATH))' 'lapack95mkl'
# intel lapack 95 multi threaded version.
#on my dual core laptop, this is about a 20% speedup for N_attard = 300,
#but actually slows it down for N=50 or N = 100.
#no changes are needed to the code at all.  just link with this library and it is parallelized automatically
#8/9/2020: having issues getting this to work with gfortran
  MATH_LIBS =    libmkl_lapack95_lp64.a  -Wl,--start-group libmkl_intel_lp64.a libmkl_intel_thread.a libmkl_core.a  -Wl,--end-group -openmp -lpthread -lm
  MATH_DEFS = -DHAVE_LAPACK95MKL
else ifeq  '$(strip $(MATH))' 'lapack95mkl_seq'
#intel lapack 95   single threaded version
# if we want to use intel math libs, here they are.
  MATH_LIBS =    libmkl_lapack95_lp64.a libmkl_intel_lp64.a -Wl,--start-group libmkl_core.a libmkl_sequential.a -Wl,--end-group
  MATH_DEFS = -DHAVE_LAPACK95mkl
else ifeq  '$(strip $(MATH))' 'lapack'
# generic lapack that comes with Debian. not the intel versions
  MATH_LIBS =  -L/usr/lib/x86_64-linux-gnu/lapack  -L/usr/lib/lapack  -llapack
  MATH_DEFS = -DHAVE_LAPACK
else
  MATH_LIBS =
  MATH_DEFS =
endif

LIBS = $(LIB_DIR) -lrappture -lrappture2  -lstdc++ -lsupc++ -lgcc_eh $(NANOHUB_LIBS)  $(MATH_LIBS)


ifeq '$(strip $(BUILD))' 'check'
  ifeq '$(strip $(F77))' 'gfortran'
#    OFLAGS =  -fcheck=all  -g -ffpe-trap=underflow,invalid,denormal
    OFLAGS =  -fcheck=all  -g 
  else
    OFLAGS =  -check bounds -check pointers -check uninit -traceback -g
  endif
  LIB_RAPPTURE = $(LIBS)
else ifeq '$(strip $(BUILD))' 'debug_singlethread'
  OFLAGS = -g	
  LIB_RAPPTURE = $(LIBS)
else ifeq '$(strip $(BUILD))' 'debug'
  OFLAGS = -g -fopenmp -DOPENMP
  LIB_RAPPTURE = $(LIBS)
else ifeq '$(strip $(BUILD))' 'noopt'
  OFLAGS = 
  LIB_RAPPTURE = $(LIBS)
else ifeq '$(strip $(BUILD))' 'prof'
  OFLAGS = -O3 -pg -fopenmp  # profiling
#  LIB_RAPPTURE = $(LIBS)  -lm_p
  LIB_RAPPTURE = $(LIBS)  
else ifeq '$(strip $(BUILD))' 'singlethread'
  OFLAGS = -O3
  LIB_RAPPTURE = $(LIBS)
else
  OFLAGS = -O3 -fopenmp  -DOPENMP
  LIB_RAPPTURE = $(LIBS)
endif


ifeq '$(strip $(F77))' 'gfortran'
  FFLAGS = -ffree-line-length-400 $(OFLAGS)  $(MATH_DEFS) -fimplicit-none    #gfortran
else ifeq '$(strip $(F77))' 'g95'
  FFLAGS = -ffree-line-length-huge -Wall $(OFLAGS)  $(MATH_DEFS) #g95
else
  FFLAGS = -implicitnone  $(OFLAGS) $(MATH_DEFS)  #ifort
endif


#separate flags for DDASKR library files
ifeq '$(strip $(BUILD))' 'prof'
  DDFFLAGS = -O3 -pg
else ifeq '$(strip $(BUILD))' 'noopt'
  DDFLAGS = -s
else
  DDFFLAGS = -O3 -s
endif

# clear out all suffixes then list only those we use
.SUFFIXES:
.SUFFIXES: .o .f .f90 .F90

ifeq '$(strip $(F77))' 'ifort'
ifeq ($(wildcard /etc/environ.sh),)
    #when this file does not exist, we are not compiling on nanohub, so do nothing
else
    #when this file does exists, we are compiling on nanohub, so use the right magic words
  F77 = . /etc/environ.sh; use -e -r intel-fc-14.0.0 ; ifort
endif
endif


ifeq '$(strip $(F77))' 'gfortran'
ifeq ($(wildcard /etc/environ.sh),)
    #when this file does not exist, we are not compiling on nanohub, so do nothing
else
    #when this file does exists, we are compiling on nanohub, so use the right magic words
    #the anaconda use is necessarily to get gfortran 7.2  without it, the default is gfortran 4.7.
    #For some unknown reason, peak force tool does not work correctly with the older compiler.  however, we do need some changes to the lirary calls     
    #there may be a more elegant way to do this instead of calling "use" everytime, but I couldnt find it.
    F77 = . /etc/environ.sh; use -e -r anaconda2-5.1 ; gfortran

endif
endif





all: $(PROGS) prony/prony_fit
	make --directory=../rappture install

install: all
	mv $(PROGS) ../bin
	make --directory=prony install

DMT_ddaskr.o: DMT_ddaskr.F90  data1.mod rappture_io.mod params.mod fft.mod 

forceviewer.o: modules.o basic_modules.o 

amp_red_averaging.o: amp_red_averaging.F90 amp_red.mod



data1.mod : modules.o

modules.o : params.mod modules.F90  lapack95.mod 

params.mod : basic_modules.o

rappture_io.mod : params.mod fft.mod rappture_io.o 

rappture_io.o: rappture_io.f90 

#this is just interface files.  always build it even if we don't use it
lapack95.mod : lapack.o
lapack.o : lapack.f90

fft.mod: fft.o
fft.o: fft.f90

amp_red.mod: modules_amp_red.o

modules_amp_red.o:  modules_amp_red.f90

attard_tester.o: attard_tester.f90

#generic rule for general files.  avoids needing specific rule for each file
.f90.o :
	$(F77) $(FFLAGS)  -c -o $@ $<  

.F90.o :
	$(F77) $(FFLAGS) -c -o $@ $<

#all the ddaskr library files have .f suffix, and we want to always use optimization for them
#anyway, so different rules
.f.o :
	$(F77) $(DDFFLAGS)  -c -o $@ $<

az_ddaskr: DMT_ddaskr.o daux.o ddaskr.o dlinpk.o rappture_io.o modules.o basic_modules.o fft.o
	$(F77) $(FFLAGS)  -o $@ $? $(LIB_RAPPTURE)

force_viewer: force_viewer.o rappture_io.o modules.o basic_modules.o fft.o
	$(F77) $(FFLAGS)  -o $@ $? $(LIB_RAPPTURE)

amp_red_averaging : amp_red_averaging.o rappture_io.o modules.o basic_modules.o modules_amp_red.o  fft.o
	$(F77) $(FFLAGS)  -o $@ $? $(LIB_RAPPTURE)

attard_tester: attard_tester.o modules.o basic_modules.o  rappture_io.o fft.o
	$(F77) $(FFLAGS)  -o $@ $?   $(LIB_RAPPTURE)

prony/prony_fit:
	make --directory=prony

clean:
	rm -f *.o *.mod $(PROGS)	
	make --directory=prony clean
	make --directory=../rappture clean

distclean: clean
	cd ../bin; rm -f az_ddaskr force_viewer
	cd ../bin; rm -f prony_fit prony_fit.ctf prony_fit.readme run_prony_fit.sh
	cd ../bin; rm -rf prony_fit_mcr

