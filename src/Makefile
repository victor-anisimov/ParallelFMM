EXE    ?= ../fmm.x
#FC      = nvfortran
#LD      = nvfortran
#FC      = ifx
#LD      = ifx
CC	= icx
FC      = mpif90
LD      = mpif90
#FC      = gfortran
#LD      = gfortran
#FC      = ftn
#LD      = ftn
#CC      = gcc-9
#FFLAGS  = -O3 -fPIE -fdefault-integer-8 -fdefault-real-8
#FFLAGS  = -O0 -g
#FFLAGS  = -O2 -s integer64 -default64
#FFLAGS  = -O2 -i8 -default64
#FFLAGS  = -O2 
#FFLAGS  = -O2 -h omp

CFLAGS	 = -O2 -g
#FFLAGS  = -O2 -fopenmp
#LDFLAGS = -fopenmp -o

# Cray PE
#FFLAGS  = -O2 -s integer64 -hpl=exe.pl -hdevelop
#LDFLAGS = -default64 -o

# NVHPC PE
#ifeq ($(FC), nvfortran)
ifeq ($(FC), ftn)
  # nVidia A100
  ACC = -DCPLR_NVFTN -DOFFL_GPU -gpu=cc80 -mp=gpu 
  #ACC = -DCPLR_NVFTN -DOFFL_CPU -mp=multicore
  FFLAGS  = -O2 -g -fPIE -Wall -Wextra $(ACC) $(COMPENV)
  LDFLAGS = -g -Wall -Wextra $(ACC) -o
  DEFINES = -DMPI
endif

# GNU PE
# -Wl,-rpath,$(HOME)/gcc/install/lib64
ifeq ($(FC), gfortran)
# FFLAGS  = -fdefault-integer-8
  FFLAGS  = -DCMPLR_GNU -DOFFL_CPU -O2 -g -fcheck=all -fopenmp -fPIE -fdefault-real-8 -Wall -Wextra \
			-Wno-conversion -Wno-unused-variable -Wno-unused-parameter -Wno-unused-dummy-argument \
			-Wno-maybe-uninitialized -fbounds-check $(COMPENV)
  LDFLAGS = -g -fcheck=all -fopenmp -fPIE -Wall -Wextra -Wconversion -o
endif
ifeq ($(FC), mpif90-gfortran)
  FFLAGS  = -DCPLR_GNU -DOFFL_CPU -O2 -g -fcheck=all -fopenmp -fPIE -fdefault-real-8 -Wall -Wextra \
			-Wno-conversion -Wno-unused-variable -Wno-unused-parameter -Wno-unused-dummy-argument \
			-Wno-maybe-uninitialized $(COMPENV)
  LDFLAGS = -g -fcheck=all -fopenmp -fPIE -Wall -Wextra -Wconversion -o
  DEFINES = -DMPI
endif


# Intel PE
#FFLAGS  = -O2 -g -qopenmp -fPIE -i4 -r8 $(COMPENV)
#LDFLAGS = -g -qopenmp -fPIE -o
#AOT
#FFLAGS  = -O2 -g -i4 -r8 -what -fiopenmp -fopenmp-targets=spir64_gen -Xopenmp-target-backend "-device pvc" $(COMPENV)
ifeq ($(FC), ifx)
  FFLAGS  = -DCPLR_IFX -DOFFL_CPU -O2 -g -i4 -r8 -what -check all -qopenmp $(COMPENV)
  LDFLAGS = $(FFLAGS) -o
endif
ifeq ($(FC), mpif90)
  # compile for GPU
  FFLAGS  = -fc=ifx -diag-disable 8291 -DCPLR_IFX -DOFFL_GPU -DUSE_DEVPTR -O2 -g -i4 -r8 -what -fiopenmp -fopenmp-targets=spir64 $(COMPENV)
  #FFLAGS  = -fc=ifx -mcmodel=medium -DCPLR_IFX -DOFFL_GPU -DUSE_DEVPTR -DUSE_ARRAY_REDUCTION -O2 -g -i4 -r8 -what -fiopenmp -fopenmp-targets=spir64 $(COMPENV)
  # compile for CPU 
  #FFLAGS  = -fc=ifx -mcmodel=medium -DCPLR_IFX -DOFFL_CPU -O2 -g -i4 -r8 -what -qopenmp $(COMPENV)
  LDFLAGS = $(FFLAGS) -lnuma -lpthread -lhwloc -o
  DEFINES = -DMPI
endif

##FFLAGS  = -O0 -g -fPIE -fdefault-integer-8 -fdefault-real-8 
##LDFLAGS = -Wl,--whole-archive -Wl,-z,muldefs -L /opt/cray/ddt/4.1.1.0_33583/lib/64 -ldmallocthcxx -Wl,--no-whole-archive -lstdc++ -static -fopenmp -o
##LDFLAGS = -static -fopenmp -o

#FFLAGS  = -O0 -g -i8 -mp=nonuma
#LDFLAGS = -default64 -o

LDLIBS  = 

#DEFINES = -DTIMER -DDEBUG -DDETAILED0 -DDYNAMICLIST -DMPI
DEFINES += -DTIMER -DDEBUG -DDETAILED0

OBJS	= hbmapi.o allocator.o common_constants.o common_data.o exit.o main.o readinp.o plan.o \
 prg.o expansion.o estimate.o energy.o timer.o utils.o rotate.o

.SUFFIXES: .o .F90 .h .c .mod

vpath % ../src/

%.mod: %.F90
	$(FC) -c $(FFLAGS) $(DEFINES) $<
.F90.o:
	$(FC) -c $(FFLAGS) $(DEFINES) $<
.c.o:
	$(CC) -c $(CFLAGS) $<

$(EXE): $(OBJS)
	@echo " "
	@echo "Linking $(EXE) ..."
	$(LD) $(LDFLAGS) $(EXE) $(OBJS) $(LDLIBS)
	@echo " "
	@ls -l $(EXE)
#	@strip $(EXE)

clean:
	rm -rf *.o *.mod $(EXE) $(EXE).dSYM ../obj

allocator.o:	allocator.F90
hbmapi.o:	hbmapi.c
plan.o:		common_data.o common_data.F90 allocator.o allocator.F90 hbmapi.o hbmapi.c
prg.o:		common_data.o common_data.F90 common_constants.o common_constants.F90 plan.o plan.F90 timer.o timer.F90
readinp.o:	common_data.o common_data.F90 timer.o timer.F90 allocator.o allocator.F90
energy.o:	common_data.o common_data.F90 timer.o timer.F90 plan.o plan.F90
#pbc.o:		common_data.o common_data.F90 timer.o timer.F90 plan.o plan.F90
utils.o:	common_data.o common_data.F90 timer.o timer.F90
expansion.o:	common_data.o common_data.F90 plan.o plan.F90
rotate.o:	common_data.o common_data.F90
main.o:		timer.o timer.F90
