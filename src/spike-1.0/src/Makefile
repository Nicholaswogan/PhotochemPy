include ./make.inc

ARCH?=x64
INSTALLDIR?=..

help:
	@echo "Usage: make ARCH=<arch> {all, clean}"
	@echo "--- <arch> is your architecture type and the name of the directory in $(INSTALLDIR)/lib. Defaults to x64"
	@echo "     where the SPIKE library will be created !!!! "
	@echo

ifeq ($(OPTION),2)
  FUNCTIONS_WRAPPER:=c
else
  FUNCTIONS_WRAPPER:=f90
endif

ifeq ($(ARCH),none)
  ARCH:=x64
endif

#==============================================
# Name of Libraries 
#==============================================

FF=$(F90)
FFFLAGS=$(F90FLAGS)
NAME=libspike.a

#sources Fortran90
CODE90_functions_wrapper = f90_functions_wrapper.f90
CODE90_aux = sclbprim.f90 dzlbprim.f90 spike_smp_utilities.f90 
CODE90 = dspike_smp.f90 sspike_smp.f90 cspike_smp.f90 zspike_smp.f90 
CODE77_aux = dpivul.f spivul.f cpivul.f zpivul.f 

#sources C
CODEC_functions_wrapper = c_functions_wrapper.c

ifeq ($(FUNCTIONS_WRAPPER),c)
	BASE_OBJECTS =  $(CODEC_functions_wrapper:.c=.o) $(CODE77_aux:.f=.o) $(CODE90_aux:.f90=.o)
	OBJECTS = $(CODE90:.f90=.o)  
else
	BASE_OBJECTS =  $(CODE90_functions_wrapper:.f90=.o) $(CODE77_aux:.f=.o) $(CODE90_aux:.f90=.o)
	OBJECTS = $(CODE90:.f90=.o)  
endif



####################################


# .SUFFIXES:
.SUFFIXES: .c .f .f90 .o
.PHONY: clean all

%.o: %.f90
	$(FF) $(FFFLAGS) $< -o $@

%.o: %.F90 
	$(FF) $(FFFLAGS) $< -o $@

%.o: %.f
	$(FF) $(FFFLAGS) $< -o $@

%.o: %.c
	$(CC) $(CFLAGS)  $< -o $@

$(OBJECTS):	$(BASE_OBJECTS)

all:dev 
	rm -f $(OBJECTS) $(BASE_OBJECTS)

dev:$(INSTALLDIR)/lib/$(ARCH)/$(NAME)

$(INSTALLDIR)/lib/$(ARCH)/$(NAME):$(OBJECTS)
	ar r $(NAME) $(OBJECTS) $(BASE_OBJECTS)
	@if test -d $(INSTALLDIR)/lib/$(ARCH); then \
	echo "directory $(ARCH) already exists"; \
	else echo "create directory $(ARCH)" && mkdir -p $(INSTALLDIR)/lib/$(ARCH); fi;
	mv -f $(NAME) $(INSTALLDIR)/lib/$(ARCH)

touch:
	@touch $(CODE90)

clean: 
	rm -f $(INSTALLDIR)/lib/$(ARCH)/$(NAME)
	rm -r -f $(INSTALLDIR)/lib/$(ARCH)
	rm -f $(OBJECTS) $(BASE_OBJECTS)



