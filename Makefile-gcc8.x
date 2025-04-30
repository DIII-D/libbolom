

#   Make to build Linux version of the Bolometer library libbolom6565.linux.so and
#   make to build Linux64 version,                       libbolom6565.linux64.so
#                                                    and libbolm.so

# NOTE: 20060614 tbt Modified this Makefile to only do libbolom6565.so.
#                      libbolom3565.so is NOT default - see other Makefile.


#
# on Omega use module load gcc/8.x
#

# Note change CFLAGS to build either 3565 or 6565 version.

# 20150508 tbt Defined OSPATH for saturn
# 20150508 tbt On saturn changed /bin/pgcc to pgcc/
# 20150508 tbt Took out /usr/local/include to see if any include files are needed on saturn.
# 20121016 tbt From prad in gaprofiles getting:
#         /c/source/libbolom_64/libbolom6565.linux64.so: undefined symbol: __fss_s
#         Take out explicit libpgc.so... allow linker to find in LD_LIBRARY_PATH
# 20090123 tbt Updated Makefile to handle linux64 in addition to linux.
#              Makefile now recognizes if this is a 64 or 32-bit make
#              automatically based on OSPATH. pgf90 is default for type.
# 20070620 tbt Make LD pgf90.... note LD is defined twice:
#  LD=$(PGI_DIR)/bin/pgcc
#  LD	      = ld
# 20070620 tbt add libblas.a to resolve scopy_
# 20060609 tbt added -defit6565
# 20060607 tbt Change $(PGI) to $(PGI_DIR)
#------------------------------------------------------------------

OSPATH = linux64


DBFLAGS =  -O -w
LDBFLAGS =  


# Note change CFLAGS to build either 3565 or 6565 version.
#CFLAGS	      = -I. -I/usr/local/include $(DBFLAGS)  -DEFIT3365 -fPIC
#20150508 tbt
#CFLAGS	      = -I. -I/usr/local/include $(DBFLAGS)  -DEFIT6565 -fPIC
CFLAGS	      = -I.                      $(DBFLAGS)  -DEFIT6565 -fPIC

# 20150508
#FFLAGS	      = -I/usr/local/include $(DBFLAGS)  -fPIC
FFLAGS	      =                      $(DBFLAGS)  -fPIC -fallow-argument-mismatch # for gcc/11.x
FFLAGS	      =                      $(DBFLAGS)  -fPIC 

# 20060607 tbt
#PGI=/usr/local/pgi
#FC=$(PGI)/linux86/bin/pgf77
#CC=$(PGI)/linux86/bin/pgcc
#LD=$(PGI)/linux86/bin/pgcc

# 20150508 tbt
#FC=$(PGI_DIR)/bin/pgf77  -fPIC
#CC=$(PGI_DIR)/bin/pgcc   -fPIC
#LD=$(PGI_DIR)/bin/pgcc   -fPIC

# 20230927 whm for Omega

FC=gfortran  -fPIC
CC=gcc   -fPIC
LD=gfortran   -fPIC

DEST	      = .

# 20150508 tbt It looks as if EXTHDRS is not used, so I commented it out.
#EXTHDRS	      = /usr/local/include/nrutil.h

HDRS	      = bolom.h \
		boundext.h \
		bound3365.h \
		gmatrix3365.h \
		gpos3365.h \
		bound6565.h \
		gmatrix6565.h \
		gpos6565.h \
		volmat.h

INSTALL	      = /etc/install

# 20070620 tbt
#LD	      = ld
#LD	      = pgf90  -fPIC

# 20060607 tbt
#LDFLAGS      =  -L $(PGI)/linux86/lib -L /usr/local/lib  -lblas  \
                    $(LDBFLAGS)  -lm  -lpgf90 -lpgftnrtl -lpgc

# 20070620 tbt
#LDFLAGS	      =  -L $(PGI_DIR)/lib    \
#			$(PGI_DIR)/lib/libpgf90.a    \
#			$(PGI_DIR)/lib/libpgftnrtl.a \
#			$(PGI_DIR)/lib/libpgc.so      \
#                        $(PGI_DIR)/lib/libpgsse1.a   \
#			$(PGI_DIR)/lib/libpgsse2.a   \
#                 -L /usr/local/lib -lblas  \
#                 $(LDBFLAGS)  -lm  -lpgf90 -lpgftnrtl -lpgc


LDFLAGS	      =


CLIBS	      =  -L /usr/local/lib -ltoms -llinpack -llapack  \
                 -lrecipes  -lpppack -lm -lcl $(LDBFLAGS) 

# 20070620 tbt
#LIBS	      = 

# 20150508 tbt PGI_DIR not defined on saturn.
#LIBS	      =  $(PGI_DIR)/lib/libblas.a 

# 20230927 whm Changed for Omega
LIBS	      =  -L./blas -lblas
LIBS64	      =  -L./blas -lblas

# 20090122 tbt As per PGI, try this library instead of libblas.a
#              which is suppose to be shared... but doesn't exist.
#LIBS64	      =  $(PGI_DIR)/libso/libacml.so
#
# 20110414 acml.so doesn't exist on venus
#LIBS64	      =  $(PGI_DIR)/lib/libblas.a $(PGI_DIR)/libso/libpgc.so
<<<<<<< HEAD
#LIBS64	      =  -L $(LD_LIBRARY_PATH) libblas.a 
=======
LIBS64	      =  -L $(LD_LIBRARY_PATH) -L./blas -lblas
>>>>>>> origin/omega


MAKEFILE      = Makefile

OBJS	      = corespline.o \
		bolomlower.o \
		bolomlowercells.o \
		bolomlowercells5pt.o \
		bolomuppercells.o \
		bolomtotal.o \
		bolomupper.o \
		bolomproj.o \
		powell.o \
		idlwrappers.o\
		splinebasis.o\
		efitgeom.o\
		bolomsetgmatrix.o\
		nrutil.o \
		spline.o \
		tomscore.o \
		linmin.o \
		brent.o \
		mnbrak.o \
		f1dim.o \
		spimage_toms_pz.o

LOWEROBJS = bolomlowercells.o \
	    efitgeom.o

UPPEROBJS = bolomuppercells.o 

COREOBJS = idlwrappers.o \
       bolomsetgmatrix.o \
       corespline.o \
       splinebasis.o\
       bolomproj.o \
<<<<<<< HEAD
       lsei.o \
=======
       587.o \
>>>>>>> origin/omega
       tomscore.o

OBJS = $(COREOBJS) $(LOWEROBJS) $(UPPEROBJS)

PRINT	      = pr

PROGRAM       = bolomfit

SHELL	      = /bin/sh

SRCS	      = corespline.c \
		bolomlower.c \
		bolomupper.c \
		bolomtotal.c \
		bolomproj.c \
		powell.c \
		idlwrappers.c

SRCSF = splinebasis.f

# 20150608 tbt SYSHDRS does not seem to be used. Commented out
#SYSHDRS	      = /usr/include/errno.h \
#		/usr/include/machine/frame.h \
#		/usr/include/machine/frame.h \
#		/usr/include/machine/save_state.h \
#		/usr/include/machine/save_state.h \
#		/usr/include/math.h \
#		/usr/include/pwd.h \
#		/usr/include/signal.h \
#		/usr/include/stdio.h \
#		/usr/include/stdlib.h \
#		/usr/include/sys/errno.h \
#		/usr/include/sys/signal.h \
#		/usr/include/sys/stdsyms.h \
#		/usr/include/sys/syscall.h \
#		/usr/include/sys/types.h

<<<<<<< HEAD
all:		libbolom.so 
=======
all:		libbolom.so tests
>>>>>>> origin/omega

$(PROGRAM):     $(OBJS)  bolomfit.o $(OBJS)
		@echo "Linking $(PROGRAM) ..."
		$(CC) bolomfit.o  $(OBJS) $(CLIBS) -o $(PROGRAM)  
		chmod u+x bolomfit
		@echo "done"
bolomfitsl:     bolomfit.o
		@echo "Linking $(PROGRAM) ..."
		$(CC) bolomfit.o  -L. -lbolom $(CLIBS) -o $(PROGRAM)  
		chmod u+x bolomfit
		@echo "done"

<<<<<<< HEAD
libbolom.so:  $(OBJS) blas/libblas.a
	$(LD) -o libbolom.so -shared $(OBJS)  $(LDFLAGS) $(LIBS64);		\
=======
blas/libblas.a:
		cd blas && make
tests:
		cd test && make
libbolom.so:  $(OBJS) blas/libblas.a
	$(LD) -o libbolom.so -shared $(OBJS)  $(LDFLAGS) $(LIBS64);		\
	   mkdir -p linux64;                                                      \
>>>>>>> origin/omega
	   cp libbolom.so libbolom6565.linux64.so;				\
	   cp libbolom.so linux64/libbolom6565.linux64.so;			\

libbolomx.so: $(OBJS) bolomx.o
		ld -b $(OBJS) bolomx.o -o libbolomx.so  $(LDFLAGS) $(LIBS)

libbolom.a: $(OBJS) 
		ar r libbolom.a $?

clean:;		@rm -f $(OBJS) core

clobber:;	@rm -f $(OBJS) $(PROGRAM) core tags

depend:;	@mkmf -f $(MAKEFILE) ROOT=$(ROOT)

echo:;		@echo $(HDRS) $(SRCS)

index:;		@ctags -wx $(HDRS) $(SRCS)

install:	$(PROGRAM)
		@echo Installing $(PROGRAM) in $(DEST)
		@-strip $(PROGRAM)
		@if [ $(DEST) != . ]; then \
		(rm -f $(DEST)/$(PROGRAM); $(INSTALL) -f $(DEST) $(PROGRAM)); fi

print:;		@$(PRINT) $(HDRS) $(SRCS)

tags:           $(HDRS) $(SRCS); @ctags $(HDRS) $(SRCS)

update:		$(DEST)/$(PROGRAM)

###
bolomcore.o: bolomcore.c spimage_dierckx.c spimage_pppack.c $(HDRS)
bolomfit.o: $(HDRS)
bolomlower.o: bolomlower.c $(HDRS)
bolomupper.o: bolomupper.c $(HDRS)
bolomlowercells.o: bolomlowercells.c $(HDRS)
bolomlowercells5pt.o: bolomlowercells5pt.c $(HDRS)
bolomuppercells.o: bolomuppercells.c $(HDRS)
bolomtotal.o: bolomtotal.c $(HDRS)
bolomproj.o: bolomproj.c $(HDRS)
corespline.o: corespline.c $(HDRS)

