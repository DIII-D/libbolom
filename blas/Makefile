LIB=libblas.a
SLIB=libblas.sl 

CFLAGS	      =

#for shared library
FFLAGS	      =+O2  +E1 +z +e +Ovectorize	
all: $(SLIB) $(LIB)

#for archive library
#FFLAGS	      =+Oall  +E1 +z +e  +Ovectorize
#all: $(LIB) 


OBJS	      = \
    dblas1.o \
		dblas2.o\
		dblas3.o\
		sblas1.o\
		sblas2.o\
		sblas3.o\
		util.o


SRCS	      = dblas1.f \
		dblas2.f\
		dblas3.f\
		sblas1.f\
		sblas2.f\
		sblas3.f\
		util.f


$(LIB): $(OBJS)
		ar r $(LIB) $? 
	
$(SLIB): $(OBJS)
		ld -v -t -b -o $(SLIB) -L/usr/lib -lf $? 

