PROG		= anacal
LIBDIR          = -L$(HOME)/lib 
INCLUDE		= -I$(HOME)/include
LIBS            = -lnr -lm 
CFLAGS          = -g -Dlinux -Wall $(DBXFLAGS) $(INCLUDE)
LDFLAGS         = $(LIBDIR) $(LIBS)

$(PROG):	$(PROG).o;	cc -o $@ $(PROG).o $(LDFLAGS)
clean:		;		rm -f *.o
