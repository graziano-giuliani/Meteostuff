LIBSCAT  = libscattering.a
BOOSTINC = /usr/include/boost-1_38

CXX      = LANG=C icpc
CXXFLAGS = -O3 -w1 -fPIC -openmp -par-report1 -vec-report0 -axSSE3 -xHost \
           -ip -ftz -fp-model precise
#CXXFLAGS = -g -openmp
CPPFLAGS = -DNDEBUG -I. -I$(BOOSTINC)
FC       = ifort
FFLAGS   = -O3 -openmp -par-report1 -vec-report0 -fPIC -xHost -axSSE3 \
					 -ip -ftz -fp-model precise
#FFLAGS   = -g
AR       = ar
ARFLAGS  = -rv

OBJS    = bessel.o gauss.o invert.o particle.o tmatrix.o lpd.o

TESTLIBS = -L. -lscattering -lifcore

all:: $(LIBSCAT) test_tmatrix

$(LIBSCAT): $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)

clean:
	rm -f $(LIBSCAT) $(OBJS) test_tmatrix

bessel.o: bessel.h
gauss.o: gauss.h
invert.o: invert.h
particle.o: gauss.h simplecomplex.h particle.h
tmatrix.o: simplecomplex.h gauss.h bessel.h invert.h \
           particle.h radiation.h tmatrix.h

test_tmatrix: $(LIBSCAT)
	$(CXX) -i-dynamic $(CPPFLAGS) $(CXXFLAGS) -o $@ -DTESTME tmatrix.cpp $(TESTLIBS)
