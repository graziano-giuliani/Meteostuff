SCATLIB = ../scattering
BOOSTINC = /usr/include/boost-1_38

CXX      = LANG=C icpc
CXXFLAGS = -O3 -w1 -fPIC -openmp -par-report1 -vec-report0 -xHost -axSSE3 -ip \
           -ftz -fp-model precise
#CXXFLAGS = -g
CPPFLAGS = -DNDEBUG -I. -I$(SCATLIB) -I$(BOOSTINC)
AR       = ar
ARFLAGS  = -rv

LIBHYDRO = libhydro.a

OBJS    = hydrometeor.o  hydropsd.o  random.o  refindex.o dSFMT.o

TESTLIBS = -L. -lhydro

all:: $(LIBHYDRO) test_hydro

$(LIBHYDRO): $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)

clean:
	rm -f $(LIBHYDRO) $(OBJS) test_hydro

dSFMT.o: dSFMT.h dSFMT-params19937.h dSFMT-params.h
refindex.o: refindex.h
random.o: dSFMT.h random.h
hydropsd.o: hydropsd.h
hydrometeor.o: refindex.h hydrometeor.h \
               random.h hydropsd.h \
               $(SCATLIB)/particle.h

test_hydro: $(LIBHYDRO)
	$(CXX) -i-dynamic $(CPPFLAGS) $(CXXFLAGS) -o $@ -DTESTME hydrometeor.cpp $(TESTLIBS)
