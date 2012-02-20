SCATLIB = ../scattering
HYDROLIB = ../hydrometeor
BOOSTINC = /usr/include/boost-1_38

CXX      = LANG=C icpc
CXXFLAGS = -O3 -w1 -openmp -par-report1 -vec-report0 -xHost -axSSE4.1 \
           -ip -fPIC -ftz -fp-model precise
#CXXFLAGS = -g -openmp
CPPFLAGS = -DNDEBUG -I. -I$(HYDROLIB) -I$(SCATLIB) -I$(BOOSTINC)
AR       = ar
ARFLAGS  = -rv

LIBRADAR = libradar.a

OBJS    = polarimetric_radar.o radarband.o radargeom.o

TESTLIBS = -L. -lradar -L$(SCATLIB) -lscattering \
					 -lifcore -L$(HYDROLIB) -lhydro

all:: $(LIBRADAR) test_radar

$(LIBRADAR): $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)

clean:
	rm -f $(LIBRADAR) $(OBJS) test_radar

radarband.o: radarband.h
radargeom.o: radargeom.h
polarimetric_radar.o: polarimetric_radar.h $(SCATLIB)/tmatrix.h \
	$(SCATLIB)/particle.h $(SCATLIB)/radiation.h \
	$(HYDROLIB)/hydropsd.h $(HYDROLIB)/hydrometeor.h

test_radar: $(LIBRADAR) $(SCATLIB)/libscattering.a $(HYDROLIB)/libhydro.a
	$(CXX) -i-dynamic $(CPPFLAGS) $(CXXFLAGS) -o $@ -DTESTME \
		polarimetric_radar.cpp $(TESTLIBS)
