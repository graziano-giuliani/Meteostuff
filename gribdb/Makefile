DBFILE = ${HOME}/data/thiresias.db
NCDBFILE = ${HOME}/data/NCthiresias.db
prefix = ${HOME}/toolchain

CXX = g++
CPPFLAGS = -I. -I/db/forecast/toolchain/include/boost-1_42 \
           -I/db/forecast/toolchain/include -Wno-long-long -Wall -pedantic
CXXFLAGS = -g -O2 -mtune=native -fbounds-check

LIBS = -L/db/forecast/toolchain/lib -lboost_iostreams -lboost_date_time -lsqlite3

CODEC_LIBS = -ljasper -lpng
NETCDF_LIBS = -L/db/forecast/toolchain/lib -lnetcdf_c++ -lnetcdf

PROGRAMS = grib2ncfcst grib2db db2fcst place2db

OBJS_g2db = earthmap.o griblib.o
OBJS_db2f = suntime.o weather.o
OBJS_p2db = places.o earthmap.o

OBJS = $(OBJS_g2db) $(OBJS_db2f) $(OBJS_p2db)

all:: $(PROGRAMS)

grib2ncfcst: grib2ncfcst.cpp $(OBJS_g2db)
	$(CXX) -DDBFILE=\"$(NCDBFILE)\" $(CPPFLAGS) $(CXXFLAGS) -o $@ $< \
               $(OBJS_g2db) $(CODEC_LIBS) $(NETCDF_LIBS) $(LIBS)

grib2db: grib2db.cpp $(OBJS_g2db)
	$(CXX) -DDBFILE=\"$(DBFILE)\" $(CPPFLAGS) $(CXXFLAGS) -o $@ $< \
               $(OBJS_g2db) $(CODEC_LIBS) $(LIBS)

db2fcst: db2fcst.cpp $(OBJS_db2f)
	$(CXX) -DDBFILE=\"$(DBFILE)\" $(CPPFLAGS) $(CXXFLAGS) -o $@ $< \
               $(OBJS_db2f) $(LIBS)

place2db: place2db.cpp $(OBJS_p2db)
	$(CXX) -DDBFILE=\"$(DBFILE)\" $(CPPFLAGS) $(CXXFLAGS) -o $@ $< \
               $(OBJS_p2db) $(LIBS)

earthmap.o: earthmap.cpp earthmap.h
griblib.o: griblib.cpp griblib.h
suntime.o: suntime.cpp suntime.h
weather.o: weather.cpp weather.h
places.o: places.cpp places.h

clean:
	rm -f $(PROGRAMS) $(OBJS)

install:
	install -m 755 $(PROGRAMS) $(prefix)/bin

db: $(DBFILE) $(NCDBFILE)

cleandb:
	rm -f $(DBFILE) $(NCDBFILE)

$(DBFILE): thiresias.sql
	sqlite3 $(DBFILE) < thiresias.sql

$(NCDBFILE): grib2ncfcst.sql
	sqlite3 $(NCDBFILE) < grib2ncfcst.sql
