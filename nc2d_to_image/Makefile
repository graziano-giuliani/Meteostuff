progs=nc_to_image image_to_nc

all:: $(progs)

CXX=g++
CXXFLAGS=-O3
LIBS=-lnetcdf_c++ -lnetcdf -lMagick++ -lMagick

nc_to_image: nc_to_image.cpp
	$(CXX) $(CXXFLAGS) -o nc_to_image $< $(LIBS)
image_to_nc: image_to_nc.cpp
	$(CXX) $(CXXFLAGS) -o image_to_nc $< $(LIBS)
