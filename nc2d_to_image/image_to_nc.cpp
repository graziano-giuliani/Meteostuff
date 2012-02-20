#include <netcdfcpp.h>
#include <Magick++.h>
#include <iostream>

#define vartype short  // variable type
#define magictype "I", Magick::ShortPixel
char *ncvar = "Bathymetry"; // change to varname
char *xdim = "longitude";   // change to dimname
char *ydim = "latitude";    // change to dimname

int main(int argc, char *argv[])
{
  vartype *vals;
  vartype *mask;

  if (argc != 2)
  {
    std::cerr << "Need input NetCDF file name where "
      << ncvar << "[" << ydim << "," << xdim << "] is !" << std::endl;
    return -1;
  }

  char *ncfname = strdup(argv[1]);
  std::cout << "Searching variable " << ncvar << " in file " << ncfname
    << std::endl;

  NcFile ncf(ncfname, NcFile::ReadOnly);
  if (!ncf.is_valid())
  {
    std::cerr << "Cannot Open NetCDF input " << ncfname << std::endl;
    return -1;
  }

  long ny = ncf.get_dim(ydim)->size();
  long nx = ncf.get_dim(xdim)->size();

  if (nx <= 0 || ny <=0)
  {
    std::cerr << "Cannot find dimensions " << xdim << " and " << ydim
      << " !" << std::endl;
  }

  vals = new vartype[nx*ny];
  mask = new vartype[nx*ny];
  ncf.get_var(ncvar)->get(vals, ny, nx);
  ncf.close();

  std::cout << "Variable " << ncvar << " found !" << std::endl;
  std::cout << "Masking...";

  Magick::Image image;
  std::string oname;
  oname = ncvar;
  oname = oname + "_mask.pgm";
  image.read(oname);

  image.write(0, 0, nx, ny, magictype, mask);

  for (int i = 0; i < nx*ny; i ++)
    if (mask[i] == 0) vals[i] = 0;

  std::cout << "Done." << std::endl;

  NcFile ncf1(ncfname, NcFile::Write);
  if (!ncf1.is_valid())
  {
    std::cerr << "Cannot Open NetCDF output " << ncfname << std::endl;
    return -1;
  }

  ncf1.get_var(ncvar)->put(vals, ny, nx);

  delete [] vals;
  ncf1.close();

  return 0;
}
