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
  ncf.get_var(ncvar)->get(vals, ny, nx);

  std::cout << "Variable " << ncvar << " found !" << std::endl;

  std::cout << "Rendering...";

  for (int i = 0; i < nx*ny; i ++)
    vals[i] = vals[i] >= 0 ? 0 : USHRT_MAX;

  Magick::Image *image = new Magick::Image(nx, ny, magictype, vals);

  std::string oname;
  oname = ncvar;
  oname = oname + ".pgm";
  image->write(oname.c_str());

  delete image;

  std::cout << "Done." << std::endl;

  delete [] vals;
  ncf.close();

  return 0;
}
