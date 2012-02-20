#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

#include <vol.h>
#include <antenna.h>
#include <netcdf.h>

static int ncstatus;

#define NOECHO (BADVAL-5)
#define BADVAL (float)0500

#define handle_error(a) if ((ncstatus = (a)) != NC_NOERR) { \
    fprintf(stderr, "NetCDF error on line %d: %s\n", __LINE__, \
      nc_strerror(ncstatus)); \
    exit(-1); \
}

static void add_string_attribute(int ncid, int fieldid,
                                 char *attname, char *value)
{
  if (value) {
    handle_error(nc_put_att_text(ncid, fieldid, attname,
                                 strlen(value), value));
  }
  return;
}

static void add_attributes(int ncid, int fieldid, char *units,
	                   char *long_name, char *standard_name)
{
  if (units) {
    add_string_attribute(ncid, fieldid, "units", units);
  }
  if (long_name) {
    add_string_attribute(ncid, fieldid, "long_name", long_name);
  }
  if (standard_name) {
    add_string_attribute(ncid, fieldid, "standard_name", standard_name);
  }
  return;
}

void add_missing_value(int ncid, int fieldid, float value)
{
  handle_error(nc_put_att_float(ncid, fieldid,
				"missing_value",
				NC_FLOAT, 1, &value));
  handle_error(nc_put_att_float(ncid, fieldid,
				"_FillValue",
				NC_FLOAT, 1, &value));
  return;
}

int main(int argc, char *argv[])
{
  char *edgefile;
  int ncid;
  int volid, sweepid, rayid, binid, time_strlenid;
  int ivol, isweep, ifield, iray, ibin;
  char *voldata = NULL;
  struct vol_struct *volume = NULL;
  time_t xtime;
  struct tm *voltime, *sweeptime;
  char tmpstring[128];
  char ncfile[PATH_MAX];
  int dimids[5]; 
  int max_range_id, lonid, latid, altid, pulsid, prfid, waveid, nyqid;
  int gatesid, aspdid, urangeid, voltid, sweeptid;
  int elevid, azimid;
  float tmp_float;
  unsigned short *sray;
  unsigned char *ucray;
  int iarg, nfields, nsweeps, nrays, nbins;
  int dzid, czid, rvid, swid, zdid;
  float missing;
  float *values;
  size_t start[5];
  size_t count[5];
  float uz,cz,rv,sw,zdr = 0.0,nyq;
  int bytes_per_bin = 0;
  int nvols;

  if (argc < 2)
  {
    fprintf(stderr, "Not enough arguments.\n");
    fprintf(stderr, "Usage: %s volfile [volfile2 volfile3 ...]\n", argv[0]);
    return -1;
  }

  fprintf(stdout, "Doing first pass of %d volumes ", argc-1);

  nvols = 0;
  nsweeps = 0;
  nbins = 0;  
  nrays = 0;
  nfields = 4;
  for (iarg = 1; iarg < argc; iarg ++)
  {
    int xns, xnr, xnb;
    edgefile = strdup(argv[iarg]);
    if (load_data(&voldata, edgefile, 0) == -1)
    {
      fprintf(stderr, "Could not load EDGE Volume File: %s\n", edgefile);
      continue;
    }
    nvols ++;
    volume = (struct vol_struct *) voldata;
    fprintf(stdout, ".");
    xns = volume->num_sweeps;
    xnr = volume->sweep[0].num_rays;
    for (isweep = 1; isweep < xns; isweep ++)
      xnr = xnr < volume->sweep[isweep].num_rays ?
               volume->sweep[isweep].num_rays : xnr;
    xnb = volume->sweep[0].rad.gates;
    for (isweep = 1; isweep < xns; isweep ++)
      xnb = xnb < volume->sweep[isweep].rad.gates ?
               volume->sweep[isweep].rad.gates : xnb;
    nsweeps = xns > nsweeps ? xns : nsweeps;
    nrays = xnr > nrays ? xnr : nrays;
    nbins = xnb > nbins ? xnb : nbins;

    bytes_per_bin = BYTES_BIN(volume);
    if (bytes_per_bin > nfields) nfields = bytes_per_bin;

    if (iarg == 1)
    {
      xtime = volume->date;
      voltime = gmtime(&xtime);
      strftime(tmpstring, 128, "%Y-%m-%d_%H:%M:%S_%Z", voltime);
      sprintf(ncfile, "%s_%s_%s.nc", volume->sweep[0].rad.site_name,
              volume->sweep[0].rad.radar_type, tmpstring);
      if (nc_create(ncfile, NC_CLOBBER, &ncid) != NC_NOERR)
      {
        fprintf(stderr, "Could not create NetCDF Volume file: %s\n", ncfile);
        return -1;
      }
      add_string_attribute(ncid, NC_GLOBAL, "Conventions", "CF-1.0");
      add_string_attribute(ncid, NC_GLOBAL, "location",
                           volume->sweep[0].rad.site_name);
      add_string_attribute(ncid, NC_GLOBAL, "system",
                           volume->sweep[0].rad.radar_type);
      add_string_attribute(ncid, NC_GLOBAL, "title", "Radar data");
      add_string_attribute(ncid, NC_GLOBAL, "source", "CETEMPS");
      add_string_attribute(ncid, NC_GLOBAL, "references",
                           "http://www.himet.it");
      add_string_attribute(ncid, NC_GLOBAL, "scantype",
                           SCAN_TYPE(volume->scan_type));
      add_string_attribute(ncid, NC_GLOBAL, "job",
                         volume->sweep[0].rad.job_name);
      add_string_attribute(ncid, NC_GLOBAL, "host",
                           volume->sweep[0].rad.arch_host);
      add_string_attribute(ncid, NC_GLOBAL, "processing",
                           volume->sweep[0].rad.proc_description);
      tmpstring[0] = 0;
      if (volume->moment_enable|ARC_MOMENT) strncat(tmpstring, " A", 3);
      if (volume->moment_enable|U_MOMENT) strncat(tmpstring, " U", 3);
      if (volume->moment_enable|Z_MOMENT) strncat(tmpstring, " Z", 3);
      if (volume->moment_enable|V_MOMENT) strncat(tmpstring, " V", 3);
      if (volume->moment_enable|W_MOMENT) strncat(tmpstring, " W", 3);
      if (volume->moment_enable|ZDR_MOMENT) strncat(tmpstring, " D", 3);
      add_string_attribute(ncid, NC_GLOBAL, "moments_enabled", tmpstring);
      add_string_attribute(ncid, NC_GLOBAL, "history",
                           "Generated from scan data in EDGE format");
    }
    free(edgefile);
  }

  fprintf(stdout, "\n");
  fprintf(stdout, "Done first pass.\n");
  fprintf(stdout, "Sweeps %d, Rays %d, Bins %d\n", nsweeps, nrays, nbins);
  fprintf(stdout, "Start processing.\n");

  handle_error(nc_def_dim(ncid, "volume", nvols, &volid));
  handle_error(nc_def_dim(ncid, "sweep", nsweeps, &sweepid));
  handle_error(nc_def_dim(ncid, "ray", nrays, &rayid));
  handle_error(nc_def_dim(ncid, "bins", nbins, &binid));
  handle_error(nc_def_dim(ncid, "time_strlen", 32, &time_strlenid));

  dimids[0] = volid;
  handle_error(nc_def_var(ncid, "max_range", NC_FLOAT, 1, dimids,
               &max_range_id));
  add_attributes(ncid, max_range_id, "km", "Max Range", NULL);
  handle_error(nc_def_var(ncid, "latitude", NC_FLOAT, 1, dimids, &latid));
  add_attributes(ncid, latid, "degrees_north",
                 "Latitude of antenna", "latitude");
  handle_error(nc_def_var(ncid, "longitude", NC_FLOAT, 1, dimids, &lonid));
  add_attributes(ncid, lonid, "degrees_east",
                 "Longitude of antenna", "longitude");
  handle_error(nc_def_var(ncid, "altitude", NC_FLOAT, 1, dimids, &altid));
  add_attributes(ncid, altid, "m",
                 "Height of instrument above mean sea level", NULL);
  dimids[1] = time_strlenid;
  handle_error(nc_def_var(ncid, "volume_time", NC_CHAR, 2,
               dimids, &voltid));
  dimids[1] = sweepid;
  handle_error(nc_def_var(ncid, "pulse_width", NC_FLOAT, 2,
               dimids, &pulsid));
  add_attributes(ncid, pulsid, "ns", "Pulse length", NULL);
  handle_error(nc_def_var(ncid, "pulse_repetition", NC_FLOAT, 2,
               dimids, &prfid));
  add_attributes(ncid, prfid, "Hz", "Pulse repetition frequency", NULL);
  handle_error(nc_def_var(ncid, "wavelength", NC_FLOAT, 2,
               dimids, &waveid));
  add_attributes(ncid, waveid, "m", "Wavelength", NULL);
  handle_error(nc_def_var(ncid, "nyquist", NC_FLOAT, 2, dimids, &nyqid));
  add_attributes(ncid, nyqid, "m/s", "Nyquist Velocity", NULL);
  handle_error(nc_def_var(ncid, "gate_size", NC_FLOAT, 1,
               dimids, &gatesid));
  add_attributes(ncid, gatesid, "m", "Data gate size", NULL);
  handle_error(nc_def_var(ncid, "antenna_speed", NC_FLOAT, 2,
               dimids, &aspdid));
  add_attributes(ncid, aspdid, "degrees/sec", "Azimuth speed", NULL);
  handle_error(nc_def_var(ncid, "unam_range", NC_FLOAT, 2,
               dimids, &urangeid));
  add_attributes(ncid, urangeid, "km", "Unambiguous range", NULL);

  dimids[2] = time_strlenid;
  handle_error(nc_def_var(ncid, "sweep_time", NC_CHAR, 3,
               dimids, &sweeptid));
  dimids[2] = rayid;
  missing = -1;
  handle_error(nc_def_var(ncid, "elevation", NC_FLOAT, 3, dimids, &elevid));
  add_missing_value(ncid, elevid, missing);
  add_attributes(ncid, elevid, "degrees",
           "Elevation of the centre of the ray above the horizon", NULL);
  handle_error(nc_def_var(ncid, "azimuth", NC_FLOAT, 3, dimids, &azimid));
  add_missing_value(ncid, azimid, missing);
  add_attributes(ncid, azimid, "degrees",
           "Azimuth of the centre of the ray clockwise from due north", NULL);
  dimids[3] = binid;
  missing = NOECHO;
  handle_error(nc_def_var(ncid, "dz", NC_FLOAT, 4, dimids, &dzid));
  add_attributes(ncid, dzid, "dB", "Uncorrected Reflectivity", NULL);
  add_missing_value(ncid, dzid, missing);
  handle_error(nc_def_var(ncid, "cz", NC_FLOAT, 4, dimids, &czid));
  add_attributes(ncid, czid, "dB", "Corrected Reflectivity", NULL);
  add_missing_value(ncid, czid, missing);
  handle_error(nc_def_var(ncid, "rv", NC_FLOAT, 4, dimids, &rvid));
  add_attributes(ncid, rvid, "m/s", "Radial Velocity", NULL);
  add_missing_value(ncid, rvid, missing);
  handle_error(nc_def_var(ncid, "sw", NC_FLOAT, 4, dimids, &swid));
  add_attributes(ncid, swid, "mm", "Spectrum Width", NULL);
  add_missing_value(ncid, swid, missing);
  if (nfields > 4)
  {
    handle_error(nc_def_var(ncid, "zdr", NC_FLOAT, 4, dimids, &zdid));
    add_attributes(ncid, zdid, "mm", "Differential Reflectivity", NULL);
    add_missing_value(ncid, zdid, missing);
  }

  handle_error(nc_enddef(ncid));

  ivol = 0;
  for (iarg = 1; iarg < argc; iarg ++)
  {
    edgefile = strdup(argv[iarg]);
    if (load_data(&voldata, edgefile, 0) == -1)
    {
      fprintf(stderr, "Could not load EDGE Volume File: %s\n", edgefile);
      continue;
    }

    volume = (struct vol_struct *) voldata;

    bytes_per_bin = BYTES_BIN(volume);
    xtime = volume->date;
    voltime = gmtime(&xtime);
    strftime(tmpstring, 128, "%Y-%m-%d %H:%M:%S %Z", voltime);
    fprintf(stdout, "Volume time: %s\n", tmpstring);
    start[0] = ivol;
    count[0] = 1;
    handle_error(nc_put_vara_float(ncid, max_range_id, start, count,
                    &(volume->sweep[0].max_range)));
    tmp_float = volume->sweep[0].rad.lat_deg;
    tmp_float = tmp_float + 
               copysign(tmp_float, volume->sweep[0].rad.lat_min) / 60.0 + 
               copysign(tmp_float, volume->sweep[0].rad.lat_sec) / 3600.0;
    handle_error(nc_put_vara_float(ncid, latid, start, count, &tmp_float));
    tmp_float = volume->sweep[0].rad.long_deg;
    tmp_float = tmp_float +
               copysign(tmp_float, volume->sweep[0].rad.long_min) / 60.0 +
               copysign(tmp_float, volume->sweep[0].rad.long_sec) / 3600.0;
    handle_error(nc_put_vara_float(ncid, lonid, start, count, &tmp_float));
    tmp_float = volume->sweep[0].rad.antenna_height;
    handle_error(nc_put_vara_float(ncid, altid, start, count, &tmp_float));

    start[1] = 0;
    count[1] = strlen(tmpstring) + 1;
    handle_error(nc_put_vara_text(ncid, voltid, start, count, tmpstring));

    values = (float *) malloc(nrays*nsweeps*sizeof(float));
    if (values == NULL)
    {
      fprintf(stderr, "Malloc error: %ld bytes\n", nrays*nsweeps*sizeof(float));
      return 1;
    }

    for (isweep = 0; isweep < nsweeps; isweep ++)
      for (iray = 0; iray < nrays; iray ++)
        values[isweep*nrays+iray] = -1;
    for (isweep = 0; isweep < nsweeps; isweep ++)
    {
      for (iray = 0; iray < volume->sweep[isweep].num_rays; iray ++)
      {
        sray = (unsigned short *) RAY_PTR(volume,isweep,iray);
        tmp_float = (BIN2IANG(sray[0]) + BIN2IANG(sray[2])) / 2.0;
        if (fabs(BIN2IANG(sray[2]) - BIN2IANG(sray[0])) > 180.0)
          tmp_float -= 180.0;
        if (tmp_float < 0.0) tmp_float += 360.0;
        values[isweep*nrays+iray] = tmp_float;
      }
    }

    start[0] = ivol;
    start[1] = 0;
    start[2] = 0;
    count[0] = 1;
    count[1] = nsweeps;
    count[2] = nrays;
    handle_error(nc_put_vara_float(ncid, azimid, start, count, values));

    for (isweep = 0; isweep < nsweeps; isweep ++)
      for (iray = 0; iray < nrays; iray ++)
        values[isweep*nrays+iray] = -1;
    for (isweep = 0; isweep < nsweeps; isweep ++)
    {
      for (iray = 0; iray < volume->sweep[isweep].num_rays; iray ++)
      {
        sray = (unsigned short *) RAY_PTR(volume,isweep,iray);
        values[isweep*nrays+iray] = ((float)(BINEL2IANG100(sray[1])) +
                  (float)(BINEL2IANG100(sray[3])))/200.0;
      }
    }

    handle_error(nc_put_vara_float(ncid, elevid, start, count, values));
    free(values);

    values = (float *) malloc(nfields*nrays*nbins*sizeof(float));
    if (values == NULL)
    {
      fprintf(stderr, "Malloc error for %ld bytes\n",
              nfields*nrays*nbins*sizeof(float));
      return 1;
    }

    for (isweep = 0; isweep < nsweeps; isweep ++)
    {
      xtime = volume->sweep[isweep].date;
      sweeptime = gmtime(&xtime);
      strftime(tmpstring, 128, "%Y-%m-%d %H:%M:%S %Z", sweeptime);
      start[0] = ivol;
      start[1] = isweep;
      start[2] = 0;
      count[0] = 1;
      count[1] = 1;
      count[2] = strlen(tmpstring) + 1;
      handle_error(nc_put_vara_text(ncid, sweeptid, start, count, tmpstring));
      tmp_float =  volume->sweep[isweep].rad.pulse_width * 1200 + 800;
      handle_error(nc_put_vara_float(ncid, pulsid, start, count, &tmp_float));
      tmp_float = volume->sweep[isweep].rad.prf1;
      handle_error(nc_put_vara_float(ncid, prfid, start, count, &tmp_float));
      tmp_float = volume->sweep[isweep].rad.wavelength;
      handle_error(nc_put_vara_float(ncid, waveid, start, count, &tmp_float));
      tmp_float = volume->sweep[isweep].rad.prf1 *
                  volume->sweep[isweep].rad.wavelength / 4.0;
      nyq = tmp_float;
      handle_error(nc_put_vara_float(ncid, nyqid, start, count, &tmp_float));
      tmp_float = volume->sweep[isweep].rad.gw1;
      handle_error(nc_put_vara_float(ncid, gatesid, start, count, &tmp_float));
      tmp_float = volume->sweep[isweep].rad.antenna_speed*0.55;
      handle_error(nc_put_vara_float(ncid, aspdid, start, count, &tmp_float));
      tmp_float = 149851.274/volume->sweep[isweep].rad.prf1;
      handle_error(nc_put_vara_float(ncid, urangeid, start, count, &tmp_float));

      for (ifield = 0; ifield < nfields; ifield ++)
        for (iray = 0; iray < nrays; iray ++)
          for (ibin = 0; ibin < nbins; ibin ++)
            values[ifield*nrays*nbins+iray*nbins+ibin] = NOECHO;

      for (iray = 0; iray < volume->sweep[isweep].num_rays; iray ++)
      {
        ucray = (unsigned char *)RAY_PTR(volume,isweep,iray);
        for (ibin = 0; ibin < volume->sweep[isweep].rad.gates; ibin ++)
        {
          uz = (float) ucray[ibin*bytes_per_bin+2];
          if (uz == 0.0) uz = NOECHO;
          else if (uz > 255.0) uz = BADVAL;
          else uz = uz/2.0-32.0;
          cz =  (float) ucray[ibin*bytes_per_bin];
          if (cz == 0.0) cz = NOECHO;
          else if (cz > 255.0) cz = BADVAL;
          else cz = cz/2.0-32.0;
          rv = (float) ucray[ibin*bytes_per_bin+1];
          if (rv == 0.0) rv = NOECHO;
          else if (rv > 255.0) rv = BADVAL;
          else rv = (rv-128.0)/128.0*nyq;
          sw = (float) ucray[ibin*bytes_per_bin+3];
          if (sw == 0.0) sw = NOECHO;
          else if (sw > 255.0) sw = BADVAL;
          else sw = sw/128.0*nyq;
          if (nfields > 4)
          {
            zdr = (float) ucray[ibin*bytes_per_bin+4];
            if (zdr == 0.0) zdr = NOECHO;
            else if (zdr > 255.0) zdr = BADVAL;
            else zdr = (zdr-128.0)/16.0;
          }
          values[iray*nbins+ibin] = uz;
          values[iray*nbins+ibin+nrays*nbins] = cz;
          values[iray*nbins+ibin+nrays*nbins*2] = rv;
          values[iray*nbins+ibin+nrays*nbins*3] = sw;
          if (nfields > 4)
            values[iray*nbins+ibin+nrays*nbins*4] = zdr;
        }
      }
      start[0] = ivol;
      start[1] = isweep;
      start[2] = 0;
      start[3] = 0;
      count[0] = 1;
      count[1] = 1;
      count[2] = nrays;
      count[3] = nbins;
      handle_error(nc_put_vara_float(ncid, dzid, start, count, values));
      handle_error(nc_put_vara_float(ncid, czid, start, count,
                   values+nrays*nbins));
      handle_error(nc_put_vara_float(ncid, rvid, start, count,
                   values+2*nrays*nbins));
      handle_error(nc_put_vara_float(ncid, swid, start, count,
                   values+3*nrays*nbins));
      if (nfields > 4)
        handle_error(nc_put_vara_float(ncid, zdid, start, count,
                     values+4*nrays*nbins));
    }
    free(edgefile);
    ivol ++;
  }

  handle_error(nc_close(ncid));

  return 0;
}
