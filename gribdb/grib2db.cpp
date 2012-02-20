/***************************************************************************
 *   Copyright (C) 2008-2010 by Graziano Giuliani                          *
 *   graziano.giuliani at aquila.infn.it                                   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details. (see COPYING)            *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 *                                                                         *
 *   LIC: GPL                                                              *
 *                                                                         *
 ***************************************************************************/
#include <iostream>
#include <fstream>
#include <vector>

#include <sqlite3.h>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <float.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>

#include <earthmap.h>
#include <griblib.h>

using namespace himet;
using namespace boost::gregorian;
using namespace boost::posix_time;

void grib1_gds_to_emap(unsigned char *gds, char *pjdef, earthmap *m);
void grib1_pds_to_reftime(unsigned char *pds, char *reftime);
void grib1_pds_to_time(sqlite3 *db, unsigned char *pds, char *reftime,
                       char *tdef, long *tid, char *time1, char *time2);
void grib1_pds_to_level(sqlite3 *db, unsigned char *pds,
                        char *levelname, long *lev_id,
                        double *level1, double *level2);
void grib1_pds_to_g1_parameter(sqlite3 *db, unsigned char *pds,
                               char *vname, long *v_id);
void grib1_to_atmodata(sqlite3 *db,
                      t_g1_pinfo &pkci,
                      unsigned char *bds,
                      long scode_id, long vtid, long gid, long vid, long *mid);

void grib2_sec3_to_emap(unsigned char *secpnt, char *pjdef, earthmap *m);
void grib2_sec1_to_reftime(unsigned char *secpnt, char *reftime);
void grib2_sec4_to_time(sqlite3 *db, unsigned char *secpnt, char *reftime,
                        char *tdef, long *tid, char *time1, char *time2);
void grib2_sec4_to_level(sqlite3 *db, unsigned char *secpnt,
                         char *level, long *lev_id,
                         double *level1, double *level2);
void grib2_sec4_to_g2_parameter(sqlite3 *db, long discip, unsigned char *secpnt,
                                char *vname, long *v_id);
void grib2_sec5_to_pkcinfo(unsigned char *secpnt, t_g2_pinfo &pkci);
void grib2_to_atmodata(sqlite3 *db,
                      t_g2_pinfo &pkci,
                      unsigned char *sec7,
                      long scode_id, long vtid, long gid, long vid, long *mid);

int check_reftime(sqlite3 *db, char *reftime, long *reft_id);
int check_time(sqlite3 *db, long reft_id, long tid,
               char *time1, char *time2, long *valt_id);
int check_grid(sqlite3 *db, char *pjdef, long *grid, earthmap *m);
int check_level(sqlite3 *db, long ldes, double l1, double l2, long *lev_id);
int check_atmovar(sqlite3 *db, int gv, long g_id, long lev_id, long *atmo_id);

static const char gribdb[64] = DBFILE;
static char sqlbuff[2048];
char *sqlerr = 0;

bool verbose = false;

typedef struct {
  int gid;
  std::vector <long> pid;
  std::vector <long> ipos;
} t_gpoints;

std::vector <t_gpoints> gp;

int main(int argc, char *argv[])
{
  if (argc < 3)
  {
    std::cerr << "Need input grib file name as argument." << std::endl;
    std::cerr << "Need input grib file source code as argument." << std::endl;
    return -1;
  }

  long scode_id;
  sscanf(argv[2], "%ld", &scode_id);

  if (sqlite3_initialize() != SQLITE_OK)
  {
    std::cerr << "Database access impossible." << std::endl;
    return -1;
  }

  sqlite3 *ppDb;
  if (sqlite3_open(gribdb, &ppDb) != SQLITE_OK)
  {
    std::cerr << "Database file " << gribdb << " not ready." << std::endl;
    return -1;
  }

  sprintf(sqlbuff, "%s", "PRAGMA journal_mode = MEMORY;");
  if (sqlite3_exec(ppDb, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    return -1;
  }
  sqlite3_free(sqlerr);
  sprintf(sqlbuff, "%s", "PRAGMA cache_size = 32768;");
  if (sqlite3_exec(ppDb, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    return -1;
  }
  sqlite3_free(sqlerr);
  sprintf(sqlbuff, "%s", "PRAGMA page_size = 16384;");
  if (sqlite3_exec(ppDb, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    return -1;
  }
  sqlite3_free(sqlerr);
  sprintf(sqlbuff, "%s", "DROP INDEX IF EXISTS atmodata_placegrid_index;");
  if (sqlite3_exec(ppDb, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    return -1;
  }
  sqlite3_free(sqlerr);
  sprintf(sqlbuff, "%s", "DROP INDEX IF EXISTS atmodata_var_index;");
  if (sqlite3_exec(ppDb, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    return -1;
  }
  sqlite3_free(sqlerr);
  sprintf(sqlbuff, "%s", "DROP INDEX IF EXISTS atmodata_vtime_index;");
  if (sqlite3_exec(ppDb, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    return -1;
  }
  sqlite3_free(sqlerr);
  sprintf(sqlbuff, "%s", "DROP INDEX IF EXISTS atmodata_source_index;");
  if (sqlite3_exec(ppDb, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    return -1;
  }
  sqlite3_free(sqlerr);
  sprintf(sqlbuff, "%s", "CREATE TEMP TABLE tmpins(a integer , b integer ,"
                   "c integer, d integer, f double);");
  if (sqlite3_exec(ppDb, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    return -1;
  }
  sqlite3_free(sqlerr);

  struct stat finfo;
  if (stat(argv[1], &finfo) < 0)
  {
    std::cerr << "GRIB file source " << argv[1] << " not ready." << std::endl;
    return -1;
  }
  size_t size = finfo.st_size;
  size_t limit = size-4;
  size_t pos = 0;

  int fhandle = open(argv[1], O_RDONLY);
  if (fhandle < 0)
  {
    std::cerr << "Error opening GRIB file " << argv[1] << std::endl;
    return -1;
  }

  unsigned char *data = 0;
  data = (unsigned char *) mmap(0, size, PROT_READ, MAP_PRIVATE, fhandle, 0);
  if (data == MAP_FAILED)
  {
    std::cerr << "Error mmapping GRIB file " << argv[1] << std::endl;
    return -1;
  }

  do
  {
    while (pos < limit)
    {
      char *cpnt = (char *) data;
      if (*(cpnt+pos)   == 'G' &&
          *(cpnt+pos+1) == 'R' &&
          *(cpnt+pos+2) == 'I' &&
          *(cpnt+pos+3) == 'B') break;
      ++pos;
    }
    if (pos >= limit)
      break;

    long grid_id;
    long reft_id;
    long valt_id;
    long var_id;
    long lvl_id;
    long atmo_id;
    long msg_id;

    char reftime[64];

    char tname[512];
    char time1[64];
    char time2[64];
    long tdesc;

    char varname[128];
    char pjdef[256];

    char levelname[128];
    double level1, level2;
    long ldesc;

    earthmap m;

    unsigned char *msg = data+pos;
    int version = msg[7];

    if (verbose)
      std::cout << "Found GRIB Version " << version << " message." << std::endl;

    long reclen = 0;

    if (version == 1)
    {
      // GRIB Version 1 messages always contain single data layer
      // File is structured as PDS+[[GDS]+[BMS]]+BDS
      // Product description is always present in PDS
      // We will not accept layers without grid description GDS
      // We will manage mask information in BMS
      // We extract and catologue data in the binary data section BDS
      t_g1_pinfo pkci;

      reclen = (msg[4]<<16) + (msg[5]<<8) + msg[6];

      // Process GRIB1 message to extract data

      // Find pointers to all sections

      unsigned char *pds = msg+8;
      try
      {
        if ((pds[7] & 128) == 0)
        {
          throw "Unsupported grid definition.";
        }
        size_t pdslen = ((size_t) ((pds[0]<<16)+(pds[1]<<8)+pds[2]));
        unsigned char *gds = pds+pdslen;
        size_t gdslen = ((size_t) ((gds[0]<<16)+(gds[1]<<8)+gds[2]));
        unsigned char *bms = 0;
        unsigned char *bds = 0;
        size_t bmslen = 0;
        unsigned char *bitmap = 0;
        if ((pds[7] & 64) != 0)
        {
          bms = gds+gdslen;
          bmslen = ((size_t) ((bms[0]<<16)+(bms[1]<<8)+bms[2]));
          if (bmslen < 8)
            throw "Bitmap not provided by the center. Unusable message";
          bds = bms+bmslen;
          bitmap = bms+6;
        }
        else
          bds = gds+gdslen;

        // Manage reference time
        grib1_pds_to_reftime(pds, reftime);
        if (check_reftime(ppDb, reftime, &reft_id) != 0)
          throw "Reference time not found.";
        if (verbose)
          std::cout << "Reference time " << reftime << " : " 
                    << reft_id << std::endl;

        // Manage time
        grib1_pds_to_time(ppDb, pds, reftime, tname, &tdesc, time1, time2);
        if (check_time(ppDb, reft_id, tdesc, time1, time2, &valt_id) != 0)
          throw "Valid time not found.";
        if (verbose)
        {
          std::cout << tname << " = ";
          if (strcmp(time1, time2))
            std::cout << time1 << " : " << lvl_id << std::endl;
          else
            std::cout << time1 << " to " << time2
                      << " : " << valt_id << std::endl;
        }

        // Manage level information
        grib1_pds_to_level(ppDb, pds, levelname, &ldesc, &level1, &level2);
        if (check_level(ppDb, ldesc, level1, level2, &lvl_id) != 0)
          throw "Valid level not found.";
        if (verbose)
        {
          std::cout << "Level " << levelname << " = ";
          if (level1 == level2)
            std::cout << level1 << " : " << lvl_id << std::endl;
          else
            std::cout << level1 << " to " << level2
                      << " : " << lvl_id << std::endl;
        }

        // Manage Grid
        grib1_gds_to_emap(gds, pjdef, &m);
        if (check_grid(ppDb, pjdef, &grid_id, &m) != 0)
          throw "Grid not found.";
        if (verbose)
          std::cout << "Grid " << pjdef << " : " << grid_id << std::endl;

        // Manage variable informations
        grib1_pds_to_g1_parameter(ppDb, pds, varname, &var_id);
        if (verbose)
          std::cout << "Variable " << varname << " : " << var_id << std::endl;

        if (check_atmovar(ppDb, version, var_id, lvl_id, &atmo_id) != 0)
          throw "Not mapped to atmovar";

        std::cout << "Processing " << varname << " : " << atmo_id << std::endl;
        
        float temp = int_power(10.0, - INT2(pds[26],pds[27]));
        pkci.ndata = ((int) ((gds[6] << 8) + gds[7])) *
                     ((int) ((gds[8] << 8) + gds[9]));
        pkci.nbits = (int) bds[10];
        pkci.reference = temp*(ibm2flt(bds+6));
        pkci.scale = temp*int_power(2.0, INT2(bds[4],bds[5]));
        pkci.bitmap = bitmap;

        // Manage data extraction, decompression and saving
        grib1_to_atmodata(ppDb, pkci, bds, scode_id, valt_id, grid_id,
                           atmo_id, &msg_id);

      }
      catch (char const* str)
      {
        std::cerr << "Skipped this message : " << str << std::endl;
      }
    }
    else if (version == 2)
    {
      t_g2_pinfo pkci;
      pkci.rlv = 0;
      pkci.mask = 0;

      // GRIB Version 2 messages can contain multiple data layer

      gbit(msg, &reclen, 96, 32);

      // Process GRIB message to extract data

      std::vector <long> seclen;
      std::vector <short> secnum;
      std::vector <unsigned char *> secpnt;

      bool endsec = false;

      // Find pointers to all the sections

      // Section 0 is always 16 bytes
      unsigned char *secp = msg+16;
      do
      {
        if (secp[0] == '7' ||
            secp[1] == '7' ||
            secp[2] == '7' ||
            secp[3] == '7')
        {
          endsec = true;
          break;
        }
        long secl = 0;
        short secn = 0;
        gbit(secp, &secl, 0, 32);
        secn = *(secp+4);
        secpnt.push_back(secp);
        seclen.push_back(secl);
        secnum.push_back(secn);
        secp += secl;
      }

      // Process all sections

      while ((! endsec) && (secp-msg < reclen));
      if (endsec)
      {
        // Process well formed GRIB2 message

        long discipline = *(msg+6);
        for (int i = 0; i < (int) seclen.size(); i ++)
        {
          try
          {
            switch(secnum[i])
            {
              case 1: // I have in the id section reference time info
                {
                  grib2_sec1_to_reftime(secpnt[i], reftime);
                  if (check_reftime(ppDb, reftime, &reft_id) != 0)
                    throw "Reference time not found.";
                  if (verbose)
                    std::cout << "Reference time " << reftime << " : " 
                              << reft_id << std::endl;
                }
                break;

              case 3: // Grid definition: geolocation
                {
                  earthmap m;
                  grib2_sec3_to_emap(secpnt[i], pjdef, &m);
                  if (check_grid(ppDb, pjdef, &grid_id, &m) != 0)
                    throw "Grid not found.";
                 if (verbose)
                    std::cout << "Grid " << pjdef << " : "
                              << grid_id << std::endl;

                  gbit(secpnt[i], &(pkci.ndata), 6*8, 32);
                }
                break;

              case 4: // Product Definition Section
                grib2_sec4_to_time(ppDb, secpnt[i], reftime,
                                   tname, &tdesc, time1, time2);
                if (check_time(ppDb, reft_id, tdesc,
                               time1, time2, &valt_id) != 0)
                  throw "Valid time not found.";
                if (verbose)
                {
                  std::cout << "Time " << tname << " = ";
                  if (strcmp(time1, time2))
                    std::cout << time1 << " : " << valt_id << std::endl;
                  else
                    std::cout << time1 << " to " << time2
                              << " : " << valt_id << std::endl;
                }

                grib2_sec4_to_level(ppDb, secpnt[i], levelname,
                                    &ldesc, &level1, &level2);
                if (check_level(ppDb, ldesc, level1, level2, &lvl_id) != 0)
                  throw "Valid level not found.";
                if (verbose)
                {
                  std::cout << "Level " << levelname << " = ";
                  if (level1 == level2)
                    std::cout << level1 << " : " << lvl_id << std::endl;
                  else
                    std::cout << level1 << " to " << level2
                              << " : " << lvl_id << std::endl;
                }

                grib2_sec4_to_g2_parameter(ppDb, discipline,
                                           secpnt[i], varname, &var_id);
                if (verbose)
                  std::cout << "Variable " << varname << " : "
                            << var_id << std::endl;

                if (check_atmovar(ppDb, version, var_id, lvl_id, &atmo_id) != 0)
                  throw "Not mapped to atmovar";
        
                std::cout << "Processing " << varname << " : " 
                          << atmo_id << std::endl;
                break;

              case 5: // Data Representation
                grib2_sec5_to_pkcinfo(secpnt[i], pkci);
                break;

              case 6: // Bitmap
                pkci.bitmap_flag = secpnt[i][5];
                pkci.mask = secpnt[i] + 6;
                break;

              case 7: // Data
                grib2_to_atmodata(ppDb, pkci, secpnt[i], scode_id, valt_id,
                         grid_id, atmo_id, &msg_id);
                break;

              default: // Nothing useful to process
                break;
            }
          }
          catch (char const* str)
          {
            if (pkci.rlv) delete [] pkci.rlv;
            std::cerr << "Skipped this message : " << str << std::endl;
            break;
          }
        }
        if (pkci.rlv) delete [] pkci.rlv;
      }
      else
      {
        std::cerr << "Skipped malformed message." << std::endl;
      }
    }
    if (verbose)
      std::cout << "Record len " << reclen << " at " << pos << std::endl;
    
    pos += (size_t) reclen;
  } while (pos < size);

  if (data)    munmap(data, size);
  if (fhandle) close(fhandle);

  sprintf(sqlbuff, "BEGIN TRANSACTION;");
  if (sqlite3_exec(ppDb, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    throw "Cannot enter transaction";
  }
  sqlite3_free(sqlerr);
  sprintf(sqlbuff, "INSERT OR IGNORE INTO "
                   "ATMODATA(source_id,plgr_id,vtime_id,var_id,value) "
                   "SELECT * FROM tmpins;");
  if (sqlite3_exec(ppDb, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    sqlite3_exec(ppDb, "ROLLBACK", NULL, NULL, &sqlerr);
    sqlite3_free(sqlerr);
    throw "Cannot enter transaction";
  }
  sqlite3_free(sqlerr);
  sprintf(sqlbuff, "COMMIT;");
  if (sqlite3_exec(ppDb, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    sqlite3_exec(ppDb, "ROLLBACK", NULL, NULL, &sqlerr);
    sqlite3_free(sqlerr);
    throw "Cannot commit transaction";
  }
  sqlite3_free(sqlerr);
  sqlite3_exec(ppDb, "DROP TABLE tmpins", NULL, NULL, &sqlerr);
  sqlite3_free(sqlerr);
  if (sqlite3_shutdown() != SQLITE_OK)
  {
    std::cerr << "Database shutdown error !!!" << std::endl;
    return -1;
  }

  std::cout << "Done." << std::endl;
  return 0;
}

void grib1_gds_to_emap(unsigned char *gds, char *pjdef, earthmap *m)
{
  switch (gds[5])
  {
    case PROJ_G1_LATLON: // Regular Latitude/Longitude
      {
        t_proj_latlon_parameters p;
        p.nx = ((long) ((gds[6] << 8) + gds[7]));
        p.ny = ((long) ((gds[8] << 8) + gds[9]));
        p.lat1 = 0.001 * INT3(gds[10],gds[11],gds[12]);
        p.lon1 = 0.001 * INT3(gds[13],gds[14],gds[15]);
        p.lat2 = 0.001 * INT3(gds[17],gds[18],gds[19]);
        p.lon2 = 0.001 * INT3(gds[20],gds[21],gds[22]);
        p.scanflag = gds[27];
        m->set(p);
        sprintf(pjdef, "regll,%ldx%ld,%8.4f,%8.4f,%8.4f,%8.4f,%d",
          p.nx, p.ny, p.lat1, p.lon1, p.lat2, p.lon2, p.scanflag);
      }
      break;
    case PROJ_G1_MERC: // Mercator
      {
        t_proj_mercator_parameters p;
        p.lat1 = 0.001 * INT3(gds[10],gds[11],gds[12]);
        p.lon1 = 0.001 * INT3(gds[13],gds[14],gds[15]);
        p.lat2 = 0.001 * INT3(gds[17],gds[18],gds[19]);
        p.lon2 = 0.001 * INT3(gds[20],gds[21],gds[22]);
        p.dx = INT3(gds[28],gds[29],gds[30]);
        p.dy = INT3(gds[31],gds[32],gds[33]);
        p.reflat = 0.001 * INT3(gds[23],gds[24],gds[25]);
        p.nx = ((long) UINT2(gds[6],gds[7]));
        p.ny = ((long) UINT2(gds[8],gds[9]));
        p.scanflag = gds[27];
        m->set(p);
        sprintf(pjdef,
          "merc,%ldx%ld,%8.4f,%8.4f,%8.4f,%8.4f,%8.4f,%8.4f,%8.4f,%d",
          p.nx, p.ny, p.lat1, p.lon1, p.lat2, p.lon2, p.dx, p.dy,
          p.reflat, p.scanflag);
      }
      break;
    case PROJ_G1_LC: // Lambert Conformal
      {
        t_proj_lambert_parameters p;
        p.nx = ((long) ((gds[6] << 8) + gds[7]));
        p.ny = ((long) ((gds[8] << 8) + gds[9]));
        p.lat1 = 0.001 * INT3(gds[10],gds[11],gds[12]);
        p.lon1 = 0.001 * INT3(gds[13],gds[14],gds[15]);
        p.orient = 0.001 * INT3(gds[17],gds[18],gds[19]);
        p.dx = INT3(gds[20],gds[21],gds[22]);
        p.dy = INT3(gds[23],gds[24],gds[25]);
        p.truelat1 = 0.001 * INT3(gds[28],gds[29],gds[30]);
        p.truelat2 = 0.001 * INT3(gds[31],gds[32],gds[33]);
        p.resflag = gds[16];
        p.projflag = gds[26];
        p.scanflag = gds[27];
        m->set(p);
        sprintf(pjdef,
          "lcc,%ldx%ld,%8.4f,%8.4f,%8.4f,%8.4f,%8.4f,%8.4f,%8.4f,%d,%d,%d",
          p.nx, p.ny, p.lat1, p.lon1, p.orient, p.dx, p.dy, p.truelat1,
          p.truelat2, p.resflag, p.projflag, p.scanflag);
      }
      break;
    case PROJ_G1_PS: // Polar Stereographic
      {
	t_proj_polar_parameters p;
        p.nx = ((long) ((gds[6] << 8) + gds[7]));
        p.ny = ((long) ((gds[8] << 8) + gds[9]));
        p.lat1 = 0.001 * INT3(gds[10],gds[11],gds[12]);
        p.lon1 = 0.001 * INT3(gds[13],gds[14],gds[15]);
        p.orient = 0.001 * INT3(gds[17],gds[18],gds[19]);
        p.dx = INT3(gds[20],gds[21],gds[22]);
        p.dy = INT3(gds[23],gds[24],gds[25]);
        p.resflag = gds[16];
        p.projflag = gds[26];
        p.scanflag = gds[27];
        m->set(p);
        sprintf(pjdef,
          "ps,%ldx%ld,%8.4f,%8.4f,%8.4f,%8.4f,%8.4f,%d,%d,%d",
          p.nx, p.ny, p.lat1, p.lon1, p.orient, p.dx, p.dy,
          p.resflag, p.projflag, p.scanflag);
      }
      break;
    case PROJ_G1_ROTLAT: // Rotated Latitude/Longitude
      {
        t_proj_rotlat_parameters p;
        p.nx = ((long) ((gds[6] << 8) + gds[7]));
        p.ny = ((long) ((gds[8] << 8) + gds[9]));
        p.lat1 = 0.001 * INT3(gds[10],gds[11],gds[12]);
        p.lon1 = 0.001 * INT3(gds[13],gds[14],gds[15]);
        p.lat2 = 0.001 * INT3(gds[17],gds[18],gds[19]);
        p.lon2 = 0.001 * INT3(gds[20],gds[21],gds[22]);
        p.polelat = INT3(gds[32],gds[33],gds[34]);
        p.polelon = INT3(gds[35],gds[36],gds[37]);
        p.rotation = ibm2flt(gds+38);
        p.scanflag = gds[27];
        m->set(p);
        sprintf(pjdef,
          "rll,%ldx%ld,%8.4f,%8.4f,%8.4f,%8.4f,%8.4f,%8.4f,%8.4f,%d",
          p.nx, p.ny, p.lat1, p.lon1, p.lat2, p.lon2, p.polelat,
          p.polelon,p.rotation,p.scanflag);
      }
      break;
    case PROJ_G1_ROTWRF: // Arakawa E 2D on Rotated Latitude/Longitude
      {
        t_proj_rotwrf_parameters p;
        p.nx = ((long) ((gds[6] << 8) + gds[7]));
        p.ny = ((long) ((gds[8] << 8) + gds[9]));
        p.lat1 = 0.001 * INT3(gds[10],gds[11],gds[12]);
        p.lon1 = 0.001 * INT3(gds[13],gds[14],gds[15]);
        p.lat0 = 0.001 * INT3(gds[17],gds[18],gds[19]);
        p.lon0 = 0.001 * INT3(gds[20],gds[21],gds[22]);
        p.resflag = gds[16];
        p.scanflag = gds[27];
        m->set(p);
        sprintf(pjdef,
          "ss2dE,%ldx%ld,%8.4f,%8.4f,%8.4f,%8.4f,%d,%d",
          p.nx, p.ny, p.lat1, p.lon1, p.lat0, p.lon0,
          p.resflag, p.scanflag);
      }
      break;
    default:
      throw "Unsupported grid projection.";
  }
  return;
}

void grib1_pds_to_reftime(unsigned char *pds, char *reftime)
{
  long year, month, day, hour, minute, second;
  year = (pds[12] + 100*(pds[24] - 1));
  month = pds[13];
  day = pds[14];
  hour = pds[15];
  minute = pds[16];
  second = 0;
  sprintf(reftime, "%04ld%02ld%02ldT%02ld%02ld%02ld",
          year, month, day, hour, minute, second);
  return;
}

void grib1_pds_to_g1_parameter(sqlite3 *db, unsigned char *pds,
                               char *vname, long *v_id)
{
  int nr, nc;
  char **result = 0;

  long centid, subcid, procid, tabid, varid;

  tabid = pds[3];
  centid = pds[4];
  procid = pds[5];
  subcid = pds[25];
  varid = pds[8];

  if (centid == 7 && tabid <= 3)
  {
    if (subcid == 1) {subcid = -1; tabid = 32;}
    if (subcid != 0 || (procid != 80 && procid!= 180) ||
        (tabid != 1 && tabid != 2)) {centid = -1; subcid = -1; tabid = -1;}
  }

  // Old Himet GRIB "fake" codes
  if (centid == 100 && subcid == 2 && tabid == 128)
  {
    centid = 81;  // Rome
    subcid = 128; //
    tabid = 2; // International exchange
  }

  sprintf(sqlbuff, "%s%ld%s%ld%s%ld%s%ld%s",
     "SELECT id,name,descr,units FROM g1_parameter WHERE center_id=", centid,
     " AND (subcen_id=", subcid, " OR subcen_id=-1) AND (table_id=", tabid,
     " OR table_id=-1) AND parameter_id=", varid,
     " ORDER BY subcen_id desc, table_id desc;");
  if (sqlite3_get_table(db, sqlbuff, &result, &nr, &nc, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    throw sqlerr;
  }
  sqlite3_free(sqlerr);
  if (nr != 1)
  {
    sqlite3_free_table(result);
    // Try NCEP Default
    sprintf(sqlbuff, "%s%s",
          "SELECT id,name,descr,units FROM g1_parameter WHERE center_id=-1",
          " AND subcen_id=-1 AND table_id=-1 AND parameter_id=-1;");
    if (sqlite3_get_table(db, sqlbuff, &result, &nr, &nc, &sqlerr) != SQLITE_OK)
    {
      std::cerr << "Cannot exec query : " << sqlerr << std::endl;
      sqlite3_free(sqlerr);
      throw sqlerr;
    }
    sqlite3_free(sqlerr);
    if (nr != 1)
    {
      sqlite3_free_table(result);
      sprintf(vname, "Var%03ld : Undefined [undefined]", varid);
      *v_id = -1;
      return;
    }
  }
  sscanf(result[4], "%ld", v_id);
  sprintf(vname, "%s : %s [%s]", result[5], result[6], result[7]);
  sqlite3_free_table(result);
  return;
}

void grib1_pds_to_time(sqlite3 *db, unsigned char *pds, char *reftime,
                       char *tdef, long *tid, char *time1, char *time2)
{
  std::string tmp = reftime;
  ptime r(from_iso_string(tmp));
  
  long tcode = pds[20];
  int t1, t2;

  switch (tcode)
  {
    case 0:
    case 1:
      t1 = t2 = pds[18];
      break;
    case 2: // From-to
    case 3: // Average from-to
    case 4: // Accumulated from-to
    case 5: // Difference from-to
      t1 = pds[18];
      t2 = pds[19];
      break;
    case 10:
      t1 = t2 = pds[18] * 256 + pds[19];
      tcode = 0;
      break;
    default:
      throw "Not a recognized time range by this program";
      break;
  }

  time_duration td1, td2;
  switch(pds[17])
  {
    case 0:
      td1 = minutes(t1);
      td2 = minutes(t2);
      break;
    case 1:
      td1 = hours(t1);
      td2 = hours(t2);
      break;
    case 2:
      td1 = hours(t1*24);
      td1 = hours(t2*24);
      break;
    case 254:
      td1 = seconds(t1);
      td2 = seconds(t2);
      break;
    default:
      throw "Not a recognized time unit by this program";
  }
  ptime vt1 = r + td1;
  ptime vt2 = r + td2;
  sprintf(time1, "%s", to_iso_string(vt1).c_str());
  sprintf(time2, "%s", to_iso_string(vt2).c_str());

  int nr, nc;
  char **result = 0;
  sprintf(sqlbuff, "%s%ld%s",
     "SELECT id,descr FROM timedef WHERE grib1key=", tcode, ";");
  if (sqlite3_get_table(db, sqlbuff, &result, &nr, &nc, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    throw sqlerr;
  }
  sqlite3_free(sqlerr);
  if (nr != 1)
  {
    sqlite3_free_table(result);
    sprintf(tdef, "Time Code %03ld", tcode);
    *tid = -1;
    return;
  }
  sscanf(result[2], "%ld", tid);
  sprintf(tdef, "%s", result[3]);
  sqlite3_free_table(result);
  return;
}

void grib1_pds_to_level(sqlite3 *db, unsigned char *pds,
                        char *levelname, long *lev_id,
                        double *level1, double *level2)
{
  long levcode;
  long value, o1, o2;
  int nr, nc;
  char **result = 0;

  levcode = pds[9];
  value = ((long) ((pds[10]<<8) + pds[11]));
  o1 = pds[10];
  o2 = pds[11];

  switch (levcode)
  {
    case 100:
      *level1 = 100*value;
      *level2 = 0.0;
      break;
    case 101:
      *level1 = o1*1000;
      *level2 = o2*1000;
      break;
    case 104:
      *level1 = o1*100;
      *level2 = o2*100;
      break;
    case 106:
      *level1 = o1*100;
      *level2 = o2*100;
      break;
    case 108:
      *level1 = o1*100;
      *level2 = o2*100;
      break;
    case 110:
      *level1 = o1;
      *level2 = o2;
      break;
    case 111:
      *level1 = ((double) value)/100.0;
      *level2 = 0.0;
      break;
    case 112:
      *level1 = ((double) o1)/100.0;
      *level2 = ((double) o2)/100.0;
      break;
    case 114:
      *level1 = 475.0-o1;
      *level2 = 475.0-o2;
      break;
    case 115:
      *level1 = 100*value;
      *level2 = 0.0;
      break;
    case 116:
      *level1 = o1*100;
      *level2 = o2*100;
      break;
    case 120:
      *level1 = ((double) o1)/100.0;
      *level2 = ((double) o2)/100.0;
      break;
    case 121:
      *level1 = (1100.0-o1)*100;
      *level2 = (1100.0-o2)*100;
      break;
    case 128:
      *level1 = 1.1-(((double)o1)/1000.0);
      *level2 = 1.1-(((double)o2)/1000.0);
      break;
    case 141:
      *level1 = o1*100;
      *level2 = (1100.0-o2)*100;
      break;
    case 235:
      *level1 = *level2 = ((double) value)/10.0;
      break;
    case 236:
      *level1 = o1*10.0;
      *level2 = o2*10.0;
      break;
    default:
      *level1 = value;
      *level2 = 0.0;
      break;
  }
  sprintf(sqlbuff, "%s%ld%s",
     "SELECT id,name,descr,units FROM leveldef WHERE grib1key=", levcode, ";");
  if (sqlite3_get_table(db, sqlbuff, &result, &nr, &nc, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    throw sqlerr;
  }
  sqlite3_free(sqlerr);
  if (nr != 1)
  {
    sqlite3_free_table(result);
    sprintf(levelname, "LEV%03ld : Undefined [undefined]", levcode);
    *lev_id = -1;
    return;
  }
  sscanf(result[4], "%ld", lev_id);
  sprintf(levelname, "%s : %s [%s]", result[5], result[6], result[7]);
  sqlite3_free_table(result);
  return;
}

void grib2_sec4_to_level(sqlite3 *db, unsigned char *secpnt,
                         char *levelname, long *lev_id,
                         double *level1, double *level2)
{
  long pdtn;
  gbit(secpnt, &pdtn, 7*8, 16);
  if (pdtn != 0 && pdtn != 8)
    throw "Not a recognized by this program product definition template";
  
  long levcode;
  int nr, nc;
  char **result = 0;

  levcode = secpnt[22];
  *level1 = scaled2flt(secpnt[23], int4(secpnt+24));
  if (secpnt[29] == 255) *level2 = *level1;
  else *level2 = scaled2flt(secpnt[29], int4(secpnt+30));
  sprintf(sqlbuff, "%s%ld%s",
     "SELECT id,name,descr,units FROM leveldef WHERE grib2key=", levcode, ";");
  if (sqlite3_get_table(db, sqlbuff, &result, &nr, &nc, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    throw sqlerr;
  }
  sqlite3_free(sqlerr);
  if (nr != 1)
  {
    sqlite3_free_table(result);
    sprintf(levelname, "LEV%03ld : Undefined [undefined]", levcode);
    *lev_id = -1;
    return;
  }
  sscanf(result[4], "%ld", lev_id);
  sprintf(levelname, "%s : %s [%s]", result[5], result[6], result[7]);
  sqlite3_free_table(result);
  return;
}

void grib2_sec1_to_reftime(unsigned char *secpnt, char *reftime)
{
  long year, month, day, hour, minute, second;
  gbit(secpnt, &year, 12*8, 16);
  month = secpnt[14];
  day = secpnt[15];
  hour = secpnt[16];
  minute = secpnt[17];
  second = secpnt[18];
  sprintf(reftime, "%04ld%02ld%02ldT%02ld%02ld%02ld",
          year, month, day, hour, minute, second);
  return;
}

void grib2_sec4_to_time(sqlite3 *db, unsigned char *secpnt, char *reftime,
                        char *tdef, long *tid, char *time1, char *time2)
{
  long pdtn;
  gbit(secpnt, &pdtn, 7*8, 16);
  if (pdtn != 0 && pdtn != 8)
    throw "Not a recognized by this program product definition template";
  
  std::string tmp = reftime;
  ptime r(from_iso_string(tmp));
  
  long fcst;
  gbit(secpnt, &fcst, 18*8, 32);
  time_duration td;

  switch(secpnt[17])
  {
    case 0:
      td = minutes(fcst);
      break;
    case 1:
      td = hours(fcst);
      break;
    case 2:
      td = hours(fcst*24);
      break;
    case 13:
      td = seconds(fcst);
      break;
    default:
      throw "Not a recognized time unit by this program";
  }

  ptime vt = r + td;
  sprintf(time1, "%s", to_iso_string(vt).c_str());
  if (pdtn == 0)
  {
    *tid = 1;
    sprintf(time2, "%s", to_iso_string(vt).c_str());
    sprintf(tdef, "%s%s", "Analysis or forecast at a horizontal level or",
            " in a horizontal layer at a point in time.");
  }
  else
  {
    long year, month, day, hour, minute, second;
    gbit(secpnt, &year, 34*8, 16);
    month = secpnt[36];
    day = secpnt[37];
    hour = secpnt[38];
    minute = secpnt[39];
    second = secpnt[40];
    sprintf(time2, "%04ld%02ld%02ldT%02ld%02ld%02ld",
            year, month, day, hour, minute, second);
    int nr, nc;
    char **result = 0;
    long statp = secpnt[46];
    sprintf(sqlbuff, "%s%ld%s",
       "SELECT id,descr FROM timedef WHERE grib2key=", statp, ";");
    if (sqlite3_get_table(db, sqlbuff, &result, &nr, &nc, &sqlerr) != SQLITE_OK)
    {
      std::cerr << "Cannot exec query : " << sqlerr << std::endl;
      sqlite3_free(sqlerr);
      throw sqlerr;
    }
    sqlite3_free(sqlerr);
    if (nr != 1)
    {
      sqlite3_free_table(result);
      sprintf(tdef, "Code %03ld.%03ld", pdtn, statp);
      *tid = -1;
      return;
    }
    sscanf(result[2], "%ld", tid);
    sprintf(tdef, "%s", result[3]);
    sqlite3_free_table(result);
  }

  return;
}

void grib2_sec3_to_emap(unsigned char *secpnt, char *pjdef, earthmap *m)
{
  if (secpnt[5] != 0)
    throw "Unsupported grid definition.";

  long gdtn;
  long basic_ang = int4(secpnt+38);
  long sub_ang = int4(secpnt+42);
  double units = basic_ang == 0 ? 0.000001 :
                 (double) basic_ang / (double) sub_ang;
  gbit(secpnt, &gdtn, 12*8, 16);
  switch (gdtn)
  {
    case PROJ_G2_LATLON: // Regular Latitude/Longitude
      {
        t_proj_latlon_parameters p;
        p.nx = uint4(secpnt+30);
        p.ny = uint4(secpnt+34);
        p.lat1 = units * int4(secpnt+46);
        p.lon1 = units * uint4(secpnt+50);
        p.lat2 = units * int4(secpnt+55);
        p.lon2 = units * uint4(secpnt+59);
        p.scanflag = secpnt[71];
        m->set(p);
        sprintf(pjdef, "regll,%ldx%ld,%8.4f,%8.4f,%8.4f,%8.4f,%d",
          p.nx, p.ny, p.lat1, p.lon1, p.lat2, p.lon2, p.scanflag);
      }
      break;
    case PROJ_G2_ROTLAT: // Rotated Latitude/Longitude
      {
        t_proj_rotlat_parameters p;
        p.nx = uint4(secpnt+30);
        p.ny = uint4(secpnt+34);
        p.lat1 = units * int4(secpnt+46);
        p.lon1 = units * uint4(secpnt+50);
        p.lat2 = units * int4(secpnt+55);
        p.lon2 = units * uint4(secpnt+59);
        p.polelat = units * int4(secpnt+72);
        p.polelon = units * int4(secpnt+76),
        p.rotation = units * int4(secpnt+80);
        p.scanflag = secpnt[71];
        m->set(p);
        sprintf(pjdef,
          "rll,%ldx%ld,%8.4f,%8.4f,%8.4f,%8.4f,%8.4f,%8.4f,%8.4f,%d",
          p.nx, p.ny, p.lat1, p.lon1, p.lat2, p.lon2, p.polelat,
          p.polelon,p.rotation,p.scanflag);
      }
      break;
    case PROJ_G2_MERC: // Mercator
      {
        t_proj_mercator_parameters p;
        p.nx = uint4(secpnt+30);
        p.ny = uint4(secpnt+34);
        p.lat1 = (int4(secpnt+38)*0.000001);
        p.lon1 = (uint4(secpnt+42)*0.000001);
        p.lat2 = (int4(secpnt+51)*0.000001);
        p.lon2 = (uint4(secpnt+55)*0.000001);
        p.dx = ((uint4(secpnt+64))*0.001);
        p.dy = ((uint4(secpnt+68))*0.001);
        p.reflat = (int4(secpnt+47)*0.000001);
        p.scanflag = secpnt[71];
        m->set(p);
        sprintf(pjdef,
          "merc,%ldx%ld,%8.4f,%8.4f,%8.4f,%8.4f,%8.4f,%8.4f,%8.4f,%d",
          p.nx, p.ny, p.lat1, p.lon1, p.lat2, p.lon2, p.dx, p.dy,
          p.reflat, p.scanflag);
      }
      break;
    case PROJ_G2_PS: // Polar Stereographic
      {
        t_proj_polar_parameters p;
        p.nx = uint4(secpnt+30);
        p.ny = uint4(secpnt+34);
        p.lat1 = (int4(secpnt+38)*0.000001);
        p.lon1 = (uint4(secpnt+42)*0.000001);
        p.orient = (uint4(secpnt+51)*0.000001);
        p.dx = (uint4(secpnt+55)*0.001);
        p.dy = (uint4(secpnt+59)*0.001);
        p.resflag = secpnt[54];
        p.projflag = secpnt[63];
        p.scanflag = secpnt[71];
        m->set(p);
        sprintf(pjdef,
          "ps,%ldx%ld,%8.4f,%8.4f,%8.4f,%8.4f,%8.4f,%d,%d,%d",
          p.nx, p.ny, p.lat1, p.lon1, p.orient, p.dx, p.dy,
          p.resflag, p.projflag, p.scanflag);
      }
      break;
    case PROJ_G2_LC: // Lambert Conformal
      {
        t_proj_lambert_parameters p;
        p.nx = uint4(secpnt+30);
        p.ny = uint4(secpnt+34);
        p.lat1 = (int4(secpnt+38) * 0.000001);
        p.lon1 = (int4(secpnt+42) * 0.000001);
        p.orient = (int4(secpnt+51) * 0.000001);
        p.dx = (int4(secpnt+55) * 0.001);
        p.dy = (int4(secpnt+59) * 0.001);
        p.truelat1 = (int4(secpnt+65) * 0.000001);
        p.truelat2 = (int4(secpnt+69) * 0.000001);
        p.resflag = secpnt[54];
        p.projflag = secpnt[63];
        p.scanflag = secpnt[71];
        m->set(p);
        sprintf(pjdef,
          "lcc,%ldx%ld,%8.4f,%8.4f,%8.4f,%8.4f,%8.4f,%8.4f,%8.4f,%d,%d,%d",
          p.nx, p.ny, p.lat1, p.lon1, p.orient, p.dx, p.dy, p.truelat1,
          p.truelat2, p.resflag, p.projflag, p.scanflag);
      }
      break;
    case PROJ_G2_ROTWRF: // Arakawa staggered E-grid
      {
        t_proj_rotwrf_parameters p;
        p.nx = uint4(secpnt+30);
        p.ny = uint4(secpnt+34);
        p.lat1 = units * int4(secpnt+46);
        p.lon1 = units * uint4(secpnt+50);
        p.lat0 = units * int4(secpnt+55);
        p.lon0 = units * uint4(secpnt+59);
        p.resflag = secpnt[54];
        p.scanflag = secpnt[71];
        m->set(p);
        sprintf(pjdef,
          "ss2dE,%ldx%ld,%8.4f,%8.4f,%8.4f,%8.4f,%d,%d",
          p.nx, p.ny, p.lat1, p.lon1, p.lat0, p.lon0,
          p.resflag, p.scanflag);
      }
      break;
    default:
      throw "Unsupported grid projection.";
  }
  return;
}

void grib2_sec4_to_g2_parameter(sqlite3 *db, long disc, unsigned char *secpnt,
                                char *vname, long *v_id)
{
  int nr, nc;
  char **result = 0;

  long category_id, parameter_id;

  category_id = secpnt[9];
  parameter_id = secpnt[10];

  sprintf(sqlbuff, "%s%ld%s%ld%s%ld%s",
     "SELECT id,name,descr,units FROM g2_parameter WHERE discipline_id=", disc,
     " AND category_id=", category_id, " AND parameter_id=", parameter_id, ";");
  if (sqlite3_get_table(db, sqlbuff, &result, &nr, &nc, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    throw sqlerr;
  }
  sqlite3_free(sqlerr);
  if (nr != 1)
  {
    sqlite3_free_table(result);
    sprintf(vname, "Var%03ld.%03ld.%03ld : Undefined [undefined]",
            disc, category_id, parameter_id);
    *v_id = -1;
    return;
  }
  sscanf(result[4], "%ld", v_id);
  sprintf(vname, "%s : %s [%s]", result[5], result[6], result[7]);
  sqlite3_free_table(result);
  return;
}

void grib2_sec5_to_pkcinfo(unsigned char *p, t_g2_pinfo &pkci)
{
  if (pkci.rlv) delete [] pkci.rlv;
  pkci.rlv = 0;
  pkci.packing = uint2(p+9);
  if (pkci.packing == 4)
    pkci.precision = p[11];
  else if (pkci.packing == 0)
  {
    pkci.reference = ieee2flt(p+11);
    pkci.bin_scale = int_power(2.0, int2(p+15));
    pkci.dec_scale = int_power(10.0, -int2(p+17));
    pkci.nbits = p[19];
  }
  else if (pkci.packing == 2 || pkci.packing == 3)
  {
    pkci.reference = ieee2flt(p+11);
    pkci.bin_scale = int_power(2.0, int2(p+15));
    pkci.dec_scale = int_power(10.0, -int2(p+17));
    pkci.nbits = p[19];
    pkci.npnts = uint4(p+5);
    pkci.ngroups = uint4(p+31);
    pkci.ctable_5_4 = p[21];
    pkci.ctable_5_5 = p[22];
    if (pkci.packing == 3)
      pkci.ctable_5_6 = (int) (p[47]);
    else
      pkci.ctable_5_6 = -1;
    if (pkci.packing == 2)
      pkci.extra_octets = p[48];
    else
      pkci.extra_octets = 0;
    pkci.ref_group_width = p[35];
    pkci.nbit_group_width = p[36];
    pkci.ref_group_length = uint4(p+37);
    pkci.group_length_factor = p[41];
    pkci.len_last = uint4(p+42);
    pkci.nbits_group_len = p[46];
    pkci.n_sub_missing = sub_missing_values(p,
                           &(pkci.missing1), &(pkci.missing2));
  }
  else if (pkci.packing == 200)
  {
    pkci.nbits = p[11];
    pkci.mv = uint2(p+12);
    pkci.mvl = uint2(p+14);
    pkci.dec_scale = p[16];
    if (pkci.dec_scale > 127) // convert signed negative values
      pkci.dec_scale = -(pkci.dec_scale - 128);
    pkci.rlv = new double[pkci.mvl];
    if (pkci.rlv == 0) throw "Error allocating memory in grib2_sec5_to_pkcinfo";
    for (int i = 0; i < (int) pkci.mvl; i++)
      pkci.rlv[i] = int2(p+17+i*2);
    double dec_factor = int_power(10.0, -pkci.dec_scale);
    for (int i = 0; i < (int) pkci.mvl; i++)
     pkci.rlv[i] = pkci.rlv[i]*dec_factor; 
  }
  else if (pkci.packing == 40 || pkci.packing == 40000)
  {
    pkci.reference = ieee2flt(p+11);
    pkci.bin_scale = int_power(2.0, int2(p+15));
    pkci.dec_scale = int_power(10.0, -int2(p+17));
    pkci.nbits = p[19];
  }
  else if (pkci.packing == 41)
  {
    pkci.reference = ieee2flt(p+11);
    pkci.bin_scale = int_power(2.0, int2(p+15));
    pkci.dec_scale = int_power(10.0, -int2(p+17));
    pkci.nbits = p[19];
  }
  else
    throw "Unsupported packing found";
  return;
}

int check_reftime(sqlite3 *db, char *reftime, long *reft_id)
{
  int nr, nc;
  char **result = 0;

  char dbdatetime[32];
  int y,m,d,H,M,S;
  sscanf(reftime, "%04d%02d%02dT%02d%02d%02d", &y, &m, &d, &H, &M, &S);
  sprintf(dbdatetime, "%04d-%02d-%02d %02d:%02d:%02d", y,m,d,H,M,S);

  sprintf(sqlbuff, "%s%s%s",
          "SELECT id FROM reftime WHERE reference='", dbdatetime, "';");
  if (sqlite3_get_table(db, sqlbuff, &result, &nr, &nc, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    return -1;
  }
  sqlite3_free(sqlerr);
  if (nr != 1)
  {
    sqlite3_free_table(result);
    sprintf(sqlbuff, "%s", "BEGIN TRANSACTION;");
    if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
    {
      std::cerr << "Cannot exec query : " << sqlerr << std::endl;
      sqlite3_free(sqlerr);
      return -1;
    }
    sqlite3_free(sqlerr);
    sprintf(sqlbuff, "%s%s%s",
            "INSERT OR IGNORE INTO reftime(reference) values ('",
            dbdatetime, "');");
    if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
    {
      std::cerr << "Cannot exec query : " << sqlerr << std::endl;
      sqlite3_free(sqlerr);
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
      return -1;
    }
    sqlite3_free(sqlerr);
    sprintf(sqlbuff, "%s%s%s",
            "SELECT id FROM reftime WHERE reference='", dbdatetime, "';");
    if (sqlite3_get_table(db, sqlbuff, &result, &nr, &nc, &sqlerr) != SQLITE_OK)
    {
      std::cerr << "Cannot exec query : " << sqlerr << std::endl;
      sqlite3_free(sqlerr);
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
      return -1;
    }
    sqlite3_free(sqlerr);
    sprintf(sqlbuff, "%s", "COMMIT;");
    if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
    {
      std::cerr << "Cannot exec query : " << sqlerr << std::endl;
      sqlite3_free_table(result);
      sqlite3_free(sqlerr);
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
      return -1;
    }
    sqlite3_free(sqlerr);
    if (nr != 1)
    {
      sqlite3_free_table(result);
      std::cerr << "Reference time NOT stored in DB : " << reftime
                << std::endl;
      return -1;
    }
  }
  if (sscanf(result[1], "%ld", reft_id) != 1)
  {
    std::cerr << "Not parsable reference time stored in DB : " << result[1]
              << std::endl;
    sqlite3_free_table(result);
    return -1;
  }
  sqlite3_free_table(result);
  return 0;
}

int check_level(sqlite3 *db, long ldes, double l1, double l2, long *lev_id)
{
  int nr, nc;
  char **result = 0;

  sprintf(sqlbuff, "%s%ld%s%f%s%f%s",
          "SELECT id FROM levels WHERE leveldef_id=", ldes,
          " AND level1=", l1, " AND level2=", l2, ";");
  if (sqlite3_get_table(db, sqlbuff, &result, &nr, &nc, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    return -1;
  }
  sqlite3_free(sqlerr);
  if (nr != 1)
  {
    sqlite3_free_table(result);
    sprintf(sqlbuff, "%s", "BEGIN TRANSACTION;");
    if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
    {
      std::cerr << "Cannot exec query : " << sqlerr << std::endl;
      sqlite3_free(sqlerr);
      return -1;
    }
    sqlite3_free(sqlerr);
    sprintf(sqlbuff, "%s%s%ld%s%f%s%f%s",
            "INSERT OR IGNORE INTO",
            " levels(leveldef_id, level1, level2) values (",
            ldes, ",", l1, ",", l2, ");");
    if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
    {
      std::cerr << "Cannot exec query : " << sqlerr << std::endl;
      sqlite3_free(sqlerr);
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
      return -1;
    }
    sqlite3_free(sqlerr);
    sprintf(sqlbuff, "%s%ld%s%f%s%f%s",
            "SELECT id FROM levels WHERE leveldef_id=", ldes,
            " AND level1=", l1, " AND level2=", l2, ";");
    if (sqlite3_get_table(db, sqlbuff, &result, &nr, &nc, &sqlerr) != SQLITE_OK)
    {
      std::cerr << "Cannot exec query : " << sqlerr << std::endl;
      sqlite3_free(sqlerr);
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
      return -1;
    }
    sqlite3_free(sqlerr);
    sprintf(sqlbuff, "%s", "COMMIT;");
    if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
    {
      std::cerr << "Cannot exec query : " << sqlerr << std::endl;
      sqlite3_free_table(result);
      sqlite3_free(sqlerr);
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
      return -1;
    }
    sqlite3_free(sqlerr);
    if (nr != 1)
    {
      sqlite3_free_table(result);
      std::cerr << "Level NOT stored in DB : " << l1 << " : "
                << l2 << std::endl;
      return -1;
    }
  }
  if (sscanf(result[1], "%ld", lev_id) != 1)
  {
    std::cerr << "Not parsable valid time stored in DB : " << result[1]
              << std::endl;
    sqlite3_free_table(result);
    return -1;
  }
  sqlite3_free_table(result);
  return 0;
}

int check_time(sqlite3 *db, long reft_id, long tid,
               char *time1, char *time2, long *valt_id)
{
  int nr, nc;
  char **result = 0;

  char dbdatetimev1[32];
  char dbdatetimev2[32];
  int y,m,d,H,M,S;
  sscanf(time1, "%04d%02d%02dT%02d%02d%02d", &y, &m, &d, &H, &M, &S);
  sprintf(dbdatetimev1, "%04d-%02d-%02d %02d:%02d:%02d", y,m,d,H,M,S);
  sscanf(time2, "%04d%02d%02dT%02d%02d%02d", &y, &m, &d, &H, &M, &S);
  sprintf(dbdatetimev2, "%04d-%02d-%02d %02d:%02d:%02d", y,m,d,H,M,S);

  sprintf(sqlbuff, "%s%ld%s%ld%s%s%s%s%s",
          "SELECT id FROM valtime WHERE ref_id=", reft_id,
          " AND timedef_id=", tid, " AND valid1='", dbdatetimev1,
          "' AND valid2='", dbdatetimev2, "';");
  if (sqlite3_get_table(db, sqlbuff, &result, &nr, &nc, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    return -1;
  }
  sqlite3_free(sqlerr);
  if (nr != 1)
  {
    sqlite3_free_table(result);
    sprintf(sqlbuff, "%s", "BEGIN TRANSACTION;");
    if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
    {
      std::cerr << "Cannot exec query : " << sqlerr << std::endl;
      sqlite3_free(sqlerr);
      return -1;
    }
    sqlite3_free(sqlerr);
    sprintf(sqlbuff, "%s%s%ld%s%ld%s%s%s%s%s",
            "INSERT OR IGNORE INTO",
            " valtime(ref_id, timedef_id, valid1, valid2) values (",
            reft_id, ",", tid, ",'", dbdatetimev1, "','", dbdatetimev2, "');");
    if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
    {
      std::cerr << "Cannot exec query : " << sqlbuff
                << " : Reason is : " << sqlerr << std::endl;
      sqlite3_free(sqlerr);
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
      return -1;
    }
    sqlite3_free(sqlerr);
    sprintf(sqlbuff, "%s%ld%s%ld%s%s%s%s%s",
            "SELECT id FROM valtime WHERE ref_id=", reft_id,
            " AND timedef_id=", tid, " AND valid1='", dbdatetimev1,
            "' AND valid2='", dbdatetimev2, "';");
    if (sqlite3_get_table(db, sqlbuff, &result, &nr, &nc, &sqlerr) != SQLITE_OK)
    {
      std::cerr << "Cannot exec query : " << sqlerr << std::endl;
      sqlite3_free(sqlerr);
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
      return -1;
    }
    sqlite3_free(sqlerr);
    sprintf(sqlbuff, "%s", "COMMIT;");
    if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
    {
      std::cerr << "Cannot exec query : " << sqlerr << std::endl;
      sqlite3_free_table(result);
      sqlite3_free(sqlerr);
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
      return -1;
    }
    sqlite3_free(sqlerr);
    if (nr != 1)
    {
      std::cerr << "Valid time NOT stored in DB : " << time1 << " : "
                << time2 << std::endl;
      sqlite3_free_table(result);
      return -1;
    }
  }
  if (sscanf(result[1], "%ld", valt_id) != 1)
  {
    std::cerr << "Not parsable valid time stored in DB : " << result[1]
              << std::endl;
    sqlite3_free_table(result);
    return -1;
  }
  sqlite3_free_table(result);
  return 0;
}

int check_grid(sqlite3 *db, char *pjdef, long *grid, earthmap *m)
{
  int nr, nc;
  char **result = 0;

  sprintf(sqlbuff, "%s%s%s",
          "SELECT id FROM grid WHERE emap='", pjdef, "';");
  if (sqlite3_get_table(db, sqlbuff, &result, &nr, &nc, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    return -1;
  }
  sqlite3_free(sqlerr);
  if (nr != 1)
  {
    sqlite3_free_table(result);
    sprintf(sqlbuff, "%s", "BEGIN TRANSACTION;");
    if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
    {
      std::cerr << "Cannot exec query : " << sqlerr << std::endl;
      sqlite3_free(sqlerr);
      return -1;
    }
    sqlite3_free(sqlerr);
    std::cout << "New Grid found !" << std::endl;
    std::cout << "Inserting grid into DB and calculating indexes." << std::endl;
    // Insert it in the DB
    sprintf(sqlbuff, "%s%s%s",
            "INSERT OR IGNORE INTO grid(emap) values ('",
            pjdef, "');");
    if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
    {
      std::cerr << "Cannot exec query : " << sqlerr << std::endl;
      sqlite3_free(sqlerr);
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
      return -1;
    }
    sqlite3_free(sqlerr);
    sprintf(sqlbuff, "%s%s%s",
            "SELECT id FROM grid WHERE emap='", pjdef, "';");
    if (sqlite3_get_table(db, sqlbuff, &result, &nr, &nc, &sqlerr) != SQLITE_OK)
    {
      std::cerr << "Cannot exec query : " << sqlerr << std::endl;
      sqlite3_free(sqlerr);
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
      return -1;
    }
    sqlite3_free(sqlerr);
    if (nr != 1)
    {
      std::cerr << "DB Error : Cannot find inserted grid !!!" << std::endl;
      sqlite3_free_table(result);
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
      return -1;
    }
    if (sscanf(result[1], "%ld", grid) != 1)
    {
      std::cerr << "Not parsable grid stored in DB : " << result[1]
                << std::endl;
      sqlite3_free_table(result);
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
      return -1;
    }
    sqlite3_free_table(result);

    // For all places, fill a record in place/grid

    sprintf(sqlbuff, "%s", "SELECT id,latitude,longitude FROM place;");
    if (sqlite3_get_table(db, sqlbuff, &result, &nr, &nc, &sqlerr) != SQLITE_OK)
    {
      std::cerr << "Cannot exec query : " << sqlerr << std::endl;
      sqlite3_free(sqlerr);
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
      return -1;
    }
    sqlite3_free(sqlerr);

    long place_id;
    double pixel, line, crot, srot, locres;
    double lat, lon, ri, rj, latp, lonp, latm, lonm;
    std::cout << "Place records are " << nr << std::endl;
    for (int i = 3; i <= nr*3; i += 3)
    {
      if (sscanf(result[i], "%ld", &place_id) != 1 ||
          sscanf(result[i+1], "%lf", &lat)    != 1 ||
          sscanf(result[i+2], "%lf", &lon)    != 1)
      {
        std::cerr << "Not parsable place stored in DB : " << result[i]
                  << ", " << result[i+1] << ", " << result[i+2] << std::endl;
        sqlite3_free_table(result);
        sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
        sqlite3_free(sqlerr);
        return -1;
      }
      m->latlon_to_ij(lat, lon, &ri, &rj);
      m->ij_to_latlon(floor(ri), floor(rj), &latm, &lonm);
      m->ij_to_latlon(ceil(ri), ceil(rj), &latp, &lonp);
      locres = sqrt(pow(latp-latm, 2.0)+pow(lonp-lonm, 2.0));
      m->ij2pl(ri, rj, &pixel, &line);
      m->rotuv_ll2gr(lat, lon, &crot, &srot);
      if (pixel > 1.0e20) continue;
      long ipos = m->ij2p((int) rintf(pixel), (int) rintf(line));
      sprintf(sqlbuff, "%s%ld,%ld,%f,%f,%ld,%f,%f,%f%s",
          "INSERT OR IGNORE INTO place_grid("
          "place_id,grid_id,pixel,line,ipos,crot,srot,locres) "
          "values (",
          place_id, *grid, pixel, line, ipos, crot, srot, locres, ");");
      if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
      {
        std::cerr << "Cannot exec query : " << sqlerr << std::endl;
        sqlite3_free(sqlerr);
        sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
        sqlite3_free(sqlerr);
        return -1;
      }
      sqlite3_free(sqlerr);
    }
    sprintf(sqlbuff, "%s", "COMMIT;");
    if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
    {
      std::cerr << "Cannot exec query : " << sqlerr << std::endl;
      sqlite3_free(sqlerr);
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
      return -1;
    }
    sqlite3_free(sqlerr);
    sqlite3_free_table(result);
    return 0;
  }
  if (sscanf(result[1], "%ld", grid) != 1)
  {
    std::cerr << "Not parsable grid stored in DB : " << result[1]
              << std::endl;
    sqlite3_free_table(result);
    return -1;
  }
  sqlite3_free_table(result);
  return 0;
}

int check_atmovar(sqlite3 *db, int gv, long g_id, long lev_id, long *atmo_id)
{
  int nr, nc;
  char **result = 0;

  switch (gv)
  {
    case 1:
      sprintf(sqlbuff, "%s%ld%s%ld%s",
            "SELECT atmovar_id FROM atmovar_g1param WHERE g1_id=", g_id,
            " AND level_id=", lev_id, ";");
      break;
    case 2:
      sprintf(sqlbuff, "%s%ld%s%ld%s",
            "SELECT atmovar_id FROM atmovar_g2param WHERE g2_id=", g_id,
            " AND level_id=", lev_id, ";");
      break;
    default:
      throw "Unknown GRIB version file";
      break;
  }
  if (sqlite3_get_table(db, sqlbuff, &result, &nr, &nc, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    return -1;
  }
  sqlite3_free(sqlerr);
  if (nr != 1)
  {
    sqlite3_free_table(result);
    return -1;
  }
  if (sscanf(result[1], "%ld", atmo_id) != 1)
  {
    std::cerr << "Not parsable atmovar stored in DB : " << result[1]
              << std::endl;
    sqlite3_free_table(result);
    return -1;
  }
  sqlite3_free_table(result);
  return 0;
}

void grib1_to_atmodata(sqlite3 *db,
                      t_g1_pinfo &pkci,
                      unsigned char *bds,
                      long scode_id, long vtid, long gid, long vid, long *mid)
{
  float *df = new float[pkci.ndata];
  if (unpk_grib1(df, bds, pkci.bitmap, pkci.nbits, pkci.ndata,
                 pkci.reference, pkci.scale) != 0)
  {
    delete [] df;
    throw "Error unpacking GRIB1 data";
  }

  if (verbose)
  {
    float maxval = FLT_MIN;
    float minval = FLT_MAX;
    for (int i = 0; i < pkci.ndata; i ++)
    {
      if (maxval < df[i]) maxval = df[i];
      if (minval > df[i]) minval = df[i];
    }
    std::cout << "Max value : " << maxval << std::endl;
    std::cout << "Min value : " << minval << std::endl;
  }

  // search gid in gp
  t_gpoints *gpp = 0;
  for (int i = 0; i < (int) gp.size(); i ++)
    if (gp[i].gid == gid) gpp = &(gp[i]);

  if (gpp == 0)
  {
    int nr, nc;
    char **result = 0;
    t_gpoints gpn;
    gpn.gid = gid;
    sprintf(sqlbuff, "SELECT id,ipos FROM place_grid WHERE grid_id=%ld;", gid);
    if (sqlite3_get_table(db, sqlbuff, &result, &nr, &nc, &sqlerr) != SQLITE_OK)
    {
      std::cerr << "Cannot exec query : " << sqlerr << std::endl;
      sqlite3_free(sqlerr);
      throw sqlerr;
    }
    sqlite3_free(sqlerr);
    if (nr == 0)
    {
      delete [] df;
      sqlite3_free_table(result);
      return;
    }
    long pid;
    long ipos;
    for (int i = 0; i < nr; i ++)
    {
      sscanf(result[2+i*2], "%ld", &pid);
      sscanf(result[2+i*2+1], "%ld", &ipos); 
      gpn.pid.push_back(pid);
      gpn.ipos.push_back(ipos);
    }
    gp.push_back(gpn);
    gpp = &(gp.back());
    sqlite3_free_table(result);
  }

  sprintf(sqlbuff, "INSERT OR IGNORE INTO tmpins VALUES(?,?,?,?,?)");
  sqlite3_stmt *stmt = 0;
  if (sqlite3_prepare_v2(db, sqlbuff, -1, &stmt, NULL) != SQLITE_OK)
  {
    std::cerr << "Cannot prepare query." << std::endl;
    throw sqlite3_errmsg(db);
  }
  sprintf(sqlbuff, "BEGIN TRANSACTION;");
  if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    throw "Cannot enter transaction";
  }
  sqlite3_free(sqlerr);

  double value;
  for (int i = 0; i < (int) gpp->pid.size(); i ++)
  {
    value = df[gpp->ipos[i]];
    if (sqlite3_bind_int(stmt, 1, scode_id) != SQLITE_OK)
    {
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
    }
    if (sqlite3_bind_int(stmt, 2, gpp->pid[i]) != SQLITE_OK)
    {
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
    }
    if (sqlite3_bind_int(stmt, 3, vtid) != SQLITE_OK)
    {
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
    }
    if (sqlite3_bind_int(stmt, 4, vid) != SQLITE_OK)
    {
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
    }
    if (sqlite3_bind_double(stmt, 5, value) != SQLITE_OK)
    {
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
    }
    if (sqlite3_step(stmt) == SQLITE_ERROR)
    {
      delete [] df;
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
      throw "Cannot insert atmodata";
    }
    if (sqlite3_reset(stmt) != SQLITE_OK)
    {
      delete [] df;
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
      throw "Cannot insert atmodata";
    }
  }
  delete [] df;
  if (sqlite3_finalize(stmt) != SQLITE_OK)
  {
    sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
    sqlite3_free(sqlerr);
    throw sqlite3_errmsg(db);
  }

  sprintf(sqlbuff, "COMMIT;");
  if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
    sqlite3_free(sqlerr);
    throw "Cannot commit !!";
  }
  return;
}

void grib2_to_atmodata(sqlite3 *db,
                      t_g2_pinfo &pkci,
                      unsigned char *sec7,
                      long scode_id, long vtid, long gid, long vid, long *mid)
{
  float *df = new float[pkci.ndata];
  try
  {
    if (unpk_grib2(df, pkci, sec7) != 0)
    {
      delete [] df;
      throw "Error unpacking GRIB2 data";
    }
  }
  catch (char const* str)
  {
    std::cerr << str << std::endl;
    delete [] df;
    throw "Error unpacking GRIB2 data";
  }

  if (verbose)
  {
    float maxval = FLT_MIN;
    float minval = FLT_MAX;
    for (int i = 0; i < pkci.ndata; i ++)
    {
      if (maxval < df[i]) maxval = df[i];
      if (minval > df[i]) minval = df[i];
    }
    std::cout << "Max value : " << maxval << std::endl;
    std::cout << "Min value : " << minval << std::endl;
  }

  // search gid in gp
  t_gpoints *gpp = 0;
  for (int i = 0; i < (int) gp.size(); i ++)
    if (gp[i].gid == gid) gpp = &(gp[i]);

  if (gpp == 0)
  {
    int nr, nc;
    char **result = 0;
    t_gpoints gpn;
    gpn.gid = gid;
    sprintf(sqlbuff, "SELECT id,ipos FROM place_grid WHERE grid_id=%ld;", gid);
    if (sqlite3_get_table(db, sqlbuff, &result, &nr, &nc, &sqlerr) != SQLITE_OK)
    {
      std::cerr << "Cannot exec query : " << sqlerr << std::endl;
      sqlite3_free(sqlerr);
      throw sqlerr;
    }
    sqlite3_free(sqlerr);
    if (nr == 0)
    {
      delete [] df;
      sqlite3_free_table(result);
      return;
    }
    long pid;
    long ipos;
    for (int i = 0; i < nr; i ++)
    {
      sscanf(result[2+i*2], "%ld", &pid);
      sscanf(result[2+i*2+1], "%ld", &ipos); 
      gpn.pid.push_back(pid);
      gpn.ipos.push_back(ipos);
    }
    gp.push_back(gpn);
    gpp = &(gp.back());
    sqlite3_free_table(result);
  }

  sprintf(sqlbuff, "INSERT OR IGNORE INTO tmpins VALUES(?,?,?,?,?)");
  sqlite3_stmt *stmt = 0;
  if (sqlite3_prepare_v2(db, sqlbuff, -1, &stmt, NULL) != SQLITE_OK)
  {
    std::cerr << "Cannot prepare query." << std::endl;
    throw sqlite3_errmsg(db);
  }
  sprintf(sqlbuff, "BEGIN TRANSACTION;");
  if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    throw "Cannot enter transaction";
  }
  sqlite3_free(sqlerr);

  double value;
  for (int i = 0; i < (int) gpp->pid.size(); i ++)
  {
    value = df[gpp->ipos[i]];
    if (sqlite3_bind_int(stmt, 1, scode_id) != SQLITE_OK)
    {
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
    }
    if (sqlite3_bind_int(stmt, 2, gpp->pid[i]) != SQLITE_OK)
    {
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
    }
    if (sqlite3_bind_int(stmt, 3, vtid) != SQLITE_OK)
    {
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
    }
    if (sqlite3_bind_int(stmt, 4, vid) != SQLITE_OK)
    {
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
    }
    if (sqlite3_bind_double(stmt, 5, value) != SQLITE_OK)
    {
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
    }
    if (sqlite3_step(stmt) == SQLITE_ERROR)
    {
      delete [] df;
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
      throw "Cannot insert atmodata";
    }
    if (sqlite3_reset(stmt) != SQLITE_OK)
    {
      delete [] df;
      sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
      sqlite3_free(sqlerr);
      throw "Cannot insert atmodata";
    }
  }
  delete [] df;
  if (sqlite3_finalize(stmt) != SQLITE_OK)
  {
    sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
    sqlite3_free(sqlerr);
    throw sqlite3_errmsg(db);
  }

  sprintf(sqlbuff, "COMMIT;");
  if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    sqlite3_exec(db, "ROLLBACK", NULL, NULL, &sqlerr);
    sqlite3_free(sqlerr);
    throw "Cannot commit !!";
  }
  return;
}
