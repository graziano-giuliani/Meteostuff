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
#include <vector>
#include <cstring>
#include <cstdio>

#include <sqlite3.h>
#include <places.h>
#include <earthmap.h>

using namespace himet;

static const char gribdb[64] = DBFILE;
static char sqlbuff[2048];
char *sqlerr = 0;

char placetable[64] = "place";

bool verbose = false;

int main(int argc, char *argv[])
{
  if (argc < 2) return -1;

  if (argc == 3)
  {
    snprintf(placetable, 64, "%s", argv[2]);
    std::cout << "Using place table " << placetable << std::endl;
  }

  placelist pl;
  pl.read(argv[1]);
  int nplaces = pl.size( );
  if (nplaces == 0)
  {
    std::cerr << "Nothing to be done. No place parsed from input file."
              << std::endl;
    return -1;
  }
  std::cout << "Found a total of " << nplaces << " places." << std::endl;
  if (verbose)
   for (int i = 0; i < nplaces; i ++)
     std::cout << pl.at(i) << std::endl;

  // Connect to DB
  if (sqlite3_initialize() != SQLITE_OK)
  {
    std::cerr << "Database access impossible." << std::endl;
    return -1;
  }

  sqlite3 *db;
  if (sqlite3_open(gribdb, &db) != SQLITE_OK)
  {
    std::cerr << "Database file " << gribdb << " not ready." << std::endl;
    return -1;
  }

  std::vector <long> new_place_id;
  try
  {
    int nr, nc;
    char **result;

    for (int ip = 0; ip < nplaces; ip ++)
    {
      nr = nc = 0;
      sprintf(sqlbuff, "SELECT id FROM place WHERE name='%s' AND"
          " adm1='%s' AND adm2='%s';", pl.at(ip).name.c_str(),
          pl.at(ip).admin1.c_str(), pl.at(ip).admin2.c_str());
      if (sqlite3_get_table(db, sqlbuff, &result, &nr, &nc, &sqlerr) !=
             SQLITE_OK)
      {
        std::cerr << "Cannot exec query : " << sqlerr << std::endl;
        sqlite3_free_table(result);
        sqlite3_free(sqlerr);
        throw "DB error";
      }
      sqlite3_free(sqlerr);
      if (nr == 1)
      {
        sqlite3_free_table(result);
        std::cerr << "Skipped place " << pl.at(ip).name
                  << " because already in the DB." << std::endl;
        continue;
      }

      if (verbose)
        std::cout << "Inserting " << pl.at(ip).name << " in the DB"
                  << std::endl;

      // Insert place and retrieve id
      sprintf(sqlbuff, "INSERT INTO %s(name, adm1, adm2, longitude, "
          "latitude, altitude) VALUES('%s','%s','%s',%f,%f,%f);",
          placetable, pl.at(ip).name.c_str(), pl.at(ip).admin1.c_str(),
          pl.at(ip).admin2.c_str(),
          pl.at(ip).lon, pl.at(ip).lat, pl.at(ip).hgt);
      if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
      {
        std::cerr << "Cannot exec query : " << sqlerr << std::endl;
        sqlite3_free(sqlerr);
        throw "DB error";
      }
      sqlite3_free(sqlerr);

      sprintf(sqlbuff, "SELECT id FROM %s WHERE name='%s' AND"
          " latitude=%f AND longitude=%f;", placetable, pl.at(ip).name.c_str(),
          pl.at(ip).lat,pl.at(ip).lon);
      if (sqlite3_get_table(db, sqlbuff, &result, &nr, &nc, &sqlerr) !=
             SQLITE_OK)
      {
        std::cerr << "Cannot exec query : " << sqlerr << std::endl;
        sqlite3_free_table(result);
        sqlite3_free(sqlerr);
        throw "DB error";
      }
      sqlite3_free(sqlerr);
      if (nr != 1)
      {
        std::cerr << "DB Error : Cannot find inserted place !!" << std::endl;
        sqlite3_free_table(result);
        throw "DB error";
      }
      long pid;
      if (sscanf(result[1], "%ld", &pid) != 1)
      {
        std::cerr << "Not parsable place stored in DB : " << result[1]
                  << std::endl;
        sqlite3_free_table(result);
        throw "DB error";
      }
      new_place_id.push_back(pid);     
    }
  }
  catch (char const* str)
  {
    std::cerr << str << std::endl;
    if (sqlite3_shutdown() != SQLITE_OK)
    {
      std::cerr << "Database shutdown error !!!" << std::endl;
      return -1;
    }
    return -1;
  }

  if (new_place_id.size() == 0)
  {
    std::cout << "No new place to insert into the DB." << std::endl;
    if (sqlite3_shutdown() != SQLITE_OK)
    {
      std::cerr << "Database shutdown error !!!" << std::endl;
      return -1;
    }
    return 0;
  }

  try
  {
    int nr, nc;
    char **result = 0;

    // Load GRID list from DB
    nr = nc = 0;
    sprintf(sqlbuff, "%s", "SELECT id,emap from grid;");
    if (sqlite3_get_table(db, sqlbuff, &result, &nr, &nc, &sqlerr) != SQLITE_OK)
    {
      std::cerr << "Cannot exec query : " << sqlerr << std::endl;
      sqlite3_free(sqlerr);
      throw sqlerr;
    }
    sqlite3_free(sqlerr);

    if (nr == 0)
    {
      sqlite3_free_table(result);
      std::cout << "No grid defined, so no need to calculate position"
                << std::endl;
      if (sqlite3_shutdown() != SQLITE_OK)
      {
        std::cerr << "Database shutdown error !!!" << std::endl;
        return -1;
      }
      return 0;
    }

    long idgrid;
    char *emap;
    for (int igrid = 0; igrid < nr; igrid ++)
    {
      sscanf(result[(igrid+1)*2], "%ld", &idgrid);
      emap = result[(igrid+1)*2+1];

      if (verbose)
        std::cout << "GRID " << idgrid << " : " << emap << std::endl;

      earthmap m;
      if (strncmp(emap, "regll", 5) == 0)
      {
        if (verbose) std::cout << "Regular lat lon" << std::endl;
        t_proj_latlon_parameters p;
        sscanf(emap, "regll,%ldx%ld,%lf,%lf,%lf,%lf,%d",
          &p.nx, &p.ny, &p.lat1, &p.lon1, &p.lat2, &p.lon2, &p.scanflag);
        m.set(p);
      }
      else if (strncmp(emap, "rll", 3) == 0)
      {
        if (verbose) std::cout << "Rotated lat lon" << std::endl;
        t_proj_rotlat_parameters p;
        sscanf(emap, "rll,%ldx%ld,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d",
          &p.nx, &p.ny, &p.lat1, &p.lon1, &p.lat2, &p.lon2, &p.polelat,
          &p.polelon,&p.rotation,&p.scanflag);
        m.set(p);
      }
      else if (strncmp(emap, "merc", 4) == 0)
      {
        if (verbose) std::cout << "Mercator" << std::endl;
        t_proj_mercator_parameters p;
        sscanf(emap, "merc,%ldx%ld,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d",
           &p.nx, &p.ny, &p.lat1, &p.lon1, &p.lat2, &p.lon2, &p.dx, &p.dy,
           &p.reflat, &p.scanflag);
        m.set(p);
      }
      else if (strncmp(emap, "ps", 2) == 0)
      {
        if (verbose) std::cout << "Polar Stereographic" << std::endl;
        t_proj_polar_parameters p;
        sscanf(emap, "ps,%ldx%ld,%lf,%lf,%lf,%lf,%lf,%d,%d,%d",
          &p.nx, &p.ny, &p.lat1, &p.lon1, &p.orient, &p.dx, &p.dy,
          &p.resflag, &p.projflag, &p.scanflag);
        m.set(p);
      }
      else if (strncmp(emap, "lcc", 3) == 0)
      {
        if (verbose) std::cout << "Lambert Conformal" << std::endl;
        t_proj_lambert_parameters p;
        sscanf(emap, "lcc,%ldx%ld,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d,%d,%d",
         &p.nx, &p.ny, &p.lat1, &p.lon1, &p.orient, &p.dx, &p.dy, &p.truelat1,
         &p.truelat2, &p.resflag, &p.projflag, &p.scanflag);
        m.set(p);
      }
      else if (strncmp(emap, "ss2dE", 5) == 0)
      {
        if (verbose) std::cout << "Rotated Arakawa E 2D" << std::endl;
        t_proj_rotwrf_parameters p;
        sscanf(emap, "ss2dE,%ldx%ld,%lf,%lf,%lf,%lf,%d,%d",
          &p.nx, &p.ny, &p.lat1, &p.lon1, &p.lat0, &p.lon0,
          &p.resflag, &p.scanflag);
        m.set(p);
      }
      for (int ip = 0; ip < nplaces; ip ++)
      {
        if (verbose)
          std::cout << "Processing place " << ip+1 << std::endl;

        long pid = new_place_id[ip];
        double pixel, line, crot, srot, locres;
        double lat, lon, ri, rj, latp, lonp, latm, lonm;
        lat = pl.at(ip).lat;
        lon = pl.at(ip).lon;
        m.latlon_to_ij(lat, lon, &ri, &rj);
        m.ij_to_latlon(floor(ri), floor(rj), &latm, &lonm);
        m.ij_to_latlon(ceil(ri), ceil(rj), &latp, &lonp);
        locres = sqrt(pow(latp-latm, 2.0)+pow(lonp-lonm, 2.0));
        m.ij2pl(ri, rj, &pixel, &line);
        m.rotuv_ll2gr(lat, lon, &crot, &srot);
        if (pixel > 1.0e20) continue;
        long ipos = m.ij2p((int) rintf(pixel), (int) rintf(line));
        sprintf(sqlbuff, "INSERT OR IGNORE INTO %s_grid("
          "place_id,grid_id,pixel,line,ipos,crot,srot,locres) "
          "values (%ld,%ld,%f,%f,%ld,%f,%f,%f);", placetable,
          pid, idgrid, pixel, line, ipos, crot, srot, locres);
        if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
        {
          std::cerr << "Cannot exec query : " << sqlerr << std::endl;
          sqlite3_free(sqlerr);
          throw "DB error";
        }
        sqlite3_free(sqlerr);
      }
    }
    sqlite3_free_table(result);
  }
  catch (char const* str)
  {
    std::cerr << str << std::endl;
    if (sqlite3_shutdown() != SQLITE_OK)
    {
      std::cerr << "Database shutdown error !!!" << std::endl;
      return -1;
    }
    return -1;
  }

  if (sqlite3_shutdown() != SQLITE_OK)
  {
    std::cerr << "Database shutdown error !!!" << std::endl;
    return -1;
  }
  std::cout << "Done." << std::endl;
  return 0;
}
