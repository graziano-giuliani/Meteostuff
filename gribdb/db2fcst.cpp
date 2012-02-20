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
#include <sstream>
#include <cstring>

#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/local_time/local_time.hpp>
#include <boost/date_time/time_facet.hpp>

#include <sqlite3.h>

#include <weather.h>
#include <suntime.h>

using namespace boost::gregorian;
using namespace boost::posix_time;
using namespace himet;

static const char gribdb[64] = DBFILE;
static char sqlbuff[2048];
char *sqlerr = 0;

const static int nvars = 10;

static char varname[nvars][32] = { "SP", "SP", "SP", "SP", "SP",
                                   "SP", "SP", "SP", "SP", "SP" };

typedef struct {
  float X[nvars];
} t_fastvar;

typedef struct
{
  char *id;
  char *name;
  char *adm1;
  char *adm2;
  char *lat;
  char *lon;
  char *hgt;
  std::vector <t_fastvar> meteo;
} t_fastplace;

bool verbose = true;
bool mega_verbose = false;

int main(int argc, char *argv[])
{
  if (argc < 5)
  {
    std::cerr << "Not enough arguments." << std::endl << std::endl
              << "Usage : " << std::endl
              << "   " << argv[0] << " place_db_table fcst_start fcst_end step"
              << std::endl << std::endl << "Times are formatted as "
              "YYYYMMDDTHHMMSS and step is in hours." << std::endl;
    return -1;
  }

  ptime tstart;
  ptime tend;
  time_duration td;
  char *place_table;
  std::vector <t_fastplace> pl;

  try
  {
    place_table = strdup(argv[1]);
    tstart = from_iso_string(argv[2]);
    tend = from_iso_string(argv[3]);
    int step;
    if (sscanf(argv[4], "%d", &step) != 1)
      throw "placeholder";
    td = hours(step);
    std::cout << "###############################################" << std::endl
              << "Read places from table : " << place_table << std::endl
              << "Forecast time start is : " 
              << to_simple_string(tstart) << std::endl
              << "Forecast time end is   : " 
              << to_simple_string(tend) << std::endl
              << "Forecast step is       : "
              << step << " hours." << std::endl
              << "###############################################" << std::endl;
  }
  catch (...)
  {
    std::cerr << "Cannot parse input arguments." << std::endl;
    return -1;
  }

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

  sprintf(sqlbuff, "%s", "PRAGMA journal_mode = MEMORY;");
  if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    return -1;
  }
  sqlite3_free(sqlerr);
  sprintf(sqlbuff, "%s", "PRAGMA cache_size = 32768;");
  if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    return -1;
  }
  sqlite3_free(sqlerr);
  sprintf(sqlbuff, "%s", "PRAGMA page_size = 16384;");
  if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    return -1;
  }
  sqlite3_free(sqlerr);

  std::cout << "Creating needed indexes. Please wait...";
  sprintf(sqlbuff, "%s",
          "CREATE INDEX atmodata_placegrid_index on atmodata(plgr_id);");
  if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    return -1;
  }
  sqlite3_free(sqlerr);
  std::cout << "...";
  sprintf(sqlbuff, "%s",
          "CREATE INDEX atmodata_var_index on atmodata(var_id);");
  if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    return -1;
  }
  sqlite3_free(sqlerr);
  std::cout << "...";
  sprintf(sqlbuff, "%s",
          "CREATE INDEX atmodata_vtime_index on atmodata(vtime_id);");
  if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    return -1;
  }
  sqlite3_free(sqlerr);
  std::cout << "...";
  sprintf(sqlbuff, "%s",
          "CREATE INDEX atmodata_source_index on atmodata(source_id);");
  if (sqlite3_exec(db, sqlbuff, NULL, NULL, &sqlerr) != SQLITE_OK)
  {
    std::cerr << "Cannot exec query : " << sqlerr << std::endl;
    sqlite3_free(sqlerr);
    return -1;
  }
  sqlite3_free(sqlerr);
  std::cout << std::endl;
  std::cout << "Done indexing." << std::endl;
  return 0;

  int nr, nc;
  char **dbpl;

  // Build up the list of place ids.
  try
  {
    sprintf(sqlbuff, "SELECT * from %s", place_table);
    if (sqlite3_get_table(db, sqlbuff, &dbpl, &nr, &nc, &sqlerr) != SQLITE_OK)
    {
      std::cerr << "Cannot exec query : " << sqlerr << std::endl;
      sqlite3_free(sqlerr);
      throw "DB error";
    }
    sqlite3_free(sqlerr);
    t_fastplace p;
    for (int i = 0; i < nr; i ++)
    {
      p.id = dbpl[7+i*7];
      p.name = dbpl[7+i*7+1];
      p.adm1 = dbpl[7+i*7+2];
      p.adm2 = dbpl[7+i*7+3];
      p.lon = dbpl[7+i*7+4];
      p.lat = dbpl[7+i*7+4];
      p.hgt = dbpl[7+i*7+4];
      if (mega_verbose)
        std::cout << p.name << std::endl;
      pl.push_back(p);
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

  if (pl.size() == 0)
  {
    std::cerr << "No place defined. Good bye!" << std::endl;
    if (sqlite3_shutdown() != SQLITE_OK)
    {
      std::cerr << "Database shutdown error !!!" << std::endl;
      return -1;
    }
    return 0;
  }

  if (verbose)
    std::cout << "Total number of places to make forecast is " << pl.size()
              << std::endl;

  std::vector <ptime> lot;
  ptime treq = tstart;
  do
  {
    lot.push_back(treq);
    treq = treq+td;
  } while (treq <= tend);

  if (verbose)
    std::cout << "Total number of time intervals to search for is "
              << lot.size() << std::endl;
  if (mega_verbose)
    for (int i = 0; i < (int) lot.size(); i ++)
      std::cout << "Time " << i+1 << " : "
                << to_simple_string(lot[i]) << std::endl;

  // retrieve data
  try
  {
    t_fastvar x;

    std::stringstream ss;
    time_facet *output_facet = new time_facet();
    time_input_facet *input_facet = new time_input_facet();
    ss.imbue(std::locale(ss.getloc(), output_facet));
    ss.imbue(std::locale(ss.getloc(), input_facet));
    output_facet->format("%Y-%m-%d %H:%M:%S");
    input_facet->format("%Y-%m-%d %H:%M:%S");

    sprintf(sqlbuff, "SELECT atmodata.value FROM atmodata "
      "INNER JOIN place_grid ON atmodata.plgr_id=place_grid.id "
      "INNER JOIN valtime ON atmodata.vtime_id=valtime.id "
      "INNER JOIN reftime ON valtime.ref_id=reftime.id "
      "INNER JOIN atmovar ON atmodata.var_id=atmovar.id "
      "WHERE atmovar.name=? AND place_grid.place_id=? "
      "AND valtime.valid1=? "
      "ORDER BY place_grid.locres ASC, reftime.reference ASC LIMIT 1;");
    
    sqlite3_stmt *stmt = 0;
    if (sqlite3_prepare_v2(db, sqlbuff, -1, &stmt, NULL) != SQLITE_OK)
    {
      std::cerr << "Cannot prepare query." << std::endl;
      throw sqlite3_errmsg(db);
    }

    const char *idp;
    const char *vtime;
    const char *vname;

    for (int iv = 0; iv < nvars; iv ++)
    {
      for (int ip = 0; ip < (int) pl.size(); ip ++)
      {
        for (int it = 0; it < (int) lot.size(); it ++)
        {
          idp = pl[ip].id;
          ss.str("");
          ss << lot[it];
          vtime = ss.str().c_str();
          vname = varname[iv];
        
          if (sqlite3_bind_text(stmt, 1, vname, -1, SQLITE_STATIC) != SQLITE_OK)
          {
            std::cerr << "Cannot bind Variable name" << std::endl;
            throw sqlite3_errmsg(db);
          }
          if (sqlite3_bind_text(stmt, 2, idp, -1, SQLITE_STATIC) != SQLITE_OK)
          {
            std::cerr << "Cannot bind point index" << std::endl;
            throw sqlite3_errmsg(db);
          }
          if (sqlite3_bind_text(stmt, 3, vtime, -1, SQLITE_STATIC) != SQLITE_OK)
          {
            std::cerr << "Binding date requested" << std::endl;
            throw sqlite3_errmsg(db);
          }

          if (sqlite3_step(stmt) == SQLITE_ERROR)
            throw sqlite3_errmsg(db);

          if (sqlite3_data_count(stmt) != 1)
          {
            x.X[iv] = -999.99;
            pl[ip].meteo.push_back(x);
            sqlite3_reset(stmt);
            continue;
          }

          x.X[iv] = sqlite3_column_double(stmt, 0);
          pl[ip].meteo.push_back(x);

          if (sqlite3_reset(stmt) != SQLITE_OK)
            throw sqlite3_errmsg(db);
        }
      }
    }
    if (sqlite3_finalize(stmt) != SQLITE_OK)
          throw sqlite3_errmsg(db);
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

  // Do the work

  // Cleanup
  sqlite3_free_table(dbpl);
  if (sqlite3_shutdown() != SQLITE_OK)
  {
    std::cerr << "Database shutdown error !!!" << std::endl;
    return -1;
  }
  std::cout << "Done." << std::endl;
  return 0;
}
