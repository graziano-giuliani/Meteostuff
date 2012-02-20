/***************************************************************************
 *   Copyright (C) 2008-2009 by Graziano Giuliani                          *
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

#ifndef __PLACES_H__
#define __PLACES_H__

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sqlite3.h>

namespace himet {

  class place {
    public:
      place( );
      int id;
      float lon;
      float lat;
      float hgt;
      std::string name;
      std::string admin1;
      std::string admin2;
      friend std::ostream& operator<< (std::ostream &out, place &p)
      {
        out << "Location ID " << p.id << " ####" << std::endl
            << "  Name                    : " << p.name << std::endl
            << "  Altitude MSL            : " << p.hgt << std::endl
            << "  Cohordinates (lon, lat) : " << p.lon << ", "
            << p.lat << std::endl
            << "  Administrative info 1   : " << p.admin1 << std::endl
            << "  Administrative info 2   : " << p.admin2;
        return out;
      }
  };

  class placelist {
    public:
      placelist( ) { }
      placelist(const char *filename);
      placelist(sqlite3 *db);
      void read(const char *filename);
      place &search(const int id);
      place &search(const std::string name);
      int size( ) { return list.size(); }
      place &at(const int i);
      void add(place &p);
    private:
      std::vector<place> list;
      std::map<int, int> item_by_id;
  };
}

#endif
