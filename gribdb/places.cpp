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

#include <places.h>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace himet;

std::vector<std::string> tokenize(const std::string& str,
                                  const std::string& delimiters);
bool compare_id(place a, place b);
bool compare_name(place a, place b);

std::vector<std::string> tokenize(const std::string& str,
                                  const std::string& delimiters)
{
  std::vector<std::string> tokens;
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  std::string::size_type pos = str.find_first_of(delimiters, lastPos);

  while (std::string::npos != pos || std::string::npos != lastPos)
  {
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    lastPos = str.find_first_not_of(delimiters, pos);
    pos = str.find_first_of(delimiters, lastPos);
  }
  return tokens;
}

placelist::placelist(const char *filename)
{
  this->read(filename);
}

void placelist::add(place &p)
{
  list.push_back(p);
  item_by_id.insert(std::pair<int,int>(p.id, list.size()));
  std::sort(list.begin(), list.end(), compare_id);
  return;
}

void placelist::read(const char *filename)
{
  int counter = 0;
  std::fstream file_op(filename, std::ios::in);
  if (file_op.good())
  {
    char str[128];
    place p;
    while (! file_op.eof())
    {
      file_op.getline(str, 128);
      std::vector<std::string> tok = tokenize(str, ",");
      if (tok.size() == 7)
      {
        std::istringstream s;
        s.str(tok[0]);
        s >> p.id;
        s.clear();
        p.name = tok[1];
        p.admin1 = tok[2];
        p.admin2 = tok[3];
        s.str(tok[4]);
        s >> p.lat;
        s.clear();
        s.str(tok[5]);
        s >> p.lon;
        s.clear();
        s.str(tok[6]);
        s >> p.hgt;
        s.clear();
        list.push_back(p);
        item_by_id.insert(std::pair<int,int>(p.id, counter));
        counter ++;
      }
    }
  }
  file_op.close();
  std::sort(list.begin(), list.end(), compare_id);
}

bool compare_id(place a, place b) { return (a.id < b.id); }
bool match_id(place a, place b) { return (a.id == b.id); }

place::place( )
{
  id = -1;
  lon = 0.0;
  lat = 0.0;
  hgt = 0.0;
  name = "UNKNOWN";
  admin1 = "UNKNOWN";
  admin2 = "UNKNOWN";
}

place &placelist::search(const int id)
{
  std::map<int, int>::iterator iter = item_by_id.find(id);
  if (iter != item_by_id.end())
    return list[iter->second];
  else
    throw "Item not found";
}

place &placelist::at(const int i)
{
  return list[i];
}

#ifdef TESTME

int main(int argc, char *argv[])
{
  place p;
  placelist pl("luoghi.himet.all");
  std::cout << "Found a total of " << pl.size( ) << " places." << std::endl;
  p = pl.search(4577);
  std::cout << p << std::endl;
}

#endif
