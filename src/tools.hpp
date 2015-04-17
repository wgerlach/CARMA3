
/******************************************************************************
 *    Copyright (C) 2009 by Wolfgang Gerlach                                  *
 *                                                                            *
 *    E-mail: wgerlachXcebitec.uni-bielefeld.de , where X=@                   *
 *            mailXwolfgang-gerlach.com                                       *
 *                                                                            *
 *    This file is part of CARMA3.                                            *
 *                                                                            *
 *    CARMA3 is free software: you can redistribute it and/or modify          *
 *    it under the terms of the GNU General Public License as published by    *
 *    the Free Software Foundation, either version 3 of the License, or       *
 *    (at your option) any later version.                                     *
 *                                                                            *
 *    CARMA3 is distributed in the hope that it will be useful,               *
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *    GNU General Public License for more details.                            *
 *                                                                            *
 *    You should have received a copy of the GNU General Public License       *
 *    along with CARMA3.  If not, see <http://www.gnu.org/licenses/>.         *
/******************************************************************************/


#ifndef GUARD_TOOLS
#define GUARD_TOOLS




#include "global.hpp"


using namespace std;



class BadConversion : public std::runtime_error {
 public:
   BadConversion(std::string const& s)
     : std::runtime_error(s)
     { }
 };


int str2int(string& text);

void checkFileWarning(string * file);

bool FileExists(const char * filename);
bool DirExists(char * dirname);

void checkFile(string * file);

void stringToLower(string& myString);
void stringToUpper(string& myString);

bool copyFile(const char * sourcefile, const char * destfile);

string stringify(double x);

string double2str(double number);
double string2double(string& text);

void deleteFile(const char * file);

double abs_double(double x);

int mycompare(string& line, const char * t);
int mycompare(string& line, int i, int j, const char * t);



int str2int (const string &str);

string int2str (int n);


vector<string> * parse_column_data_line(string& line, char seperator);
vector<string> * parse_column_data_line(string& line, string seperator_string);

size_t getFreeDiskSpace(char * mountpoint);

struct delete_object
{
  template <typename T>
  void operator()(T *ptr){ delete ptr;}
};



template <class T>
void DeleteContainerWithPointers (T * container) {


    if (container == 0) return;
	for_each( container->begin(), container->end(), delete_object() );
	delete container;

	return;
}


// Functor for deleting pointers in map.
template<class A, class B>
struct DeleteMapFntor
{
    // Overloaded () operator.
    // This will be called by for_each() function.
    bool operator()(pair<A,B> x) const
    {
        // Assuming the second item of map is to be
        // deleted. Change as you wish.
        delete x.second;
        return true;
    }
};

template <class A, class B>
void DeleteMapWithPointers (map< A, B> * my_map) {


    if (my_map == 0) return;
	for_each( my_map->begin(), my_map->end(), DeleteMapFntor<A, B>() );
	delete my_map;

	return;
};

#endif
