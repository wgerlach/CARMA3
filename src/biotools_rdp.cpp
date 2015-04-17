
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

#include "biotools_rdp.hpp"



void parse_RDP_description_line(string& line, string ** identifier, tax_id* taxid, int * offset, map<string, tax_id> * NameToNCBI ) {

	size_t space = line.find_first_of(' ');
	if (space == string::npos) {
			cerr << "parsing error, no space found: \"" << line << "\"" << endl;
			exit(1);
	}


	if (space > 3) {
		if ( identifier != 0) {
			*identifier = new string(line.substr(1, space - 1));
		}
		//cout << "(" << *identifier << ")"  << endl;
		//exit(0);
	} else {
			cerr << "error parsing alignments" << endl;
			exit(1);
	}


	size_t semi = line.find_first_of(';');
	if (semi == string::npos) {
			cerr << "parsing error, no semicolon found."<< endl;
			exit(1);
	}


	if (taxid != 0) {
		string taxon_string = line.substr(space+1, semi-space -1);
		stringToLower(taxon_string);

		//cout << "(" << taxon_string << ")" << endl;

	if (false) {
		int found_uncultured = taxon_string.find("uncultured", 0);
		if (found_uncultured != string::npos) {
			*taxid=0;
			return;
			//uncultured_count++;

			//continue;
		}

		int found_unidentified = taxon_string.find("unidentified", 0);
		if (found_unidentified != string::npos) {
			*taxid=0;
			return;
			//unidentified_count++;

			//continue;
		}
	}

		*taxid = taxonname2taxid(taxon_string, NameToNCBI);

	}

	if (offset != 0) {
		size_t last_sep = line.find_last_of(';');
		if (last_sep == string::npos) {
			cerr << "error parsing description: "<< "\""<< line << "\"" << endl;
			exit(1);
		}

		if (line.compare(last_sep, 8 , "; OFFSET") != 0 ) {
			cerr << "error parsing description: " << "\""<< line << "\"" << endl;
			cerr << "OFFSET information not found" << endl;
			exit(1);
		}
		string offset_str = line.substr(last_sep+8);
		if (offset_str.length() == 0) {
			cerr << "error: offset_str.length() == 0" << endl;
			exit(1);
		}

		*offset = str2int(offset_str);
	}
}
