
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

#include "tools.hpp"

int str2int(string& text){

	int number = std::atoi( text.c_str() );

	return number;

}


void checkFileWarning(string * file) {
	if (file == 0) {
		return;
	}

	if (! FileExists((char *) file->c_str())) {
		cerr << "WARNING: File " << *file << " not found. If you don't need this file, ignore this warning."<< endl;
		return;
	}

}

bool copyFile(const char * sourcefile, const char * destfile) {

	namespace bfs=boost::filesystem;
//	cerr << "copy: " << sourcefile << " to " << destfile << endl;

//	if (FileExists(sourcefile)) {
//		cerr << "source exists" << endl;
//	} else {
//		cerr << "source does not exist" << endl;
//	}
//
//	if (FileExists(destfile)) {
//		cerr << "destfile exists" << endl;
//	} else {
//		cerr << "destfile does not exist" << endl;
//	}

//	string cp_command = *cp_bin;
//	cp_command.append(sourcefile);
//	cp_command.append(" ");
//	cp_command.append(destfile);
//	int cp_ret = system(cp_command.c_str());

	bfs::path source (sourcefile);
	bfs::path dest (destfile);
	try	{
		bfs::copy_file(source, dest, bfs::copy_option::overwrite_if_exists);
	} catch(std::exception e) {
		cerr << "exception:" << e.what() << endl;
		return false;
	}

//	if (FileExists(destfile)) {
//		cerr << "2destfile exists" << endl;
//	} else {
//		cerr << "2destfile does not exist" << endl;
//	}

//	if (cp_ret != 0) {
//		return false;
//	}
	return true;
}

bool copyFile_old(const char * sourcefile, const char * destfile) {

	ifstream srce( sourcefile, ios::binary );

	if( !srce ) {
		cerr << "could not open " << sourcefile << " for reading\n" ;
		return false ;
	}

	ofstream dest( destfile, ios::binary ) ;

	if( !dest ){
		cerr << "could not open " << destfile << " for writing\n" ;
		return false ;
	}

	dest << srce.rdbuf() ;

	if( !dest ){

		cerr << "error while copying\n" ;

		return false ;

	}
	return true;
}


bool DirExists(char * dirname){

	struct stat st;
	if(stat(dirname,&st) == 0) {
		return true;
	}

    return false;
}

bool FileExists(const char * filename) {
  struct stat stFileInfo;
  bool blnReturn;
  int intStat;

  // Attempt to get the file attributes
  intStat = stat(filename,&stFileInfo);
  if(intStat == 0) {
    // We were able to get the file attributes
    // so the file obviously exists.
    blnReturn = true;
  } else {
    // We were not able to get the file attributes.
    // This may mean that we don't have permission to
    // access the folder which contains this file. If you
    // need to do that level of checking, lookup the
    // return values of stat which will give you
    // more details on why stat failed.
    blnReturn = false;
  }

  return(blnReturn);
}

void checkFile(string * file) {
	if (file == 0) {
		return;
	}

	if (! FileExists((char *) file->c_str())) {
		cerr << "error reading file: " << *file << endl;
		exit(EXIT_FAILURE);
	}

}

void stringToLower(string& myString) {
	std::transform(myString.begin(), myString.end(), myString.begin(),
			   (int(*)(int)) std::tolower);
}

void stringToUpper(string& myString) {
	std::transform(myString.begin(), myString.end(), myString.begin(),
			   (int(*)(int)) std::toupper);
}

string double2str(double number) {


	stringstream NumberString;
	NumberString << number;
	string number_str = NumberString.str();

	return number_str;
}


double string2double(string& text){
      std::istringstream stm;
      stm.str(text);
      double d;
      stm >>d;
      return d;
}


string stringify(double x) {
   std::ostringstream o;
   if (!(o << x))
     throw BadConversion("stringify(double)");
   return o.str();
 }


void deleteFile(const char * file) {
	if( remove( file ) != 0 ) {
		cerr << "Error deleting file " << file << endl;
	}

	return;
}



double abs_double(double x){
	return x<0?-x:x;
}


int mycompare(string& line, const char * t) {
	int r;
	try {
		r = line.compare(t);

	} catch (std::out_of_range) {
		return 1;
	}

	return r;

}

int mycompare(string& line, int i, int j, const char * t) {
	int r;
	try {
		r = line.compare(i,j, t);
//cout << "r: " << r << endl;
//cout << line.substr(i, j) << endl;
//cout << t << endl;
	} catch (std::out_of_range) {
		return 1;
	}

	return r;
}


vector<string> * parse_column_data_line(string& line, char seperator) {

	vector<string> * results = new vector<string>;

	size_t left = 0;
	size_t right = 0;


	right = line.find_first_of(seperator);
	if (right == string::npos) {
		results->push_back(line);
		return results;
		//cerr << "error: right == string::npos" << endl;
		//cerr << "line:" << line << endl;
		//exit(EXIT_FAILURE);
	}

	string entry = line.substr(left, right-left);
	results->push_back(entry);

	bool loop = true;
	while ( loop ) {

		// next start
		left = line.find_first_not_of(seperator, right);
		if (left == string::npos) {
			cerr << "error: left == string::npos" << endl;
			exit(EXIT_FAILURE);
		}

		// next end
		right = line.find_first_of(seperator, left);
		if (right == string::npos) {
			right = line.length();
			loop = false;
		}
		entry = line.substr(left, right-left);
		//cout << "got:_" << entry << "_"<< endl;
		results->push_back(entry);
	}


	//cout << "got:_" << results->at(0) << "_"<< endl;


//			// end of tlen:
//			left = line.find_first_of(' ', left);
//			if (left == string::npos) {
//					cerr << "blubb" << endl;
//					exit(EXIT_FAILURE);
//			}


	return results;
}


vector<string> * parse_column_data_line(string& line, string seperator_string) {

	vector<string> * results = new vector<string>;

	size_t prev = 0;
	size_t next = 0;

	size_t sep_len = seperator_string.length();


	string entry;
	//results->push_back(entry);


	while ( true ) {
		prev = next;

		// next start
		next = line.find(seperator_string, prev+1);
		if (next == string::npos) {
			break;
		}


		entry = line.substr(prev, next-prev);
		//cout << "got:_" << entry << "_"<< endl;
		results->push_back(entry);
	}

	entry = line.substr(prev);
	results->push_back(entry);

	//cout << "got:_" << results->at(0) << "_"<< endl;


//			// end of tlen:
//			left = line.find_first_of(' ', left);
//			if (left == string::npos) {
//					cerr << "blubb" << endl;
//					exit(EXIT_FAILURE);
//			}


	return results;
}

size_t getFreeDiskSpace(char * mountpoint) {

	struct statvfs fiData;
	struct statvfs *fpData;
	//char fnPath[128];
   // int i;

//        if( argc < 2 ) {
//                printf("Usage, webdisk DEVICE0 ..... DEVICEX\n");
//                return(2);
//        }

	//Lets loopyloop through the argvs
	//for( i = 1 ; i<argc; i++ ) {
	//strcpy(fnPath, argv[i]);
	if((statvfs(mountpoint,&fiData)) < 0 ) {
			//printf("Failed to stat %s:\n", mountpoint);
			cerr << "Failed to stat: " << mountpoint << endl;
			exit(1);
	} else {
			//printf("Disk %s: \n", fnPath);

			//printf("\tblock size: %u\n", fiData.f_bsize);
			//cout << "block size: " << fiData.f_bsize << endl;
			//printf("\ttotal no blocks: %i\n", fiData.f_blocks);
			//cout << "total no blocks: " << fiData.f_blocks << endl;
			//printf("\tfree blocks: %i\n", fiData.f_bfree);
			//cout << "free blocks: " <<fiData.f_bfree << endl;

			size_t disc_size = ( fiData.f_bsize * fiData.f_bfree );
			return disc_size;
			//cout << size <<  " MB" << endl;
	}
	//}

}


/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain.
 */

typedef double elem_type ;
#define ELEM_SWAP(a,b) { register elem_type t=(a);(a)=(b);(b)=t; }

elem_type quick_select(elem_type arr[], int n)
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return arr[median] ;

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]) ;
            return arr[median] ;
        }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]) ;

    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]) ;

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
        do ll++; while (arr[low] > arr[ll]) ;
        do hh--; while (arr[hh]  > arr[low]) ;

        if (hh < ll)
        break;

        ELEM_SWAP(arr[ll], arr[hh]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]) ;

    /* Re-set active partition */
    if (hh <= median)
        low = ll;
        if (hh >= median)
        high = hh - 1;
    }
}

int str2int (const string &str) {
	stringstream ss(str);
	int n;
	ss >> n;
	return n;
}

string int2str (int n) {
	stringstream ss;
	ss << n;
	return ss.str();
}
