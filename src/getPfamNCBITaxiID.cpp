/***************************************************************************
 *   Copyright (C) 2008 by Wolfgang Gerlach						  *
 *														   *
 *   E-mail: wgerlachXcebitec.uni-bielefeld.de , where X=@			  *
 ***************************************************************************/


#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "limits.h"
#include <vector>
#include <map>



using namespace std;


int str2int(string text){
	int i;
	istringstream myStream(text);
	myStream>>i;
	return i;
}

vector<string> * parseTabSepLine(string& line) {

	vector<string> * results = new vector<string>;

	// cout << "a" << endl;

	string entry;
	size_t left=1;
	size_t right=0;
	while (1) {
		right = line.find ( "\'\t\'", left);
		if (right == string::npos) {
			entry = line.substr(left);
			//cout << "last: " << entry << "." << endl;
			break;
		}
		right += 3;

		entry = line.substr(left, right-left-3);
		//cout << "got: " << entry << "." << endl;
		results->push_back(entry);

		left = right;

	}
	//cout << "b" << endl;
	//exit(0);
	return results;
}

void getPfamToNCBIFrompfamseq(const char * file) {


	int id_column = 1;
	int species_column = 9;
	int ncbi_taxid_column = 15;

	//map<string, int> * PfamToNCBI = new map<string, int>;
	FILE *fp;
	int status;
	int max_line_length = 100000;
	char path[max_line_length+1];

	string zcat_command = "zcat -f ";
	zcat_command.append(file);



	int * intbuffer = new int[250];
	int buf_pos;

	int linecount =0;
	int column = 0;


	fp = popen(zcat_command.c_str(), "r");
	if (fp == NULL) {
		cerr << "error" << endl;
		exit(EXIT_FAILURE);
	}


	map<string, int > * name_to_taxid = new map<string, int >;
	map<string, bool > * amb = new map<string, bool >;
	map<string, bool >::iterator amb_it;
	map<string, int >::iterator name_to_tax_id_it;

	string line;
	int collision = 0;

	//ifstream myfile (file);
	//if (myfile.is_open())
	//{
		//while (! myfile.eof() )	{
		while (fgets(path, max_line_length, fp) != NULL) {
			linecount++;
			//column = 1;

			//cout << "read: [X]" << path << "[X]" <<endl;
			//getline (myfile,line, '\n');
			int i = 0;
			while (1) {
				if (path[i] == '\0') {
					break;
				}
				if (path[i] == '\n') {
					path[i] = '\0';
					break;
				}
				i++;
			}
			line = path;
			int line_len = line.length();
			//cout << "last: " << line[line_len-1] << endl;
			if (line[line_len-1] == '\\') {
				//cout << "short?: " << line << "X" << endl;

				fgets(path, max_line_length, fp);
				int i = 0;
				while (1) {
					if (path[i] == '\0') {
						break;
					}
					if (path[i] == '\n') {
						path[i] = '\0';
						break;
					}
					i++;
				}
				line.append("\n");
				line.append(path);
				//cout << "appended: " << line << endl;
			}




			if (line.length() == 0) {
				continue;
			}

			vector<string> * cols = parseTabSepLine(line);

			if (cols->size() < 16 ) {
				cout << "_" << line << "_" << endl;
				cout << "error, not enough columns..." << endl;
				exit(1);
			}


		 	string name = (*cols)[id_column];

		 	if (false) {
			 	size_t underscore = name.find_first_of ( '_');
				if (underscore == string::npos) {
					cerr << "line: " << path << endl;
					cerr << "parsing error 2" << endl;
					exit(1);
				}
				underscore++;
			 	name = name.substr(underscore);
		 	}

		 	string ncbi = (*cols)[ncbi_taxid_column];
		 	string species = (*cols)[species_column];

		 	//cout << species << endl;


		 	int ncbi_int = str2int(ncbi);

			if (ncbi_int != 0) {
		 		cout << name << "," << ncbi << endl;
		 	}


		 	name_to_tax_id_it = name_to_taxid->find(name);

		 	if (name_to_tax_id_it == name_to_taxid->end()) {
		 		name_to_taxid->insert( pair<string, int >(name, ncbi_int) );
		 	} else {
		 		//cout << name_to_tax_id_it->second << endl;
		 		if (name_to_tax_id_it->second != ncbi_int) {

		 			cout << "warning, unequal! " << collision <<  endl;
		 			cout << name << "," << ncbi << endl;
		 			//exit(1);
		 			//amb_it = amb->find(name);

		 			//if (amb_it == amb->end()) {
			 		//	collision++;
		 			//	amb->insert( pair<string, bool >(name, true) );
		 			//}
			 		//exit(1);
		 		}

		 	}

		 	//buf_pos=0;
		 	//intbuffer[0]=ncbi_int;


//		 	while (ncbi_int > 1) {
//		 		buf_pos++;
//		 		ncbi_int = (*parent_taxid)[ncbi_int];
//		 		intbuffer[buf_pos]=ncbi_int;
//		 		//cout << ncbi_int << endl;
//		 	}
//
//		 	while (buf_pos >= 0) {
//		 		cout << "," << intbuffer[buf_pos];
//			 	buf_pos--;
//		 	}
//		 	cout << endl;

		 	//exit(0);
		 	//} else {
		 	//	cerr << "error" << endl;
		 	//	exit(1);
		 	//}


		 	//if (linecount % 100000 == 0) {

		 	//cout << PfamToNCBI->size() << endl;
		 	//}
	 		delete cols;


		}
	//myfile.close();
	//}

	//else cout << "Unable to open file";

	delete [] intbuffer;


}



map<int, int> * init_NCBI_NODES(const char * nodes, const char * merged) {

	map<int, int> * parent_taxid = new map<int, int>;
	string line;

	ifstream myfile (nodes);
	if (myfile.is_open())
	{
		while (! myfile.eof() )	{
			getline (myfile,line, '\n');


			if (line.length() == 0) {
				continue;
			}


			int l_start = line.find_first_of ( '\t', 0);

			string tax_id = line.substr(0, l_start);
			int tax_id_int = str2int(tax_id);

			l_start += 3; // tab,pipe,tab

			int l_end = line.find_first_of ( '\t', l_start);

			string parent = line.substr(l_start, l_end - l_start);
			int parent_int = str2int(parent);

			//cout << "p: \"" << parent << "\"" << endl;
			(*parent_taxid)[tax_id_int]=parent_int;


			//int l_end = line.find_first_of ( '\t', l_start);



		}

		myfile.close();
	} else {
		cerr << "Error: Unable to open file " << nodes << endl;
		exit(1);
	}



	ifstream myfile2 (merged);
	if (myfile2.is_open())
	{
		while (! myfile.eof() )	{
			getline (myfile,line, '\n');


			if (line.length() == 0) {
				continue;
			}


			int l_start = line.find_first_of ( '\t', 0);

			string old_id = line.substr(0, l_start);
			int old_id_int = str2int(old_id);

			l_start += 3; // tab,pipe,tab

			int l_end = line.find_first_of ( '\t', l_start);

			string new_id = line.substr(l_start, l_end - l_start);
			int new_id_int = str2int(new_id);

			//cout << "p: \"" << parent << "\"" << endl;
			(*parent_taxid)[old_id_int]=(*parent_taxid)[new_id_int];


			//int l_end = line.find_first_of ( '\t', l_start);



		}

		myfile2.close();
	} else {
		cerr << "Error: Unable to open file " << merged << endl;
		exit(1);
	}



	return parent_taxid;
}

int main(int argc, char *argv[])
{




	//map<int, int> * parent_taxid = init_NCBI_NODES(nodes_dmp, merged_dmp);
	char * file;

	if (argc == 2) {
		file = argv[1];
		//cerr << file << endl;
		getPfamToNCBIFrompfamseq(file);

	} else {
		cerr << "ARGC != 2" << endl;
		exit(1);
	}






}
