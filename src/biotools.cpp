
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


#include "biotools.hpp"


void reverseComplementDNA(string& sequence) {
	int len = sequence.length();

	// A-T
	// G-C

	// reverse:
	int half = len/2;
	for (int i=0; i<half; i++) {
		swap(sequence[i], sequence[len-1-i]);
	}

	// complement:
	char c;
	for (int i=0; i<len; i++) {
		c = sequence[i];
		toupper(c);
		switch (c) {
			case 'A':
				c='T';
				break;
			case 'C':
				c='G';
				break;
			case 'G':
				c='C';
				break;
			case 'T':
				c='A';
				break;
			case 'U':
				c='C';
				break;
			default:
				//cerr << "error at " << i << " with " << sequence << endl;
				//exit(1);
				c='N';
		}
		sequence[i]=c;
	}


}

map<string, tax_id> * init_NCBI_NAMES(const char * names) {

	map<string, tax_id> * NameToNCBI = new map<string, tax_id>;


	ifstream myfile (names);
	if (myfile.is_open())
	{

		string line;

		while (! myfile.eof() )	{

			std::getline(myfile,line );


			if (line.length() == 0) {
				continue;
			}



			size_t pipe1 = line.find_first_of ( '|', 0);
			if (pipe1==string::npos) {
				cerr << "Error: init_NCBI_NAMES" << endl;
				exit(EXIT_FAILURE);
			}

			size_t pipe2 = line.find_first_of ( '|', pipe1+1);
			if (pipe2==string::npos) {
				cerr << "Error: init_NCBI_NAMES" << endl;
				exit(EXIT_FAILURE);
			}

			size_t pipe3 = line.find_first_of ( '|', pipe2+1);
			if (pipe3==string::npos) {
				cerr << "Error: init_NCBI_NAMES" << endl;
				exit(EXIT_FAILURE);
			}

			size_t pipe4 = line.find_first_of ( '|', pipe3+1);
			if (pipe4==string::npos) {
				cerr << "Error: init_NCBI_NAMES" << endl;
				exit(EXIT_FAILURE);
			}


			//if ((pipe4-pipe3) >= 17 && line.substr(pipe3, 17).compare("|\tscientific name")==0 ){

				//cout << line << endl;
				tax_id tax_id = str2int(line);
				//cout << tax_id << endl;

				int name_start = pipe1+2;
				int name_end = pipe2-1;

				string name = line.substr(name_start, name_end-name_start);
				stringToLower(name);
				(*NameToNCBI)[name]=tax_id;
				//cout << "-" << name << "-"<<endl;
			//}

			//string tax_id = line.substr(0, l_start);



		}

		myfile.close();
	}

	else {
		cerr << "Error: Unable to open file " << names << endl;
		exit(1);
	}

	return NameToNCBI;
}


tax_id taxonname2taxid(string taxon_string, map<string, tax_id> * NameToNCBI) {

	assert(NameToNCBI);

	map<string, tax_id>::iterator map_it;
	tax_id taxid;

	stringToLower(taxon_string);

	while (taxon_string.length() > 3) {


		map_it = NameToNCBI->find(taxon_string);

		taxid = 0;
		if (map_it != NameToNCBI->end()) {
			taxid = (*map_it).second;
			return taxid;
			//cout << "tax: " <<  taxid << endl;
		} else {

			size_t last_space = taxon_string.find_last_of(' ');
			if (last_space != string::npos) {
				//cout << "taxon_string: (" << taxon_string << ")"<<endl;
				taxon_string = taxon_string.substr(0, last_space);
				//cout << "taxon_string: (" << taxon_string << ")"<<endl;
				//exit(0);
			} else {
				return 0;
			}

			//cout << line<< endl;
			//cout << "(" << taxon_string << ")" << endl;
			//cout << "tax: not found !!!! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
			//unknown_count++;
			//exit(1);
		}

	}


	return 0;

}


void init_NCBI_NODES(const char * nodes, const char * merged, tax_id *& parent_taxid,  rank_t *& ranks, map<string, rank_t> *& rank_to_id, int NCBI_MAX) {




	parent_taxid = new tax_id[NCBI_MAX];
	for (int i = 0 ; i < NCBI_MAX; i++) parent_taxid[i] = 0;

	ranks = new rank_t[NCBI_MAX];
	for (int i = 0 ; i < NCBI_MAX; i++) ranks[i] = 0;

	rank_to_id = new map<string, rank_t>;


	//(*rank_to_id)["no rank"] = 0;
	rank_t last_used_rankid = 0;


	(*rank_to_id)["no rank"] 		= last_used_rankid++; // "0"
	(*rank_to_id)["superkingdom"]	= last_used_rankid++; // "1" and so on...
	(*rank_to_id)["kingdom"] 		= last_used_rankid++;
	(*rank_to_id)["subkingdom"] 	= last_used_rankid++;
	(*rank_to_id)["superphylum"] 	= last_used_rankid++;
	(*rank_to_id)["phylum"] 		= last_used_rankid++; // 5
	(*rank_to_id)["subphylum"] 		= last_used_rankid++;
	(*rank_to_id)["superclass"] 	= last_used_rankid++;
	(*rank_to_id)["class"] 			= last_used_rankid++;
	(*rank_to_id)["subclass"] 		= last_used_rankid++;
	(*rank_to_id)["infraclass"] 	= last_used_rankid++; //10
	(*rank_to_id)["cohort"] 		= last_used_rankid++;
	(*rank_to_id)["subcohort"] 		= last_used_rankid++;
	(*rank_to_id)["superorder"] 	= last_used_rankid++;
	(*rank_to_id)["order"]			= last_used_rankid++;
	(*rank_to_id)["suborder"] 		= last_used_rankid++; //15
	(*rank_to_id)["infraorder"] 	= last_used_rankid++;
	(*rank_to_id)["parvorder"] 		= last_used_rankid++;
	(*rank_to_id)["superfamily"] 	= last_used_rankid++;
	(*rank_to_id)["family"] 		= last_used_rankid++;
	(*rank_to_id)["subfamily"] 		= last_used_rankid++; //20
	(*rank_to_id)["tribe"] 			= last_used_rankid++;
	(*rank_to_id)["subtribe"] 		= last_used_rankid++;
	(*rank_to_id)["genus"] 			= last_used_rankid++;
	(*rank_to_id)["subgenus"] 		= last_used_rankid++;
	(*rank_to_id)["species group"] 	= last_used_rankid++; //25
	(*rank_to_id)["species subgroup"] = last_used_rankid++;
	(*rank_to_id)["species"] 		= last_used_rankid++;
	(*rank_to_id)["subspecies"] 	= last_used_rankid++;
	(*rank_to_id)["varietas"]		= last_used_rankid++;
	(*rank_to_id)["forma"] 			= last_used_rankid++; //30



	ifstream myfile (nodes);
	if (myfile.is_open())
	{

		string line;

		while (! myfile.eof() )	{

			std::getline(myfile,line );


			if (line.length() == 0) {
				continue;
			}

		//cout << line << endl;

			// get tax_id
			size_t l_start = line.find_first_of ( '\t', 0);
			if (l_start==string::npos) {
				cerr << "Error: init_NCBI_NODES" << endl;
				exit(EXIT_FAILURE);
				}
			string tax_id = line.substr(0, l_start);
			//cout << "XXXXXXXXXXXXXXXXXXXXx_TAXID: " << tax_id << endl;
			int tax_id_int = str2int(tax_id);



//if (tax_id_int > 40000) exit(0);
			// get parent
			l_start += 3; // tab,pipe,tab

			size_t l_end = line.find_first_of ( '\t', l_start);
			if (l_end==string::npos) {
				cerr << "Error: init_NCBI_NODES" << endl;
				exit(EXIT_FAILURE);
				}

			string parent = line.substr(l_start, l_end - l_start);
			int parent_int = str2int(parent);

			//cout << "p: \"" << parent << "\"" << endl;

			// get rank:
			l_start = l_end+3;
			l_end = line.find_first_of ( '\t', l_start);
			if (l_end==string::npos) {
				cerr << "Error: init_NCBI_NODES" << endl;
				exit(EXIT_FAILURE);
				}

			string rank_string = line.substr(l_start, l_end - l_start);
			//cout << "(" << rank_string << ")" << endl;

			map<string, rank_t>::iterator map_it;
			map_it =  rank_to_id->find(rank_string);

			rank_t rankid;
			if (map_it == rank_to_id->end()) {
					cerr << "not found: "<< rank_string << endl;

					exit(EXIT_FAILURE);

					last_used_rankid++;
					(*rank_to_id)[rank_string] = last_used_rankid;
					rankid=last_used_rankid;

			} else {
				//cout << "yes, found." << endl;
				rankid = (*map_it).second;

			}
//exit(0);
			if (tax_id_int > NCBI_MAX) {
					cerr << "sorry, NCBI_MAX is too small..." << endl;
					exit(EXIT_FAILURE);
			}

			if (parent_int > NCBI_MAX) {
					cerr << "sorry, NCBI_MAX is too small..." << endl;
					exit(EXIT_FAILURE);
			}

			parent_taxid[tax_id_int]=parent_int;
			ranks[tax_id_int]=rankid;
		}

		myfile.close();
	}

	else {
		cerr << "Error: Unable to open file " << nodes << endl;
		exit(1);
	}



	ifstream myfile2 (merged);
	if (myfile2.is_open())
	{
		string line;
		while (! myfile2.eof() )	{
			getline (myfile2,line, '\n');

			if (line.length() == 0) {
				continue;
			}


			size_t l_start = line.find_first_of ( '\t', 0);
			if (l_start==string::npos) {
				cerr << "Error: init_NCBI_NODES" << endl;
				exit(EXIT_FAILURE);
			}
			string old_id = line.substr(0, l_start);
			int old_id_int = str2int(old_id);

			l_start += 3; // tab,pipe,tab

			size_t l_end = line.find_first_of ( '\t', l_start);
			if (l_end==string::npos) {
				cerr << "Error: init_NCBI_NODES" << endl;
				exit(EXIT_FAILURE);
			}
			string new_id = line.substr(l_start, l_end - l_start);
			int new_id_int = str2int(new_id);

			if (old_id_int > NCBI_MAX) {
					cerr << "sorry, NCBI_MAX is too small..." << endl;
					exit(EXIT_FAILURE);
			}

			if (new_id_int > NCBI_MAX) {
					cerr << "sorry, NCBI_MAX is too small..." << endl;
					exit(EXIT_FAILURE);
			}

			parent_taxid[old_id_int]=parent_taxid[new_id_int];
			ranks[old_id_int]=ranks[new_id_int];

		}

		myfile2.close();
	}

	else {
		cerr << "Error: Unable to open file " << merged << endl;
		exit(1);
	}

	return;
}




map<tax_id, string> * init_taxid_to_name(const char * names) {

	map<tax_id, string> * taxid_to_name = new map<tax_id, string>;


	ifstream myfile (names);
	if (myfile.is_open())
	{

		string line;

		while (! myfile.eof() )	{

			std::getline(myfile,line );


			if (line.length() == 0) {
				continue;
			}



			size_t pipe1 = line.find_first_of ( '|', 0);
			if (pipe1==string::npos) {
				cerr << "Error: init_NCBI_NAMES" << endl;
				exit(EXIT_FAILURE);
			}

			size_t pipe2 = line.find_first_of ( '|', pipe1+1);
			if (pipe2==string::npos) {
				cerr << "Error: init_NCBI_NAMES" << endl;
				exit(EXIT_FAILURE);
			}

			size_t pipe3 = line.find_first_of ( '|', pipe2+1);
			if (pipe3==string::npos) {
				cerr << "Error: init_NCBI_NAMES" << endl;
				exit(EXIT_FAILURE);
			}

			size_t pipe4 = line.find_first_of ( '|', pipe3+1);
			if (pipe4==string::npos) {
				cerr << "Error: init_NCBI_NAMES" << endl;
				exit(EXIT_FAILURE);
			}


			if ((pipe4-pipe3) >= 17 && line.substr(pipe3, 17).compare("|\tscientific name")==0 ){

				//cout << line << endl;
				tax_id tax_id = str2int(line);
				//cout << tax_id << endl;

				int name_start = pipe1+2;
				int name_end = pipe2-1;

				string name = line.substr(name_start, name_end-name_start);
				stringToLower(name);
				(*taxid_to_name)[tax_id]=name;
				//cout << "-" << name << "-"<<endl;
			}

			//string tax_id = line.substr(0, l_start);



		}

		myfile.close();
	}

	else {
		cerr << "Error: Unable to open file " << names << endl;
		exit(1);
	}

	return taxid_to_name;
}


int mapRankToSimpleRank(rank_t rank){

	if ((int)rank >= 27) {
		// species
		return 0;
	}
	if ((int)rank >= 23) {
		// genus
		return 1;
	}
	if ((int)rank >= 19) {
		// family
		return 2;
	}
	if ((int)rank >= 14) {
		// order
		return 3;
	}
	if ((int)rank >= 8) {
		// class
		return 4;
	}
	if ((int)rank >= 5) {
		// phylum
		return 5;
	}

	if ((int)rank >= 1) {
		// superkingdom
		return 6;
	}

	// unknown, other
	return 7;

}

// if buffer is a null pointer, a new int array is created
int getTaxa(tax_id ncbi_tax_id, tax_id ** buffer_location, tax_id * parent_taxid, int NCBI_MAX) {
	tax_id * taxa;

	assert(parent_taxid);

	if (buffer_location == 0) {
		cerr << "Error: buffer location must be a valid pointer"<< endl;
		exit(EXIT_FAILURE);
	}

	if (*buffer_location == 0) {
		*buffer_location = new tax_id[256];
	}

	taxa = *buffer_location;

	if (ncbi_tax_id > NCBI_MAX) {
			cerr << "Error (getTaxa): ncbi_tax_id is too big..." << endl;
			exit(1);
	}

	tax_id temp_tax_id = ncbi_tax_id;

	taxa[0] = ncbi_tax_id;
	int last_element_position = 1;
//cout << endl;
	while (temp_tax_id > 1) {
//	cout << "reverse: " << (int)temp_tax_id << endl;

		temp_tax_id = parent_taxid[temp_tax_id];

		taxa[last_element_position] = temp_tax_id;
		last_element_position++;
	}
	taxa[last_element_position] = -1;
	last_element_position--;

	// wrong order, reverse:
	int half = (last_element_position+1)/2;

//cout << endl;

	for (int i = 0; i < half; i++) {
		tax_id tempx;
		tempx = taxa[last_element_position-i];
		taxa[last_element_position-i] = taxa[i];
		taxa[i] = tempx;

		//cerr << i << " " << (int)tempx << endl;
	}

	return last_element_position;

}

int get_LCP_Length(tax_id * a, tax_id * b) {
		int i = 0;
		while (1) {
				if (a[i] < 0 || b[i] < 0) {
						return i;
				}
				if (a[i] != b[i]) {
					return i;

				}
				i++;
		}

}

rank_t getLowestRank(tax_id taxid, tax_id * parent_taxid, rank_t * ranks, int NCBI_MAX){


	while (taxid > 1) {

		if (taxid > NCBI_MAX) {
				cerr << "Error (getLowestRank): taxid is too big" << endl;
				exit(1);
		}

		if (ranks[taxid] != 0) {
			return ranks[taxid];
		}

		taxid = parent_taxid[taxid];
		//last_element_position++;
	}
	return 0;
}


tax_id getTaxonAtRank(tax_id taxid, rank_t rank, tax_id * parent_taxid, rank_t * ranks, int NCBI_MAX){

	tax_id last_taxid = taxid;

	//cerr << endl <<  "start: " << (int) rank<< endl;

	while (taxid > 1) {

		rank_t current_rank = ranks[taxid];

		//cerr << "taxid: " << (int) taxid << endl;
		//cerr << "rank: " << (int) current_rank << endl;

		if (taxid > NCBI_MAX) {
				cerr << "Error (getLowestRank): taxid is too big" << endl;
				exit(1);
		}

		if (current_rank == rank) { // try to get correct rank
//			cerr << "A" << endl;
			return taxid;
		} else if ((current_rank < rank) && (current_rank != (rank_t)0)) { // if not possible, take this one
//			cerr << "B" << endl;
			return last_taxid;
		}
//cerr << "C" << endl;

		if (current_rank != (rank_t)0) {
				last_taxid = taxid;
		}

		taxid = parent_taxid[taxid];
		//last_element_position++;
	}
	return 0;
}



tax_id getLowestCommonTaxID(tax_id * taxonomy_A, tax_id * taxonomy_B, int NCBI_MAX) {

	int i = 0;
	while ((taxonomy_A[i] == taxonomy_B[i]) && (taxonomy_A[i] != -1)) {
		i++;
	}
	i--;

	if (i < 0) {
			return 0;
	}

	if (taxonomy_A[i] > NCBI_MAX ) {
		cerr << "Error (getLowestCommonTaxID): taxonomy_A[i] is too big" << endl;
		cerr << "i: " << i << endl;
		exit(1);
	}

	return taxonomy_A[i];

}


int reduceTaxonByRank(tax_id * buffer_taxid, rank_t * ranks, rank_t assignment_rank_id, int NCBI_MAX) {

		//int found_rank = 0;
		int last_valid_rank=0;
		int last_valid_node=-1;
		tax_id ncbi_tax_id;
		int node_pos = 1;


		if (buffer_taxid[0] == -1) {
				return -1;
		}

		if (buffer_taxid[1] > NCBI_MAX) {
			cerr << "buffer_taxid[1] > NCBI_MAX" << endl;
			exit(1);
		}

		while ((ncbi_tax_id = buffer_taxid[node_pos]) != -1) { // walking from superkingdom to species/strain....

			if (ncbi_tax_id > NCBI_MAX) {
				cerr << "ncbi_tax_id > NCBI_MAX" << endl;
				exit(1);
			}

			rank_t current_rank = ranks[ncbi_tax_id];
			//cout << "ncbi_tax_id: " << ncbi_tax_id << endl;
			//cout << "rank: " << (int) current_rank << endl;

			if (current_rank <= assignment_rank_id && current_rank > last_valid_rank)  {
				last_valid_rank = current_rank;
				last_valid_node=node_pos;
				//found_rank=1;
			}

			if (assignment_rank_id == current_rank)  {
				//cout << "GOT IT! " << (int) current_rank << endl;
				//found_rank=1;
				break;
			}
			node_pos++;
		}

		buffer_taxid[last_valid_node+1] = -1;


		return last_valid_node;

}


string getTaxonStringByTaxId(tax_id ncbi_tax_id, map<tax_id, string> ** taxid_to_name, tax_id * parent_taxid, string * names_dmp) {

	if (*taxid_to_name == 0) {
		*taxid_to_name = init_taxid_to_name(names_dmp->c_str());
	}

	string taxon = "";
	string name;

	if (ncbi_tax_id <= 1 ) {
			taxon = "unknown";
			return taxon;
	}


	tax_id temp_tax_id = ncbi_tax_id;
	map<tax_id, string>::iterator it;
	while (temp_tax_id > 1) {

		it = (*taxid_to_name)->find(temp_tax_id);

		if (it == (*taxid_to_name)->end()) {
				cerr << "error mapping taxid to name: " << (int)temp_tax_id << endl;
				name = "unknown";
		} else {
				name = it->second;
		}

		if (temp_tax_id == ncbi_tax_id) {
			taxon = name ;
		} else {
			taxon = name + "!" + taxon;
		}

		//cout << "found: " << it->second << endl;

		temp_tax_id = parent_taxid[temp_tax_id];

	}
	return taxon;
}


set<tax_id> * init_filtered_taxa(vector<tax_id> * filtered_species_vec, rank_t filter_rank , map<string, rank_t> * rank_to_id, tax_id * parent_taxid, rank_t * ranks, int NCBI_MAX){


	//size_t size_of_array = (sizeof filtered_species)/(sizeof filtered_species[0]);



	set<tax_id> * filtered_taxa = new set<tax_id>;

	vector<tax_id>::iterator vec_it;

	//for (int i = 0; i < size_of_array; i++) {
	for(vec_it=filtered_species_vec->begin(); vec_it != filtered_species_vec->end(); vec_it++) {
		tax_id species = *vec_it;
		tax_id filtered_taxon = getTaxonAtRank(species, filter_rank, parent_taxid, ranks, NCBI_MAX);
	//cerr << "fil:" << (int) filtered_taxon << endl;
		if (filtered_taxon == 0) {
				cerr << "error: (filtering) rank for species " << (int) species << " not found.." << endl;
				exit(1);
		}
		//cout << (int) species << " ___ " << (int) filtered_taxon << endl;
		filtered_taxa->insert(filtered_taxon);
	}

//exit(0);
	return filtered_taxa;
}


bool isFilteredTaxon(rank_t filter_rank, tax_id taxid, set<tax_id> ** filtered_taxa, map<string, rank_t> * rank_to_id, tax_id * parent_taxid, rank_t * ranks, int NCBI_MAX){


	//rank_t filter_rank = (*rank_to_id)["order"];
	//rank_t filter_rank = (*rank_to_id)["species"];

//cerr << "filter_rank: " << (int) filter_rank << endl;

	if (*filtered_taxa == 0) {
		cerr << "error: *filtered_taxa == 0" << endl;
		exit(1);
//		*filtered_taxa = init_filtered_taxa(filter_rank , rank_to_id, parent_taxid, ranks, NCBI_MAX);
	}




	while (taxid > 1) {

		if (taxid > NCBI_MAX) {
				cerr << "Error (getLowestRank): taxid is too big" << endl;
				exit(1);
		}

		if (ranks[taxid] >= filter_rank) {

			set<tax_id>::iterator it;

			it=(*filtered_taxa)->find(taxid);

			if (it != (*filtered_taxa)->end() ) {
					return true;
			}

		} else {
			if (ranks[taxid] != 0) {
				return false;
			}
		}

		taxid = parent_taxid[taxid];
		//last_element_position++;
	}


	// rank not found... assume "false".


	return false;
}

FASTA_Parser::~FASTA_Parser(){
	delete rna_stream;
	delete fpstream;
}


FASTA_Parser::FASTA_Parser(string command, bool getOnlyIdentfier) {

	this->start(command, getOnlyIdentfier);

}

FASTA_Parser::FASTA_Parser(string file, bool getOnlyIdentfier, string * zcat_bin) {

	string zcat_command = *zcat_bin;
	zcat_command.append(" -f ");
	zcat_command.append(file);

	this->start(zcat_command, getOnlyIdentfier);

}

void FASTA_Parser::start(string command, bool getOnlyIdentfier) {

	this->getOnlyIdentfier=getOnlyIdentfier;

	expected_sequence_length=0;


	fp = popen(command.c_str(), "r");
	if (fp == NULL) {
		cerr << "error: " << command << endl;
		exit(EXIT_FAILURE);
	}


	namespace io = boost::iostreams;

	io::file_descriptor_flags flags= io::close_handle;



	fpstream = new io::stream_buffer<io::file_descriptor_source>(fileno(fp), flags);

	rna_stream = new std::istream(fpstream);


	string line;


	while (1) {

		std::getline(*rna_stream, line);
		if (rna_stream->eof()) {
				return;
		}

		if (line.length() == 0) {
			continue;
		}

		if (line[0] == '#' ) {
			continue;
		}
		
		if (line[line.size() - 1] == '\r') {
			line.resize(line.size() - 1);
		}
		
		break;
	}

	if (line.length() < 2) {
			cerr << "FASTA_Parser: error parsing" << endl;
	}
	if (line[0] != '>') {
			cerr << "FASTA_Parser: error parsing" << endl;
	}
	if (getOnlyIdentfier) {
		this->description_line = getFASTADescription(line);
	} else {
		this->description_line = line;
	}
	have_read_descr = false;

}


string getFASTADescription (string & line) {
	size_t start = line.find_first_not_of(" \t", 1);

	if (start == string::npos)  {
			cerr << "error: no char found" << endl;
			exit(1);
	}

	size_t end = line.find_first_of(" \t", start);

	string identifier;

	if (end == string::npos)  {
		identifier = line.substr(start);
	} else {
		identifier = line.substr(start, end-start);
	}
	return identifier;
}

bool FASTA_Parser::getNextDescriptionLine(string& descr) {
	//cerr << "getNextDescriptionLine" << endl;
	if ( *rna_stream == false ) {
		return false;
	}

	// jump over sequence to next description line
	if (have_read_descr == true) {
		bool found_descr=false;
		string line;
		while (*rna_stream) {

			std::getline(*rna_stream, line);
			//cerr << "x:" << line << endl;

			if (line.length() == 0) {
				continue;
			}

			if (line[0] == '#' ) {
				continue;
			}
			
			if (line[line.size() - 1] == '\r') {
				line.resize(line.size() - 1);
			}

			if (line[0] == '>' ) {
				found_descr=true;
				//exit(1);
				this->description_line = getFASTADescription(line);
				if (getOnlyIdentfier) {
					this->description_line = getFASTADescription(line);
				} else {
					this->description_line = line;
				}
				break;
			}

		}
		if (! found_descr) {
			return false;
		}


	}

	have_read_descr = true;
	descr = this->description_line;
	return true;
}



string * FASTA_Parser::getSequence() {
	string line;



	string * rna_sequence = new string();

	if (expected_sequence_length > 0) {
		rna_sequence->reserve(expected_sequence_length);
	}

//cerr << "getSequence" << endl;

	//cout <<"line: "<< line << endl;
	while (*rna_stream) {
		std::getline(*rna_stream, line);
//cerr << "y:" << line << endl;
		if (line.length() == 0) {
			continue;
		}

		if (line[0] == '#' ) {
			continue;
		}

		if (line[line.size() - 1] == '\r') {
			line.resize(line.size() - 1);
		}
		
		if (line[0] == '>' ) {

				this->description_line = getFASTADescription(line);
				if (getOnlyIdentfier) {
					this->description_line = getFASTADescription(line);
				} else {
					this->description_line = line;
				}
				have_read_descr = false;
				rna_sequence->reserve(rna_sequence->length());
				return rna_sequence;
		} else {
			int lastpos = line.length()-1;
			if (line[lastpos] == '\n') {
				line.erase(lastpos); // remove trailing newline
			}
			//rna_sequence->append("@");
			rna_sequence->append(line);
		}
	}

	//file_closed = true;
	if (fp != 0) {
		pclose(fp);
		fp = 0;
	}
	rna_sequence->reserve(rna_sequence->length());
	return rna_sequence;

}

char Sequence_with_taxid::operator[](size_t i){
	size_t real_pos = i - offset;
	if (real_pos < 0) {
		return '?'; // this are leading and trailing gaps
	}
	if (real_pos >= string::length()) {
		return '?';
	}


	return this->at(real_pos);

}


string Sequence_with_taxid::substr(size_t pos, size_t len) {
	size_t real_pos = pos - offset;
	if (real_pos < 0) {
		cerr << "real_pos < 0" << endl;
		exit(1);
	}
	if (real_pos >= string::length()) {
		string fail_str = "substr fail";
		return fail_str;
		//cerr << "real_pos >= this->length(); wrong alignment area" << endl;
		//exit(1);
	}

	return string::substr(real_pos, len);
}

size_t Sequence_with_taxid::length() const{

	return string::length() + this->offset;
}

void Sequence_with_taxid::setAlignmentSequence(string& seq) {
	size_t start = seq.find_first_not_of('-');
	if (start == string::npos) {
		cerr << "error: find_first_not_of" << endl;
		exit(1);
	}

	size_t end = seq.find_last_not_of('-');
	if (end == string::npos) {
		cerr << "error: find_last_not_of" << endl;
		exit(1);
	}

	if (end < start) {
		cerr << "error: end < start" << endl;
		exit(1);
	}

	this->assign(seq.substr(start, end-start+1));
	this->offset=start;
}

Sequence::Sequence(string description, string data_sequence){
                this->description = description;
                this->data_sequence = data_sequence;

}



