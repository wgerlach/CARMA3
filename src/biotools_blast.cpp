
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

#include "biotools_blast.hpp"

Match_hit * getFirstUnmaskedHit(vector<Match_hit * > * hits) {

	assert(hits);

	vector<Match_hit * >::iterator hit_it;
//cerr<< "x --" << endl;
	for (hit_it = hits->begin(); hit_it != hits->end();  ++hit_it) {
		Match_hit * match_hit = *hit_it;
//cerr<< "x " << match_hit->ncbi_tax_id << endl;
//cerr<< "x " << match_hit->masked << endl;
		if (match_hit->tail) {
			break;
		}

		if (! match_hit->masked) {
			return match_hit;
		}

	}

	return 0;
}

string * extractSequenceID(string& complete_line) {
	// example: >ref|YP_026187.1| predicted dehydro....

	size_t id_start = complete_line.find_first_not_of(" >");
	if (id_start==string::npos) {
		cerr << "error extractSequenceID" << endl;
		exit(EXIT_FAILURE);
	}

	size_t id_end = complete_line.find_first_of(" ", id_start);
	if (id_end==string::npos) {
		cerr << "error extractSequenceID" << endl;
		exit(EXIT_FAILURE);
	}

	string * id = new string(complete_line.substr(id_start, id_end-id_start));

	size_t pipe=0;

	while ( (pipe = id->find_first_of('|', pipe)) != string::npos ) {
		id->insert(pipe, 1, '\\');
		pipe+=2;
	}

	return id;

}


int extractGiFromBlastDB(string * blast_id, map<string, tax_id> * NameToNCBI, int match_type, string * blast_database, string * fastacmd_script) {


	FILE *fp;
	int status;
	char path[PATH_MAX];

	int gi=-1;

	//string fastacmd_command = "fastacmd -d /vol/biodb/asn1/nr -T T -p T -s ";

	string fastacmd_command = *fastacmd_script;
	fastacmd_command.append(" -d ");

	fastacmd_command.append(*blast_database);
	//fastacmd_command.append("/vol/biodb/asn1/nr");

	if (match_type == MATCH_TYPE_BLASTX || match_type == MATCH_TYPE_BLASTP) {
		fastacmd_command.append(" -p T ");
	} else if (match_type == MATCH_TYPE_BLASTN) {
		fastacmd_command.append(" -p F ");
	} else {
		cerr << "error: unknown blast type.."	 << endl;
		exit(EXIT_FAILURE);
	}

	fastacmd_command.append(" -T T -s ");

	fastacmd_command.append(*blast_id);
//cerr << "call: " << fastacmd_command << endl;
	if (PRINT_EXTERNAL_CALLS == 1) {
		cerr << "call: " << fastacmd_command << endl;
	}
	fp = popen(fastacmd_command.c_str(), "r");
	if (fp == NULL) {
			cerr << "error: " << fastacmd_command << endl;
			exit(EXIT_FAILURE);
	}

	while (fgets(path, PATH_MAX, fp) != NULL) {
		//printf("%s", path);
		//cout << "get: " << strcmp( path, "NCBI taxonomy id: ") << endl;

		if (strncmp( path, "NCBI sequence id: ", 18) == 0 )  {
			//cout <<"GOOD! "  << path << endl;

			string test = path+18+3; // contains a pipe-sign, but this is no problem for str2int
			// example: NCBI sequence id: "gi|67921655|ref|ZP_00515173.1|" ---> "67921655|ref|ZP_00515173.1|"
			//cout <<"_"  << test << "_" << endl;
			gi = str2int(test);

			//cout <<"_"  << taxid << "_" << endl;
			break;
		}
	}

	// in order to avoid "broken pipe":
	while (fgets(path, PATH_MAX, fp) != NULL) {}

	//hits_fasta_file_stream.flush();
	status = pclose(fp);
	if (status == -1) {
		/* Error reported by pclose() */
		cerr << "-1" << endl;
		exit(EXIT_FAILURE);
	}

	if (gi == -1) {
			cerr << "failed: " << fastacmd_command << endl;
			cerr << "no gi found.." << endl;
			return 0;
			exit(EXIT_FAILURE);
	}

	//delete id;
	//exit(0);
	return gi;

}


tax_id extractTaxonomyInformationById(string * id, int match_type, string * blast_database, string * fastacmd_script) {



	tax_id taxid = 0;
	FILE *fp;
	int status;
	char path[PATH_MAX];

	// put backslashes in front of pipes
	string formatted_id = *id;

	size_t pipe = formatted_id.find_first_of('|');
	while (pipe != string::npos) {
		formatted_id.insert(pipe, 1 ,'\\');
		pipe = formatted_id.find_first_of('|', pipe+2);
	}

	//cerr << "formatted_id: " << formatted_id << endl;

//cerr << "hallo" << endl;
//cerr << *blast_database << endl;
	//string fastacmd_command = "fastacmd -d /vol/biodb/asn1/nr -T T -p T -s ";
	string fastacmd_command = *fastacmd_script;

	fastacmd_command.append(" -d ");
	fastacmd_command.append(*blast_database);

	if (match_type == MATCH_TYPE_BLASTX || match_type == MATCH_TYPE_BLASTP) {
		fastacmd_command.append(" -p T ");
	} else if (match_type == MATCH_TYPE_BLASTN) {
		fastacmd_command.append(" -p F ");
	} else {
		cerr << "error: unknown blast type.."	 << endl;
		exit(EXIT_FAILURE);
	}

	fastacmd_command.append(" -T T -s ");
	fastacmd_command.append(formatted_id);
	fastacmd_command.append(" 2>/dev/null ");
//cerr << "call: " << fastacmd_command << endl;
	if (PRINT_EXTERNAL_CALLS == 1) {
		cerr << "call: " << fastacmd_command << endl;
	}
	fp = popen(fastacmd_command.c_str(), "r");
	if (fp == NULL) {
			cerr << "error (extractTaxonomyInformationById): "<< fastacmd_command << endl;
			return 0;
			exit(EXIT_FAILURE);
	}

	while (fgets(path, PATH_MAX, fp) != NULL) {
		//printf("%s", path);
		//cout << "get: " << strcmp( path, "NCBI taxonomy id: ") << endl;

		if (strncmp( path, "NCBI taxonomy id: ", 18) == 0 )  {
			//cout <<"GOOD! "  << path << endl;

			string test = path+18; // contains a linebreak, but this is no problem for str2int
			//cout <<"_"  << test << "_" << endl;
			taxid = (tax_id) str2int(test);

			//cout <<"_"  << taxid << "_" << endl;
			break;
		}
	}

	while (fgets(path, PATH_MAX, fp) != NULL) {}

	//hits_fasta_file_stream.flush();
	status = pclose(fp);
	if (status == -1) {
		/* Error reported by pclose() */
		cerr << "-1" << endl;
		exit(EXIT_FAILURE);
	}
//exit(1);
	return taxid;
}


tax_id extractTaxonomyInformation(string& complete_line, map<string, tax_id> * NameToNCBI, bool use_fastacmd, int match_type, string * blast_database, string * fastacmd_script) {


	//cerr << "cline: "<< complete_line << endl;

	//bool use_fastacmd = true;

	if (use_fastacmd) {

		//fastacmd -d /vol/biodb/asn1/nr -T T -p T -s
		tax_id taxid = 0;
		//cout << complete_line << endl;

		string * id = extractSequenceID(complete_line);

		taxid = extractTaxonomyInformationById(id, match_type, blast_database, fastacmd_script);


		delete id;
		//exit(0);
		if (taxid == 0) {


				taxid = extractTaxonomyInformation(complete_line, NameToNCBI, false, match_type, blast_database, fastacmd_script);

				if (taxid == 0) {
					cerr << "failed: " << complete_line << endl;
					//cerr << "failed: " << fastacmd_command << endl;
				}
		}


		return taxid;


	} else {

		string subject_name;

//		if ((SwissProtToNCBI != 0) && complete_line.substr(1,3).compare("sp|")==0 ) { // SwissProt
//
//			size_t name_start = complete_line.find_first_of('_', 1);
//			size_t name_end = complete_line.find_first_of(' ', 1);
//
//			if (name_start!=string::npos && name_end!=string::npos) {
//				name_start++;
//				subject_name=complete_line.substr(name_start, name_end-name_start);
//				stringToLower(subject_name);
//				//cerr << "-" << subject_name << "-"<<endl;
//
//
//				map<string, tax_id>::iterator map_it;
//				map_it = SwissProtToNCBI->find(subject_name);
//
//				if (map_it == SwissProtToNCBI->end()) {
//					if (PRINT_UNKNOWN_SWISSPROT == 1) {
//						cerr << complete_line << endl;
//						cerr << "(sp) not found: "<< subject_name << endl;
//					}
//					return 0;
//				} else {
//
//					return (*map_it).second;
//				}
//			}
//		} else
		if (complete_line.substr(1,4).compare("pir|")==0) {



			size_t name_start = complete_line.find(" - ");

			if (name_start == string::npos) {
				return 0;

			}
			name_start+=3;

			string subject_name = complete_line.substr(name_start);
			stringToLower(subject_name);



			map<string, tax_id>::iterator map_it;
			map_it = NameToNCBI->find(subject_name);

			if (map_it != NameToNCBI->end()) {
				return (*map_it).second;
			}


	//cerr << "\""<<subject_name << "\"" << endl;

				//cerr << complete_line << endl;
				//cerr << subject_name << endl;

				//int name_start = complete_line.find(" - ");


			size_t last_close_bracket = subject_name.find_last_of(')');
			size_t last_open_bracket = subject_name.find_last_of('(');

			// no ()-brackets:
			if (last_close_bracket == string::npos || last_open_bracket == string::npos) {

				size_t last_space = subject_name.find_last_of(' ');

				while (last_space != string::npos ) {
					subject_name.erase(last_space);

					map_it = NameToNCBI->find(subject_name);
					if (map_it != NameToNCBI->end()) {

						return (*map_it).second;
					}
					last_space = subject_name.find_last_of(' ');
				}
				return 0;
			}

			string bracket_content = subject_name.substr(last_open_bracket+1, last_close_bracket-last_open_bracket-1);

			map_it = NameToNCBI->find(bracket_content);
			if (map_it != NameToNCBI->end()) {
				return (*map_it).second;
			}

			subject_name.erase(last_close_bracket, 1);
			subject_name.erase(last_open_bracket, 1);

			size_t substrain = subject_name.find(", substrain");
			if (substrain != string::npos) {
				subject_name.erase(substrain, 11);
			}

			size_t strain = subject_name.find("strain");
			if (strain != string::npos) {
				subject_name.erase(strain, 6);
			}


			// remove "  "
			size_t spaces = subject_name.find("  ");
			while (spaces != string::npos) {
				subject_name.erase(spaces, 1);
				spaces = subject_name.find("  ");
			}

			map_it = NameToNCBI->find(subject_name);
			if (map_it != NameToNCBI->end()) {
				return (*map_it).second;
			}



			size_t last_space = subject_name.find_last_of(' ');

			while (last_space != string::npos ) {
				subject_name.erase(last_space);

				map_it = NameToNCBI->find(subject_name);
				if (map_it != NameToNCBI->end()) {

					return (*map_it).second;
				}
				last_space = subject_name.find_last_of(' ');
			}


			return 0;

		} else {

			// "[name]" for: ref, gb, dbj, emb, tpg, tpd,
			//int name_start = complete_line.find_last_of ( '[');
			size_t name_end = complete_line.find_last_of ( ']');
			map<string, tax_id>::iterator map_it;

			if (name_end==string::npos) {

				stringToLower(complete_line);
				size_t from = complete_line.find(" from ");
				if (from != string::npos) {
					subject_name=complete_line.substr(from+6);


					map_it = NameToNCBI->find(subject_name);
					if (map_it != NameToNCBI->end()) {
						return (*map_it).second;
					}
				}


				size_t first_space = complete_line.find_first_of(' ');
				if (first_space == string::npos ) {
					return 0;
				}

				subject_name = complete_line.substr(first_space+1);

				size_t last_space = subject_name.find_last_of(' ');

				// suitable for blastn !
				while (last_space != string::npos ) {
					subject_name.erase(last_space);
			//cerr << subject_name << endl;
					map_it = NameToNCBI->find(subject_name);
					if (map_it != NameToNCBI->end()) {

						return (*map_it).second;
					}
					last_space = subject_name.find_last_of(' ');
				}

				return 0;
			}

			int height = 1;
			int pos = name_end-1;

			while (height > 0 && pos > 0) {
				if (complete_line[pos] == ']') {
					height++;
				} else if (complete_line[pos] == '[') {
					height--;
				}
				pos--;
			}

			int name_start;
			if (height==0) {
				name_start = pos+2;
			} else {
				return 0;
			}


			try {
				subject_name=complete_line.substr(name_start, name_end-name_start);
			} catch (out_of_range& oor) {
				cerr << "Out of Range error: " << oor.what() << endl;
				cerr << complete_line << endl;
				exit(EXIT_FAILURE);
			}
			stringToLower(subject_name);



			map_it = NameToNCBI->find(subject_name);

			if (map_it == NameToNCBI->end()) {
				//cerr << "-" << complete_line << "-"<<endl;
				//cerr << "-" << subject_name << "-"<<endl;

				// try to cut at comma
				size_t comma = subject_name.find_first_of (',', 0);
				if (comma!=string::npos) {
					subject_name=subject_name.substr(0, comma);
					map_it = NameToNCBI->find(subject_name);
					//cerr<< complete_line << endl;
					//cerr<< subject_name << endl;
					//exit(0);
				}

			}

			if (map_it == NameToNCBI->end()) {
				if (PRINT_UNPARSABLE_BLAST_DESCRIPTIONS == 1) {
					cerr<< complete_line << endl;
					cerr<< "(other)not found: "<< subject_name << endl;
				}
				return 0;

			} else {


				//cout << "found: " << ncbi_tax_id << endl;
				return (*map_it).second;
			}


		}
	} // end else use_fastacmd

	return 0;
}


double extractBLASTBitScore(string& line){

	size_t startofbitscore = line.find_first_not_of ( ' ', 8);
	if (startofbitscore==string::npos) {
		cerr << "Error (extractBLASTBitScore): find_first_of " << endl;
		exit(EXIT_FAILURE);
	}

	size_t endofbitscore = line.find_first_of ( ' ', startofbitscore);
	if (endofbitscore==string::npos) {
		cerr << "Error (extractBLASTBitScore): find_first_of " << endl;
		exit(EXIT_FAILURE);
	}
	string bitscore_string=line.substr(startofbitscore, endofbitscore-startofbitscore);

	double bitscore_double = atof(bitscore_string.c_str());

	return bitscore_double;
}


double extractBLASTEvalue(string& line){

	size_t startofevalue = line.find_last_of (' ');

	if (startofevalue==string::npos) {
		cerr << "Error (extractBLASTEvalue): find_first_of " << endl;
		exit(EXIT_FAILURE);
	}
	startofevalue++;

	string evalue_string  = line.substr(startofevalue);

	double evalue_double = atof(evalue_string.c_str());

	return evalue_double;

}

size_t extractIdentities(string& line) {

	// "Identities = 22/34 (64%), Positives = 26/34 (76%)"

	size_t startofidentities = line.find_first_of('=');
	if (startofidentities==string::npos) {
		cerr << "Error (extractIdentities 1): find_first_of " << endl;
		cerr << line << endl;
		exit(EXIT_FAILURE);
	}
	startofidentities+=2;

	size_t identities_slash = line.find_first_of('/');
	if (identities_slash==string::npos) {
		cerr << "Error (extractIdentities 2): find_first_of " << endl;
		cerr << line << endl;
		exit(EXIT_FAILURE);
	}

	string identities_string = line.substr(startofidentities, identities_slash-startofidentities);
	size_t identities_int = atof(identities_string.c_str());


	return identities_int;
}



size_t extractAlignmentLength(string& line) {

	// "Identities = 22/34 (64%), Positives = 26/34 (76%)"

	size_t startofalignmentlength = line.find_first_of('/');
	if (startofalignmentlength==string::npos) {
		cerr << "Error (extractAlignmentLength): find_first_of " << endl;
		exit(EXIT_FAILURE);
	}
	startofalignmentlength+=1;

	size_t endofalignmentlength = line.find_first_of(' ', startofalignmentlength);
	if (endofalignmentlength==string::npos) {
		cerr << "Error (extractAlignmentLength): find_first_of " << endl;
		exit(EXIT_FAILURE);
	}

	string alignmentlength_string = line.substr(startofalignmentlength, endofalignmentlength-startofalignmentlength);
	size_t alignmentlength_int = atof(alignmentlength_string.c_str());


	return alignmentlength_int;
}


size_t extractPercentIdentities(string& line) {
	// "Identities = 22/34 (64%), Positives = 26/34 (76%)"
	// **** Percent Identities
	size_t identities_bracket_open = line.find_first_of('(');
	if (identities_bracket_open==string::npos) {
		cerr << "Error (extractPercentIdentities): find_first_of " << endl;
		exit(EXIT_FAILURE);
	}

	size_t identities_bracket_close = line.find_first_of(')', identities_bracket_open);
	if (identities_bracket_close==string::npos) {
		cerr << "Error (extractPercentIdentities): find_first_of " << endl;
		exit(EXIT_FAILURE);
	}

	string percent_intentity_string = line.substr(identities_bracket_open+1, identities_bracket_close-identities_bracket_open-2);
	//cout << percent_intentity_string << endl;
	size_t percent_intentity_int = atof(percent_intentity_string.c_str());
	//cout << percent_intentity_int << endl;

	return percent_intentity_int;
}

size_t extractPositives(string& line) {
	// "Identities = 22/34 (64%), Positives = 26/34 (76%)"
	// **** Positives
	size_t startofpositives = line.find("Positives = ", 0);
	if (startofpositives==string::npos) {
		cerr << "Error: find_last_of  1" << endl;
		exit(EXIT_FAILURE);
	}

	startofpositives+=12;

	size_t slash = line.find_first_of ('/', startofpositives);


	if (slash==string::npos) {
		cerr << "Error: find_last_of 2" << endl;
		exit(EXIT_FAILURE);
	}

	string positives_string = line.substr(startofpositives, slash-startofpositives);
//cout << "positives_string: " << positives_string << endl;
	size_t positives_int = atof(positives_string.c_str());

//cout << "line: " << line << endl;
//cout << "positives_int: " << positives_int << endl;
	return positives_int;
}


size_t extractGaps(string& line) {

	// " Identities = 32/59 (54%), Positives = 39/59 (66%), Gaps = 1/59 (1%)"

	// **** Positives
	size_t startofgaps = line.find("Gaps = ", 0);
	if (startofgaps==string::npos) {
		return 0;
	}
//cout << "startofgaps" << (int)startofgaps<< endl;
	startofgaps+=7;

	size_t slash = line.find_first_of ('/', startofgaps);


	if (slash==string::npos) {
		cerr << "Error: find_last_of 3" << endl;
		exit(EXIT_FAILURE);
	}

	string gaps_string = line.substr(startofgaps, slash-startofgaps);
	//cout << "gaps_string: -" << gaps_string << "-" <<endl;
	size_t gaps_int = atof(gaps_string.c_str());


//cout << "line: " << line << endl;

//cout <<"gaps_int : "<<gaps_int << endl;
//exit(0);
	return gaps_int;
}


// index blast-entry only if a Hit/HSP with some good score exists.
// using index_non_hits==true will "index" bad hits with "-1" as position
map<string, ifstream::pos_type > * indexBlastFile(const char * blast_file) {

	map<string, ifstream::pos_type > * index = new map<string, ifstream::pos_type >;
//int count = 0;
	bool already_added;
	string query;
	ifstream::pos_type pos=-1;
	ifstream::pos_type line_start=-1;

	ifstream myfile (blast_file);
	if (myfile.is_open())
	{
		string line;

		while (! myfile.eof() )	{
			line_start = myfile.tellg();
			std::getline(myfile,line );

			if (mycompare(line, 0,5, "BLAST")== 0) {
				pos = line_start;
			}
//cerr << line << endl;
			if (mycompare(line, 0,6, "Query=")== 0) {

				already_added=false;

				query = line.substr(7);

				index->insert(pair<string, ifstream::pos_type>(query, pos));
//cerr << query << " " << pos << endl;

//count++;
//if (count > 1000) {
//	return index;
//	}
			}

		}
//exit(1);
	} else {
		cerr << "Error: Unable to open file " << blast_file << endl;
		exit(1);
	}

	return index;
}



int extractFrame(string& line){
	// exmaple: Frame = -2

	size_t start_frame = line.find_last_of(' ');
	if (start_frame == string::npos) {
		cerr << "error extractFrame: " << line << endl;
		exit(EXIT_FAILURE);
	}

	string frame_string = line.substr(start_frame);

	int frame = str2int(frame_string);
	return frame;
}

int extractSequenceStartNumber(string& line){
	//example: Sbjct: 69 LLSFRLNDM TVIEGGPGGATWWD KVPSKFEG 99

	size_t start_number = line.find_first_of("0123456789");
	if (start_number == string::npos) {
		cerr << "error extractSequenceStartNumber 1:" << line << endl;
		exit(EXIT_FAILURE);
	}

	size_t end_number = line.find_first_not_of("0123456789", start_number);
	if (end_number == string::npos) {
		cerr << "error extractSequenceStartNumber 2:" << line << endl;
		exit(EXIT_FAILURE);
	}


	string number_string = line.substr(start_number, end_number-start_number+1);

	int position = str2int(number_string);


	return position;
}

int extractSequenceEndNumber(string& line){
	//example: Sbjct: 69 LLSFRLNDM TVIEGGPGGATWWD KVPSKFEG 99

	size_t end_number = line.find_last_of("0123456789");
	if (end_number == string::npos) {
		cerr << "errorA: " << line << endl;
		exit(EXIT_FAILURE);
	}

	size_t start_number = line.find_last_not_of("0123456789", end_number);
	if (start_number == string::npos) {
		cerr << "errorB: " << line << endl;
		exit(EXIT_FAILURE);
	}


	string number_string = line.substr(start_number+1, end_number-start_number+1);

	int position = str2int(number_string);


	return position;
}


string extractSequence(string& line){
//example: Sbjct: 69 LLSFRLNDM TVIEGGPGGATWWD KVPSKFEG 99

	size_t start_seq = line.find_first_of(':', 1);
	if (start_seq == string::npos) {
		cerr << "errorA: " << line << endl;
		exit(EXIT_FAILURE);
	}

	// search start of number
	start_seq = line.find_first_not_of(" 0123456789" , start_seq+1);
	if (start_seq == string::npos) {
		cerr << "errorB: " << line << endl;
		exit(EXIT_FAILURE);
	}

	size_t end_seq = line.find_last_not_of(" 0123456789");
	if (end_seq == string::npos) {
		cerr << "errorC: " << line << endl;
		exit(EXIT_FAILURE);
	}

	string sequence = line.substr(start_seq, end_seq-start_seq+1);

	size_t pos = sequence.find_first_of(" -/\\");
	while (pos != string::npos) {
		sequence.erase(pos, 1);
		pos = sequence.find_first_of(" -/\\", pos);
	}
	return sequence;
}

int getFrameShiftNum(string& line){
//example: Sbjct: 69 LLSFRLNDM TVIEGGPGGATWWD KVPSKFEG 99

	size_t start_seq = line.find_first_of(':', 1);
	if (start_seq == string::npos) {
		cerr << "errorA: " << line << endl;
		exit(EXIT_FAILURE);
	}

	// search start of number
	start_seq = line.find_first_not_of(" 0123456789" , start_seq+1);
	if (start_seq == string::npos) {
		cerr << "errorB: " << line << endl;
		exit(EXIT_FAILURE);
	}

	size_t end_seq = line.find_last_not_of(" 0123456789");
	if (end_seq == string::npos) {
		cerr << "errorC: " << line << endl;
		exit(EXIT_FAILURE);
	}

	string sequence = line.substr(start_seq, end_seq-start_seq+1);

	int frameShiftNum = 0;

	size_t pos = sequence.find_first_of("/\\");
	while (pos != string::npos) {
		//sequence.erase(pos, 1);
		frameShiftNum++;
		pos = sequence.find_first_of("/\\", pos+1);
	}
	return frameShiftNum;
}


string extractSequenceFromBlastDB(Match_hit * match_hit, Blast_HSP * blast_hsp, int match_type, int best_query_start, int best_query_end, string * blast_nr_database, string * blast_nt_database, char * path, int max_blast_description_length, string * fastacmd_script){



	int subject_start = blast_hsp->subject_start;
	int subject_end = blast_hsp->subject_end;

	if (match_type == MATCH_TYPE_BLASTN) {
			if (blast_hsp->frame < 0) {
				subject_end = blast_hsp->subject_start;
				subject_start = blast_hsp->subject_end;
			}
	}



	int query_start = blast_hsp->query_start;
	int query_end = blast_hsp->query_end;

	int query_length = abs(query_end-query_start);
	//cout << endl;
	//cout << "subject_start: " << subject_start << endl;
	//cout << "subject_end: " << subject_end << endl;

	//cout << "query_start: " << query_start << endl;
	//cout << "query_end: " << query_end << endl;

	// compute overlap:
	int overlap_start = (best_query_start>query_start)?best_query_start:query_start; // max
	int overlap_end = (best_query_end<query_end)?best_query_end:query_end; // min

	int overlap_size = abs(overlap_end-overlap_start);

//	cout << "overlap: " << overlap_size << endl;

	if ((double)overlap_size/(double)query_length < 0.25 ) {

		// TODO try next HSP !!!!!

		//cerr  <<  "q: " << matching_query->query << endl;
		//cerr  <<  "overlap_size: " << overlap_size << endl;
		//cerr  <<  "query_length: " << query_length << endl;
		//cerr  <<  "no enough overlap... but continue... " << endl;

		//cerr << *(blast_hsp->subject_sequence) << endl;
		//cerr << blast_hsp->subject_start << endl;

		//cerr << rm_command << endl;
		//exit(1);
	}
	//cout << blast_hsp->frame << endl;
//cout << subject_start << endl;
//cout << subject_end << endl;

	if (subject_start > subject_end) {
		cerr << "error: subject_start > subject_end"  << endl;
		cerr << "frame: " << blast_hsp->frame << endl;
		cerr << "subject_start:" << subject_start << endl;
		cerr << "subject_end:" << subject_end << endl;
		cerr << "blast_id:" << *(match_hit->blast_id) << endl;
		exit(1);
	}

	subject_start-= (1+5);
	if (subject_start < 0) {
		subject_start = 0;
	}

	subject_end+=5;

//subject_start=5;
//subject_end=2100;





	string fastacmd_command;

	if (match_hit->blast_id == 0) {
			cerr << "match_hit->blast_id == 0" << endl;
			exit(1);
	}

	//fastacmd_command = "fastacmd -d /vol/biodb/asn1/nr -p T -s ";
	fastacmd_command.append(*fastacmd_script);

	fastacmd_command.append(" -d ");

	if (match_type == MATCH_TYPE_BLASTX || match_type == MATCH_TYPE_BLASTP) {

		fastacmd_command.append(*blast_nr_database);
		fastacmd_command.append(" -p T -s ");
	} else if (match_type == MATCH_TYPE_BLASTN) {
		fastacmd_command.append(*blast_nt_database);
		fastacmd_command.append(" -p F -s ");
	} else {
			cerr << "unknown blast type" << endl;
			exit(1);
	}

	string formatted_id = *match_hit->blast_id;

//	size_t pipe = formatted_id.find_first_of('|');
//	while (pipe != string::npos) {
//		formatted_id.insert(pipe, 1 ,'\\');
//		pipe = formatted_id.find_first_of('|', pipe+2);
//	}

	fastacmd_command.append(formatted_id);




	//fastacmd_command.append(" 2>&1");
//cerr << fastacmd_command << endl;

	if (PRINT_EXTERNAL_CALLS == 1) {
		cerr << "call: " << fastacmd_command << endl;
	}
	FILE * fp = popen(fastacmd_command.c_str(), "r");
	if (fp == NULL) {
			cerr << "error (extractSequenceFromBlastDB): " << fastacmd_command << endl;
			exit(EXIT_FAILURE);
	}
//cerr << "call done... " << endl;
//cout << "f3" << endl;
	//long_sequence_buffer[0] = '\0';

	string extracted_sequence = "";
	int current_position = 0;
	int db_sequence_length = 0;
//cout  <<  "q: " << matching_query->query << endl;
	while (true) {
		char * c = fgets(path, max_blast_description_length, fp);

		if (c == NULL) {

			if (ferror(fp) != 0) {
				cerr << "error reading from stream" << endl;
			}
			if ( feof(fp) != 0) {
				//cerr << "EOF reached" << endl;
				break; // eof reached
			} else {
				//cerr << "no EOF continue" << endl;
				continue;

			}


		}

		if (path[0] != '>') {
			int line_length = strlen(path) - 1 ; // remove one for endline char
			db_sequence_length += line_length;
			//cout << "path: " << path << endl;
			//cout << "current_position: " << current_position << endl;
			//cout << "db_sequence_length: " << db_sequence_length << endl;

			if ( (subject_start >= current_position) && (subject_start < db_sequence_length)) {
				// start
				if (subject_end <= db_sequence_length) {
					extracted_sequence.append(path+(subject_start-current_position), (subject_end-subject_start));
					//cout  <<  "extr: " << extracted_sequence << endl;
					break;
				}
//cout << "b" << endl;
				extracted_sequence.append(path+(subject_start-current_position), line_length-(subject_start-current_position));

				//cout  <<  "prefix: " << extracted_sequence << endl;

			} else if ( (subject_end >= current_position) && (subject_end <= db_sequence_length) ) {
				extracted_sequence.append(path, (subject_end-current_position+1));
				//cout  <<  "suffix: -" << extracted_sequence << "-"<<endl;
				break;
			} else if ((subject_start < current_position) && (subject_end > db_sequence_length)) {
				//cout << "d" << endl;
				extracted_sequence.append(path, line_length);
				//cout  <<  "center: " << extracted_sequence << endl;
			}
			current_position = db_sequence_length+1;

		}
	}

	while (feof(fp) == 0 ) {
		char * c = fgets(path, max_blast_description_length, fp); // TODO need better solution for flush
	}

	int status = pclose(fp);
	if (status == -1) {
		/* Error reported by pclose() */
		cerr << "fastacmd exited with -1, continue..." << endl;
		//exit(EXIT_FAILURE);
	}


	if (extracted_sequence[0] == ']') {
			cerr << "error: ] " << endl;
			cerr << fastacmd_command << endl;
			exit(1);
	}

	return extracted_sequence;
}

string extractExactSubstringFromBlastDB(string * sequence_id, int start, int end,  int match_type, string * blast_database, string * fastacmd_script){

	bool reverse_complement = false;

	if (start > end) {
		reverse_complement=true;;
		swap(start, end);
	}


	string fastacmd_command;


	//fastacmd_command = "fastacmd -d /vol/biodb/asn1/nr -p T -s ";
	fastacmd_command.append(*fastacmd_script);

	fastacmd_command.append(" -d ");
	fastacmd_command.append(*blast_database);

	if (match_type == MATCH_TYPE_BLASTX || match_type == MATCH_TYPE_BLASTP) {
		fastacmd_command.append(" -p T ");
	} else if (match_type == MATCH_TYPE_BLASTN) {
		fastacmd_command.append(" -p F ");
	} else {
			cerr << "unknown blast type" << endl;
			exit(1);
	}

	fastacmd_command.append(" -L ");
	fastacmd_command.append(int2str(start));
	fastacmd_command.append(",");
	fastacmd_command.append(int2str(end));

	fastacmd_command.append(" -s ");
	string formatted_id = *sequence_id;

	size_t pipe = formatted_id.find_first_of('|');
	while (pipe != string::npos) {
		formatted_id.insert(pipe, 1 ,'\\');
		pipe = formatted_id.find_first_of('|', pipe+2);
	}

	fastacmd_command.append(formatted_id);




	//fastacmd_command.append(" 2>&1");
//cerr << fastacmd_command << endl;

	if (PRINT_EXTERNAL_CALLS == 1) {
		cerr << "call: " << fastacmd_command << endl;
	}



	namespace io = boost::iostreams;
	io::file_descriptor_flags flags= io::close_handle;

	FILE * fp = popen(fastacmd_command.c_str(), "r");
	if (fp == NULL) {
		cerr << "error(extractExactSubstringFromBlastDB): " << fastacmd_command << endl;
		exit(EXIT_FAILURE);
	}

	io::stream_buffer<io::file_descriptor_source> * fpstream = new io::stream_buffer<io::file_descriptor_source>(fileno(fp), flags);
	std::istream * fastacmd_stream = new std::istream(fpstream);

	int len = end-start+1;

	string line;
	string sequence="";
	sequence.reserve(len);

	int count =0;
	while (*fastacmd_stream) {

			std::getline(*fastacmd_stream, line);

			//cerr << "line: " << line << endl;
			if (line.length() == 0) {
				continue;
			}

			if (line[0] == '>') {
				continue;
			}

			//cerr << "line: " << line << endl;
			count++;



			sequence.append(line);

	}


	if (line[0] == ']') {
			cerr << "error]" << endl;
			cerr << fastacmd_command << endl;
			exit(1);
	}


	if (reverse_complement) {
		//cerr << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX: "<< sequence << endl;
		reverseComplementDNA(sequence);
		//cerr << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX: " << sequence << endl;
	}

	int status = pclose(fp);
	if (status == -1) {
		/* Error reported by pclose() */
		cerr << "-1" << endl;
		exit(EXIT_FAILURE);
	}
	return sequence;
}


Fastacmd_Class::Fastacmd_Class() {
	data = new map<string, pair<tax_id, string * > >();
	count_new = 0;
	count_old = 0;

}

Fastacmd_Class::~Fastacmd_Class() {

}

string Fastacmd_Class::fastacmd_sequence(string id, int start, int end, string * blast_database, string * fastacmd_script){

	if (id.compare(0, 3, "gi|") != 0) {
		cerr << "error: prefix wrong" << endl;
		exit(1);
	}

	size_t first_pipe = id.find_first_of('|');
	if (first_pipe == string::npos) {
		cerr << "error: first_pipe == string::npos" << endl;
		exit(1);
	}

	size_t second_pipe = id.find_first_of('|', first_pipe+1);
	if (second_pipe == string::npos) {
		cerr << "error: second_pipe == string::npos" << endl;
		exit(1);
	}

	string gi_identifier = id.substr(first_pipe+1, second_pipe-first_pipe-1);

	//cerr << id << endl;
	//cerr << gi_identifier << endl;
	//exit(1);

	string * database_sequence = 0;

	bool reverse_complement = false;

	if (start > end) {
		reverse_complement=true;
		swap(start, end);
	}


	map<string, pair<tax_id, string * > >::iterator data_it;


	data_it = this->data->find(gi_identifier);

	if ((count_new+count_old) % 100 == 0) {
		cerr << count_new << endl;
		cerr << count_old << endl;
		cerr << "example: " << gi_identifier << endl;
	}

	if (data_it != data->end()) {
		count_old++;
		//database_sequence
		pair<tax_id, string *> pp  = data_it->second;
		database_sequence = pp.second;
		//cerr << "found" << endl;
		//exit(1);
	} else {
		count_new++;




		string fastacmd_command;


		//fastacmd_command = "fastacmd -d /vol/biodb/asn1/nr -p T -s ";
		fastacmd_command.append(*fastacmd_script);

		fastacmd_command.append(" -d ");
		fastacmd_command.append(*blast_database);

		fastacmd_command.append(" -p T ");


		fastacmd_command.append(" -s ");
		string formatted_id = gi_identifier;

		size_t pipe = formatted_id.find_first_of('|');
		while (pipe != string::npos) {
			formatted_id.insert(pipe, 1 ,'\\');
			pipe = formatted_id.find_first_of('|', pipe+2);
		}

		fastacmd_command.append(formatted_id);




		//fastacmd_command.append(" 2>&1");
	//cerr << "actual: "<< fastacmd_command << endl;

		if (PRINT_EXTERNAL_CALLS == 1) {
			cerr << "call: " << fastacmd_command << endl;
		}



		namespace io = boost::iostreams;
		io::file_descriptor_flags flags= io::close_handle;

		FILE * fp = popen(fastacmd_command.c_str(), "r");
		if (fp == NULL) {
			cerr << "error(extractExactSubstringFromBlastDB): " << fastacmd_command << endl;
			exit(EXIT_FAILURE);
		}

		io::stream_buffer<io::file_descriptor_source> * fpstream = new io::stream_buffer<io::file_descriptor_source>(fileno(fp), flags);
		std::istream * fastacmd_stream = new std::istream(fpstream);

		int len = end-start+1;

		string line;
		string sequence="";
		sequence.reserve(len);

		int count =0;
		while (*fastacmd_stream) {

				std::getline(*fastacmd_stream, line);

				//cerr << "line: " << line << endl;
				if (line.length() == 0) {
					continue;
				}

				if (line[0] == '>') {
					continue;
				}

				//cerr << "line: " << line << endl;
				count++;



				sequence.append(line);

		}


		if (line[0] == ']') {
				cerr << "error]" << endl;
				cerr << fastacmd_command << endl;
				exit(1);
		}


		if (reverse_complement) {
			//cerr << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX: "<< sequence << endl;
			reverseComplementDNA(sequence);
			//cerr << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX: " << sequence << endl;
		}

		int status = pclose(fp);
		if (status == -1) {
			/* Error reported by pclose() */
			cerr << "-1" << endl;
			exit(EXIT_FAILURE);
		}

		//cerr << sequence << endl;
		database_sequence = new string(sequence);
		pair<tax_id, string * > pp = pair<tax_id, string * >(0, database_sequence);
		data->insert(pair<string, pair<tax_id, string * > >(gi_identifier,  pp));




	}


		string final = database_sequence->substr(start, end-start);
		if (reverse_complement) {
		reverseComplementDNA(final);

		}
		return final;

}

int getQueryScore(vector<string> * line_data, int score_match, int score_mismatch, int score_gapopen, int score_gapextension){
	// Fields: (0) Query id, (1)Subject id, (2)% identity, (3)alignment length, (4)mismatches, (5)gap openings, (6)q. start, (7)q. end, (8)s. start, (9)s. end, (10)e-value, (11)bit score
	double ali_identity = string2double((*line_data)[2]) / 100;

	//cerr << "ali_identity: " << ali_identity << endl;
	int ali_len = str2int((*line_data)[3]);
	//cerr << "ali_len: " << ali_len << endl;

	int ali_matches = (int)(ali_identity*(double)ali_len);
	//cerr << "computed ali_matches: " << ali_matches << endl;

	int ali_mismatches = str2int((*line_data)[4]);
	//cerr << "ali_mismatches: " << ali_mismatches << endl;

	int ali_gapopenings = str2int((*line_data)[5]);
	//cerr << "ali_gapopenings: " << ali_gapopenings << endl;

	int ali_gaps = ali_len-(ali_matches+ali_mismatches);
	//cerr << "computed ali_gaps: " << ali_gaps << endl;

	int query_score =  ali_matches*score_match + ali_mismatches*score_mismatch + ali_gapopenings*score_gapopen + ali_gaps*score_gapextension;
	//cerr << "query_score: " << query_score << endl;
	return query_score;
}


int runBLAST(string * blastall_script, string * blast_database, double evalue, string * input, string temp_output_file, char blast_type , bool use_gzip, string * gzip_bin){

	//cout << "hello"	 << endl;
	string evalue_str = double2str(evalue);

	//char blast_output_file_name [L_tmpnam];
	//char * tmpnam_res = tmpnam ( blast_output_file_name );

//cout << "blast_output_file_name: " << blast_output_file_name << endl;



//exit(0);

	string blast_command;

	blast_command.append(*blastall_script);

	if (blast_type == 'x') {
		blast_command.append(" -p blastx");
	} else if (blast_type == 'n' || blast_type == 'r') {
		blast_command.append(" -p blastn");
	} else if (blast_type == 'p') {
		blast_command.append(" -p blastp");
	} else {
		cerr << "error: unkown blast type: " << blast_type << endl;
		exit(1);
	}


	blast_command.append(" -d ");
	blast_command.append(*blast_database);
	blast_command.append(" -e ");
	blast_command.append(evalue_str);


// disable gaps:
//blast_command.append(" -g F ");


	if (blast_type == 'x') {
		blast_command.append(" -F \"m S\"");
		blast_command.append(" -w 15");
	} else if (blast_type == 'n'){
		blast_command.append(" -F \"m D\" -m9");
	} else if (blast_type == 'r') {
		blast_command.append(" -m9"); // -m 9 = tabular output
	} else if (blast_type == 'p') {
		blast_command.append(" -F \"m S\" -m9");
	} else {
		cerr << "error: unkown blast type: " << blast_type << endl;
		exit(1);
	}



	blast_command.append(" -i \"");
	blast_command.append(*input);
	blast_command.append("\"");

	if (use_gzip) {
		if (gzip_bin == 0) {
			cerr << "error: gzip_bin == 0" << endl;
			exit(1);
		}

		blast_command.append(" | ");
		blast_command.append(*gzip_bin);
		blast_command.append(" > ");
		blast_command.append(temp_output_file);
	} else {

		if (temp_output_file.length() > 0) {
			blast_command.append(" -o ");
			blast_command.append(temp_output_file);
			//blast_command.append(blast_output_file_name);
		}
	}


	if (PRINT_EXTERNAL_CALLS == 1) {
		cerr << "call: " << blast_command << endl;
	}
cerr << "call: " << blast_command << endl;
	int blast_ret = system(blast_command.c_str());

	if (blast_ret != 0) {

		cerr << "blast exited with error:" << blast_ret << endl;

			//if (FileExists(blast_output_file_name)) {
			//	remove(blast_output_file_name);
			//}

		return 1;
	}




	// **** ugly stuff to check if the blast file is ok...*****


	bool checkBLASToutput = true;

	if (checkBLASToutput && temp_output_file.length() > 0 && blast_type=='x' && (not use_gzip)) {
		// check last line of BLAST output
		ifstream myfile;
//cerr << "check blast output: " << temp_output_file << endl;
		myfile.open(temp_output_file.c_str(), ios::binary);
		//myfile.open(blast_output_file_name, ios::binary);

		char c=0;
		int i=0;
		//char lastline[256];
		string lastline;
		if (myfile.is_open()) {

			myfile.seekg(i, ios::end);
			int length = myfile.tellg();

			if (length < 5) {
					cerr << "error: BLAST file empty !" << endl;
					exit(1);
			}

			i = length -2;
			myfile.clear();

			while ((c != '\n') ||  (i > (length -3)) ) {

				if (i < 10) {
					cerr << "error: BLAST file corrupt ??" << endl;
					exit(1);
				}

				myfile.seekg(i);

	//			if (myfile.failbit) {
	//					cerr << "failbit" << endl;
	//					//exit(1);
	//			}
				c = myfile.peek();
				//cout << "c: \"" << c << "\"" << endl;
				i--;

			}
			c = myfile.get();
			//myfile.read(lastline, 5);
			getline(myfile, lastline);
		//while ( myfile.good() )
		//{
		//  getline (myfile,line);
		 // cout << line << endl;
		} else {
				cerr << "error: could not open file: " <<temp_output_file << " for validation" << endl;
				exit(1);
		}

		//cout << "lastline: \"" << lastline << "\"" << endl;

		if (lastline.substr(0,4).compare("S2: ") !=0 ) {
				cerr << "error: last line of BLAST file seems to be wrong..." << endl;
				cerr << "error: it is expected that the last line starts with \"S2: \"..." << endl;
				exit(1);
		}

		// so far so good, now search for runs of zeros
		myfile.clear();
		myfile.seekg(0, ios::beg);
		myfile.clear();

		int zerocounter = 0;
		while (myfile.good())     // loop while extraction from file is possible
		{
			c = myfile.get();       // get character from file
			if (myfile.good()) {
				if (c == 0x00) {
					zerocounter++;
					if (zerocounter == 2) {
						cerr << "error: BLAST file contains x00..." << endl;
						exit(1);
					}

				} else {
			//cout << c ;
					zerocounter = 0;
				}
			}
		}
		myfile.close();
	} // end of checkBLASToutput


	return 0;
}

BLAST_M9_Parser::~BLAST_M9_Parser(){
	//if (fp != 0) {
	//	int s = pclose(fp);
	//}
}

BLAST_M9_Parser::BLAST_M9_Parser(string file, string * zcat_bin){
	first_line_data=0;
	end_of_file=false;

	string zcat_command = *zcat_bin;
	zcat_command.append(" -f ");
	zcat_command.append(file);

	fp = popen(zcat_command.c_str(), "r");
	if (fp == NULL) {
		cerr << "error: "<< zcat_command << endl;
		exit(EXIT_FAILURE);
	}


	namespace io = boost::iostreams;

	io::file_descriptor_flags flags= io::close_handle;


	fpstream = new io::stream_buffer<io::file_descriptor_source>(fileno(fp), flags);

	blast_stream = new std::istream(fpstream);

	this->current_query = "";

}

bool BLAST_M9_Parser::getNextBLASTHits(vector<vector<string> * > ** blast_hits, string& query_id){
	//cerr << "first_line_data: " << first_line_data << endl;
	//if (first_line_data == 0) {
	if (end_of_file==true) {
		return false;
	}


	if ((*blast_hits) != 0) {
		cerr << "*blast_hits != 0" << endl;
		exit(1);
	}

	//this->current_query = "";

	*blast_hits = new vector<vector<string> * > ;
	if (first_line_data != 0) {
		//current_query = (*first_line_data)[0];
		(*blast_hits)->push_back(first_line_data);
	}



	int query_count = 0;


	first_line_data=0;
	string line;
	string query;
	while (*blast_stream) {
		std::getline(*blast_stream, line);
//cerr << "l: " << line << endl;
		if (line.length() == 0) {
			continue;
		}
		if (line[0] == '#') {

			if (line.compare(0, 9, "# Query: ") == 0) {
				query_count++;

				query_id = this->current_query;

				size_t start_q = line.find_first_not_of(" \t", 9);
				if (start_q == string::npos) {
						cerr << "error: start_q == string::npos" << endl;
				}

				size_t end_q = line.find_first_of(" \t", start_q);
				if (end_q == string::npos) {
					this->current_query =  line.substr(start_q);
				} else {
					this->current_query =  line.substr(start_q, end_q-start_q);
				}
				//cerr << "line: " << line << endl;
				//cerr << "\""<<this->current_query <<"\"" << endl;

//cerr << "current_query: " << this->current_query << endl;
//cerr << "query_count: " << query_count << endl;

//cerr << "XXXXXXXXXXXXXXXXXXXXXXX: \"" << this->current_query << "\""<<endl;

				if ((*blast_hits)->size() > 0 || (query_count > 1)) {
					return true;
				}

				if (query_id.size() > 0) {
					if (query_id.compare(this->current_query) != 0) {
							return true;
					}
				}

//cerr << "did not return" << endl;
			}

			continue;
		}
		//cerr << "do parse_column_data_line: " << line << endl;
		vector<string> * line_data = parse_column_data_line(line, '\t');

		while (line_data->size() > 12) { // TODO if fasta description contains tabulator, blast taublar output also contains this
			line_data->erase(line_data->begin()+1); // erase second element
		}

		//cerr << "done" << endl;
		query = line_data->at(0);
		//cerr << "query: " << query << endl;
		if ( query.compare(this->current_query)==0 || ((*blast_hits)->size() == 0)){

			if ((*blast_hits)->size() == 0) {
				this->current_query = query;
			}

			(*blast_hits)->push_back(line_data);

		} else {
			first_line_data = line_data;

			this->current_query = (*first_line_data)[0];

			query_id = query;

			return true;
		}

	}

	end_of_file = true;

	if (fp != 0) {
		int s = pclose(fp);
	}
	query_id = query;
	return true;
}
