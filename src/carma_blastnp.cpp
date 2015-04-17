
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


#include "carma_blastnp.hpp"


CARMA_BLASTNP::CARMA_BLASTNP(string * config_file) : CARMA_BASE(config_file, 0){
		this->query_sequences=0;
		this->rdp_unaligned=0;
		//this->rdp2taxid_mapping = 0;

		setVariables();
}

void CARMA_BLASTNP::setVariables(){




	this->rdp_unaligned=getConfigStringValue(parameters, "rdp_unaligned");
	//checkFileWarning(this->rdp_unaligned);

}


CARMA_BLASTNP::~CARMA_BLASTNP(){
	if (query_sequences != 0) {
		delete query_sequences;
	}

	if (rdp_unaligned != 0) {
		delete rdp_unaligned;
	}

}


// not used anymore
map<string, tax_id> * CARMA_BLASTNP::get_rdp2taxid_mapping(string * file) {
	map<string, tax_id> * rdp2taxid_mapping = new map<string, tax_id>();


	FASTA_Parser * rdp_parser = new FASTA_Parser(*file, false, zcat_bin);
	string line;
	string * identifier = 0;
	tax_id taxid;
	map<string, tax_id>::iterator map_it;
	while ( rdp_parser->getNextDescriptionLine(line) ) {
		taxid = 0;
//cerr << "line" << line << endl;
		parse_RDP_description_line(line, &identifier, &taxid, NULL, this->NameToNCBI );


		map_it = rdp2taxid_mapping->find(*identifier);
		if (map_it != rdp2taxid_mapping->end()) {
			cerr << "Identifier already in map: " << *identifier << endl;
			exit(1);
		}

		rdp2taxid_mapping->insert( pair<string, tax_id > (*identifier, taxid) );

		delete identifier;

	}



	return rdp2taxid_mapping;
}


map<int,  pair<tax_id, string * > * > *  CARMA_BLASTNP::initBlastDump() {

	map<int,  pair<tax_id, string * > * > * blast_dump = new map<int, pair<tax_id, string * > * >;

	string fastacmd_dump = *fastacmd_script;
	fastacmd_dump.append(" -d ");
	fastacmd_dump.append(*(this->blast_nr_database));
	fastacmd_dump.append(" -D 1 ");

	//cerr << fastacmd_dump << endl;
	string sep_string = " >gi|";
	FASTA_Parser * fastacmd_dump_parser = new FASTA_Parser(fastacmd_dump, false);
	fastacmd_dump_parser->expected_sequence_length=1500;
	string line;
	int count =0;
	while ( fastacmd_dump_parser->getNextDescriptionLine(line) ) {
		count++;
		if (count%1000000 == 0) {
			cerr << "count: " << count << endl;
		}
		//cerr << "line: " << line << endl;
		set<int> gi_list;
		tax_id taxid = 0;
		bool found_taxid=false;
		vector<string> * blastdb_data =  parse_column_data_line(line, sep_string);
		for (int i = 0; i < blastdb_data->size(); i++) {
			//cerr << i<<": "<< blastdb_data->at(i) << endl;

			// GI -numbers
			size_t gi_start = blastdb_data->at(i).find_first_of('|');
			if (gi_start == string::npos) {
				cerr << "error: gi_start == string::npos" << endl;
				exit(1);
			}
			size_t gi_end = blastdb_data->at(i).find_first_of('|', gi_start+1);
			if (gi_end == string::npos) {
				
				cerr << "warning: gi_end == string::npos" << endl;
				continue;
				
				cerr << "blastdb_data->at(i): " << blastdb_data->at(i) << endl;
				cerr << "line: " << line << endl;
				exit(1);
			}

			string gi_string = blastdb_data->at(i).substr(gi_start+1, gi_end-gi_start-1);
			//cerr << gi_string << endl;
			int gi_number = str2int(gi_string);
			gi_list.insert(gi_number);


			// taxonomy
			if (!found_taxid) {
				size_t tax_start = blastdb_data->at(i).find_first_of('[');
				if (tax_start == string::npos) {
					//cerr << "error: tax_start == string::npos" << endl;
					//exit(1);
				}
				size_t tax_end = blastdb_data->at(i).find_first_of(']', tax_start+1);
				if (tax_end == string::npos) {
					//cerr << "error: tax_end == string::npos" << endl;
					//exit(1);
				}
				if (tax_start != string::npos && tax_end != string::npos) {
					string tax_string = blastdb_data->at(i).substr(tax_start+1, tax_end-tax_start-1);
					//cerr << "tax: "<< tax_string << endl;
					taxid = taxonname2taxid(tax_string, this->NameToNCBI);
					//cerr << "tax_i: "<< (int) taxid << endl;
					found_taxid = true;
				}
			}

		}
		delete blastdb_data;

		if (gi_list.size()==0) {
			continue;
		}

		if (taxid == 0) {
			continue;
		}
//continue;
		string * seq = fastacmd_dump_parser->getSequence();

		if (seq->length() <= 5) {
				cerr << "line: " << line << endl;
				cerr << "error: database sequence length: " << seq->length() << endl;
				exit(1);
		}

		pair<tax_id, string * > * pair_tax_seq = new pair<tax_id, string * >(taxid, seq);
//cerr << seq->length() << endl;
//cerr << seq->capacity() << endl << endl;;
		map<int,  pair<tax_id, string * >  * >::iterator blast_dump_it;
		pair<map<int,  pair<tax_id, string * > * >::iterator,bool> dump_ret;

		set<int>::iterator set_it;
		//set_it=gi_list.begin();
		for (set_it=gi_list.begin(); set_it!=gi_list.end(); set_it++) {
			//cerr << *set_it << endl;
			//exit(1);
			dump_ret = blast_dump->insert(pair<int, pair<tax_id, string * > * >(*set_it, pair_tax_seq));
			if (dump_ret.second==false) {
				cerr << "element already existed: " << *set_it << endl;
				exit(1);
			}

			//gi_list
		}
		//cerr << *seq << endl;


		//exit(1);
	}

	return blast_dump;
}


void CARMA_BLASTNP::parse_blast_M9(string * blast_file, string * fasta_input_file, char database, int match_type){
//cerr <<"huhu: " << *fasta_input_file << endl;
//this->query_sequences=0;
//cerr <<"huhu: " << this->query_sequences << endl;

	if (fasta_input_file == 0) {
		cerr << "error: fasta_input_file == 0" << endl;
		exit(1);
	}

	this->query_sequences = readFastaFileIntoMap(fasta_input_file);


	blast_dump = 0;
	if (match_type == MATCH_TYPE_BLASTP) {
		blast_dump = initBlastDump();
		//exit(1);
	}



//cerr <<"huhu2" << endl;
	string * blast_database;

	switch (match_type) {
		case MATCH_TYPE_BLASTN:
			switch (database) {
				case 'n':
					blast_database = this->blast_nt_database;
					break;
				case 'r':
					blast_database = this->blast_rdp_database;
					//cerr << "read rdp2taxid_mapping...   ";
					//this->rdp2taxid_mapping = this->get_rdp2taxid_mapping(this->rdp_unaligned);
					cerr << "done." << endl << endl;
					break;
				default:
					cerr << "error: database \"" << database << "\" not supported." << endl;
					exit(1);
			}
			break;
		case MATCH_TYPE_BLASTP:
			blast_database = this->blast_nr_database;
			break;
		default:
			cerr << "match_type not supported" << endl;
			exit(1);

	}

	BLAST_M9_Parser * parser = new BLAST_M9_Parser(*blast_file, zcat_bin);

	vector<vector<string> * > * blast_hits = 0;
	string query_id;
	while (parser->getNextBLASTHits(&blast_hits, query_id)) {
//cerr << "+++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		//cerr << "q: "<< blast_hits->at(0)->at(0) << endl;
		//cerr << "number of hits: "<< blast_hits->size() << endl;
		if ( blast_hits->size() > 0 ) {
			classify_blastnp_results(blast_hits, blast_database, database, match_type);
		} else {

			string prop="carma_";
			prop.append(MatchTypeToString(match_type));

			string evalue_none = "-";

			ClassificationResult * classification_result = new ClassificationResult(query_id, prop, 0, 0, evalue_none );
			printClassificationResult(classification_result);
			delete classification_result;
		}


		DeleteContainerWithPointers(blast_hits);
		blast_hits=0;
//cerr << "-----------------------------------" << endl;
	}

}

MatchingQuery * CARMA_BLASTNP::blastn_m8_to_CARMA3(vector<vector<string> * > * blast_hits, string * blast_database, char database, int match_type) {
	MatchingQuery * matching_query = new MatchingQuery();

	Match_hit * blasthit = 0;
	Blast_HSP * blasthsp = 0;

	if (blast_hits->size() == 0) {
		cerr << "error: blast_hits->size() == 0" << endl;
		exit(1);
	}

	matching_query->hits = new vector<Match_hit * >;

	matching_query->match_type=match_type;
	matching_query->query = blast_hits->at(0)->at(0);
	//cerr << "matching_query->query: " << matching_query->query << endl;

	vector<vector<string> * >::iterator blast_it;

	string last_subject = "";
	map<string, tax_id>::iterator map_it;

	for (blast_it = blast_hits->begin(); blast_it != blast_hits->end(); ++blast_it) {

		vector<string> * line_data = *blast_it;
		string current_subject = line_data->at(1);

		if (current_subject.compare(last_subject) == 0) {
			// another HSP of previous hit, skip because we use only first HSP
			continue;
		}

		blasthit = new Match_hit();
		blasthit->blast_id = new string(current_subject);
		blasthsp = new Blast_HSP();
		blasthit->HSPs = new vector<Blast_HSP * >;
		blasthit->HSPs->push_back(blasthsp);

		tax_id taxid=0;
		string database_sequence;




		if (blast_dump != 0) {

			//cerr << *(blasthit->blast_id) << endl;

			size_t start_pipe = current_subject.find_first_of('|');
			if (start_pipe == string::npos) {
				cerr << "error: start_pipe == string::npos" << endl;
				exit(1);
			}

			size_t end_pipe = current_subject.find_first_of('|', start_pipe+1);
			if (end_pipe == string::npos) {
				cerr << "error: end_pipe == string::npos" << endl;
				exit(1);
			}

			string gi_string = current_subject.substr(start_pipe+1, end_pipe-start_pipe-1);
			//cerr << gi_string << endl;
			int gi_int = str2int(gi_string);

//gi_int = 12725253;

			map<int, pair<tax_id, string * > * >::iterator dump_it;
			dump_it = blast_dump->find(gi_int);

			if (dump_it == blast_dump->end()) {
				cerr << "warning: GI " << gi_int << " not found.." << endl;
				taxid = 0;
				database_sequence = "";

			} else {
				//cerr << "found" << endl;

				pair<tax_id, string * > * pp = dump_it->second;
				taxid = pp->first;
				database_sequence=*(pp->second);

				//cerr << taxid << endl;
				//cerr << database_sequence << endl;

				//exit(1);
			}

 		} else {

		//cerr << "*(blasthit->blast_id): "<< *(blasthit->blast_id) << endl;



			switch (match_type) {
				case MATCH_TYPE_BLASTP:
					taxid = extractTaxonomyInformationById(blasthit->blast_id, match_type, blast_database, fastacmd_script);
					break;
				case MATCH_TYPE_BLASTN:

					switch (database) {
						case 'n':
							{
							taxid = extractTaxonomyInformationById(blasthit->blast_id, match_type, blast_database, fastacmd_script);
							break;
							}
						case 'r':
							{


							vector<string> * data_field = parse_column_data_line(*(blasthit->blast_id), '_');
							if (data_field->size() != 2) {
								cerr << "error: data_field->size() != 2" << endl;
								exit(1);
							}
							taxid = str2int(data_field->at(1));

							delete data_field;

							break;
							}
						default:
							cerr << "error: database \"" << database << "\" not supported." << endl;
							exit(1);
					}
					break;

				default:
					cerr << "error: match_type" << endl;
					exit(1);
			}
 		}
		//cerr << "taxid: " << taxid << endl;

		blasthit->ncbi_tax_id = taxid;


		// Fields: (0) Query id, (1)Subject id, (2)% identity, (3)alignment length, (4)mismatches, (5)gap openings, (6)q. start, (7)q. end, (8)s. start, (9)s. end, (10)e-value, (11)bit score
		blasthsp->query_start = str2int(line_data->at(6));
		blasthsp->query_end = str2int(line_data->at(7));

		blasthsp->subject_start = str2int(line_data->at(8));
		blasthsp->subject_end = str2int(line_data->at(9));


		blasthsp->evalue = string2double(line_data->at(10));
		blasthsp->bitscore = string2double(line_data->at(11));

		blasthsp->alignment_length = str2int(line_data->at(3));
//		blasthsp->gapopen = str2int(line_data->at(5));
//		blasthsp->percent_identities = str2int(line_data->at(2));
//		blasthsp->mismatches = str2int(line_data->at(4));


		if (blast_dump != 0) {
			if (taxid !=0 ) {
				bool reverse = false;

				int s = blasthsp->subject_start;
				int e = blasthsp->subject_end;


				//cerr << s << endl;
				//cerr << e << endl;

				if ( s > e) {
					reverse = true;
					swap(s,e);
				}

				database_sequence = database_sequence.substr(s, e-s-1);
				if (reverse) {
					reverseComplementDNA(database_sequence);
				}
				blasthsp->subject_sequence = new string(database_sequence);
			}
			//cerr << database_sequence << endl;
			//exit(1);
		} else {
			blasthsp->subject_sequence = new string(
				extractExactSubstringFromBlastDB(blasthit->blast_id, blasthsp->subject_start, blasthsp->subject_end, matching_query->match_type, blast_database, fastacmd_script)
			);
		}

		matching_query->hits->push_back(blasthit);
		last_subject = current_subject;
	}

	blasthit = getFirstUnmaskedHit(matching_query->hits);
	blasthsp= blasthit->HSPs->front();


	map<string, string * >::iterator query_it;

	query_it = this->query_sequences->find(matching_query->query);
	if (query_it == this->query_sequences->end()) {
		cerr << "error: did not find query sequence for blast query \"" << matching_query->query<< "\""<< endl;
		cerr << "error: do the FASTA file and the BLAST file match?" << endl;
		exit(1);
	}
	string * query_sequence = query_it->second;

	blasthsp->query_sequence = new string(*query_sequence);
	matching_query->query_sequence = blasthsp->query_sequence;

	return matching_query;
}


void CARMA_BLASTNP::classify_blastnp_results(vector<vector<string> * > * blast_hits, string * blast_database, char database, int match_type){

	//string evalue_string = blast_hits->at(0)->at(10);

	//cerr << "--1" << endl;
	MatchingQuery * matching_query = blastn_m8_to_CARMA3(blast_hits, blast_database, database, match_type);
	//cerr << matching_query->hits->size() << endl;
//cerr << "--2" << endl;
	ClassificationResult * classification_result=0;
	this->CARMA_top_percent = 0;

	while (1){


		if (classification_result != 0) {
			delete classification_result;
			classification_result = 0;
		}

		if (this->CARMA_top_percent == 0) {
			this->CARMA_top_percent = 0.1;
		} else if (this->CARMA_top_percent == 0.1) {
			this->CARMA_top_percent = 0.3;
		} else if (this->CARMA_top_percent == 0.3) {
			this->CARMA_top_percent = 0.6;
		} else {
			this->CARMA_top_percent = 1;
		}
// do something
//cerr << "#this->CARMA_top_percent: " << this->CARMA_top_percent << endl;
		classification_result = filterMatchingQuery(matching_query, CARMA3);
	//Match_hit * mh = getFirstUnmaskedHit(matching_query->hits);
//cerr <<"front: " <<  matching_query->hits->front()->ncbi_tax_id << endl;
//cerr <<"front masked: " <<  matching_query->hits->front()->masked << endl;

		if (classification_result == 0) {

			bool doReduceTaxon=false;
			if (match_type == MATCH_TYPE_BLASTP) {
				doReduceTaxon = true;
			}

			classification_result = computeClassification(matching_query, doReduceTaxon);

			if (classification_result->error_no_worse_hit == false) {
				//good
				break;
			}
				// not good, repeat
		}

		if (this->CARMA_top_percent >= 1) {
			break;
		}

	}

	if (classification_result == 0) {


		//string prop = "BLASTN";
		string prop="carma_";
		prop.append(MatchTypeToString(match_type));

		string evalue_string = double2str(getFirstUnmaskedHit(matching_query->hits)->HSPs->front()->evalue);


		classification_result = new ClassificationResult(matching_query->query, prop, 0, 0, evalue_string);
	}

//cerr << "top_003: " << top_003 << endl;
//cerr << "top_01: " << top_01 << endl;
//cerr << "top_03: " << top_03 << endl;
//cerr << "top_06: " << top_06 << endl;
//cerr << "top_1: " << top_1 << endl;

	printClassificationResult(classification_result);


	delete classification_result;

	delete matching_query;

}


