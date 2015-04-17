
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

#include "carma_base.hpp"




int dotToGap(int c) {
		return (c=='.')?'-':c;
}

bool compare_hit_by_percent_identity ( const Match_hit * a, const Match_hit * b ) {
	return a->HSPs->front()->percent_identities > b->HSPs->front()->percent_identities;
}

bool compare_hit_by_bitscore ( const Match_hit * a, const Match_hit * b ) {
	return a->HSPs->front()->bitscore > b->HSPs->front()->bitscore;
}

bool compare_hit_by_evalue ( const Match_hit * a, const Match_hit * b ) {
	return a->HSPs->front()->evalue < b->HSPs->front()->evalue;
}




bool compare_hsp_by_bitscore ( const Blast_HSP * a, const Blast_HSP * b ) {
	return a->bitscore > b->bitscore;
}





// map "." to "-" for stockholm-alignment
void stringDotToGap(string& myString) {
	std::transform(myString.begin(), myString.end(), myString.begin(),
			   dotToGap);
}



// file from www.uniprot.org/docs/speclist
map<string, tax_id> * init_SPECLIST_TXT(const char * speclist_txt) {

	map<string, tax_id> * SwissProtToNCBI = new map<string, tax_id>;


	ifstream myfile (speclist_txt);
	if (myfile.is_open())
	{

		string line;

		// search start of data:
		while (! myfile.eof() )	{

			std::getline(myfile,line );
			if (line.length() >= 7 && line.substr(0,7).compare("_____ _") == 0) {
				break;
			}
		}

		if ( myfile.eof() ) {
			cerr<< "error: no data found in file" << endl;
			exit(EXIT_FAILURE);
		}

		while (! myfile.eof() )	{

			std::getline(myfile,line );

//cout << line << endl;
			if (line.length() == 0) {
				continue;
			}

			if (line.length() >= 5 && line.substr(0,5).compare("     ") == 0) {
				continue;
			}

			if (line.length() >= 5 && line.substr(0,5).compare("=====") == 0) {
				continue;
			}

			if (line.length() >= 1 && line.substr(0,1).compare("(") == 0) {
				continue;
			}
			// example:
			// "AADNV V 648330: N=Aedes albopictus densovirus (isolate Boublik/1994)"

			size_t space1 = line.find_first_of ( ' ', 0);
			if (space1==string::npos) {
				cerr << "Error: init_SPECLIST_TXT" << endl;
				exit(EXIT_FAILURE);
			}

			size_t colon = line.find_first_of ( ':', space1);
			if (colon==string::npos) {
				cerr << "Error: init_SPECLIST_TXT" << endl;
				exit(EXIT_FAILURE);
			}

			size_t space2 = line.find_last_of ( ' ', colon);
			if (space2==string::npos) {
				cerr << "Error: init_SPECLIST_TXT" << endl;
				exit(EXIT_FAILURE);
			}



			//if ((pipe4-pipe3) >= 17 && line.substr(pipe3, 17).compare("|\tscientific name")==0 ){

			  // cout << line << endl;
				string tax_string = line.substr(space2+1, colon-space2);
				tax_id my_tax_id = str2int(tax_string);
				//cout << tax_id << endl;

				//int name_start = pipe1+2;
				//int name_end = space1-1;

				string name = line.substr(0, space1);

				stringToLower(name);
				(*SwissProtToNCBI)[name]=my_tax_id;
				//cout << "-" << name << "-"<< my_tax_id << endl;
				//exit(0);
			//}

			//string tax_id = line.substr(0, l_start);



		}

		myfile.close();
	}

	else cout << "Unable to open file";

	return SwissProtToNCBI;
}


map<string, tax_id> * init_EMBL_CDS(const char * by_species_gr, map<string, tax_id> * NameToNCBI) {

	map<string, tax_id> * EMBL_CDS = new map<string, tax_id>;

	int total = 0;
	int known = 0;

	ifstream myfile (by_species_gr);
	if (myfile.is_open())
	{

		tax_id tax_node = 0;

		string line;

		while (! myfile.eof() )	{

			std::getline(myfile,line );

//cout << line << endl;
			if (line.length() == 0) {
				continue;
			}

			if (line.length() >= 2 && line.substr(0,2).compare("//") == 0) {
				tax_node = 0;
				continue;
			}

			if (line.length() >= 5 && line.substr(0,5).compare("ID   ") == 0) {
				string name = line.substr(5);
			   // cout << name << endl;
				total++;
				stringToLower(name);

				map<string, tax_id>::iterator map_it;
				map_it =  NameToNCBI->find(name);

				if (map_it == NameToNCBI->end()) {

					size_t bracket1 = name.find_first_of('(', 0);

					if (bracket1 != string::npos) {
					   //cerr << "origin: " << name << endl;

						name = name.substr(0,bracket1-1);
						map_it =  NameToNCBI->find(name);
						if (map_it == NameToNCBI->end()) {
							if (PRINT_UNKNOWN_EMBL_CDS_SPECIES == 1) {
								cerr << "Error parsing by_species.gr, species unknown: \"" << name << "\"" << endl;
							}

							tax_node = 0;
							continue;
						}

					}

				}

				known++;
				tax_node = (*map_it).second;

			}

			if (tax_node == 0) {
				continue;
			}


			//cout << "(" << line << ")" << endl;
			(*EMBL_CDS)[line]=tax_node;


		}
		if (PRINT_UNKNOWN_EMBL_CDS_SPECIES == 1) {
			cerr << known << " of " << total << " species in EMBL_CDS also found in NCBI taxonomy." << endl;
		}

		myfile.close();
	}

	else cout << "Unable to open file";

	return EMBL_CDS;
}






map<string, MatchingQuery * > * mergeHits(map<string, MatchingQuery * > * matching_query_A, map<string, MatchingQuery * > * matching_query_B){


//	double top_percent = TOP_PERCENT / 100;

	map<string, MatchingQuery * > * mq_merged = new map<string, MatchingQuery * >;

	map<string, MatchingQuery * >::iterator mq_A_it;
	map<string, MatchingQuery * >::iterator mq_B_it;

int a=0;
int b=0;

	for (mq_A_it = matching_query_A->begin(); mq_A_it != matching_query_A->end(); mq_A_it++) {

//		vector<Match_hit * > * hits_merged = 0;


		vector<Match_hit * > * hits_A = (*mq_A_it).second->hits;



		mq_B_it = matching_query_B->find( (*mq_A_it).second->query );

		if (mq_B_it == matching_query_B->end()) {
			cerr << "error, query not found..." << endl;
			exit(EXIT_FAILURE);
		}

		vector<Match_hit * > * hits_B = (*mq_B_it).second->hits;

		int a_num=0;
		int b_num=0;
		if (hits_A != 0) {
			a_num = hits_A->size();
			cerr << "A " << a_num << endl;
		} else {
			cerr << "A --" << endl;
		}
		if (hits_B != 0) {
			b_num = hits_B->size();
			cerr << "B " << hits_B->size() << endl;
		} else {
			cerr << "B --" << endl;
		}

		if (a_num <= 1&& b_num > 1) {
				b++;
		}

		if (b_num <= 1 && a_num > 1) {
				a++;
		}
//		if (hits_A == 0 && hits_B != 0) {
//
//			hits_merged = new vector<Match_hit * >;
//			vector<Match_hit * >::iterator hit_it;
//
//			for (hit_it = hits_B->begin(); hit_it != hits_B->end(); hit_it++){
//				Match_hit * hit = new Match_hit((*hit_it)->ncbi_tax_id, (*hit_it)->bitscore, (*hit_it)->evalue, (*hit_it)->positives);
//				hits_merged->push_back(hit);
//			}
//
//		} else if (hits_A != 0 && hits_B == 0) {
//
//			hits_merged = new vector<Match_hit * >;
//			vector<Match_hit * >::iterator hit_it;
//
//			for (hit_it = hits_A->begin(); hit_it != hits_A->end(); hit_it++){
//				Match_hit * hit = new Match_hit((*hit_it)->ncbi_tax_id, (*hit_it)->bitscore, (*hit_it)->evalue, (*hit_it)->positives);
//				hits_merged->push_back(hit);
//			}
//		} else if (hits_A == 0 && hits_B == 0) {
//
//
//		} else {
//			hits_merged = new vector<Match_hit * >;
//
//
//			if (hits_A==0) {exit(EXIT_FAILURE);}
//			if (hits_A->size()==0) {exit(EXIT_FAILURE);}
//			double max_score = hits_A->front()->bitscore;
//
//
//			if (hits_B==0) {exit(EXIT_FAILURE);}
//			if (hits_B->size()==0) {exit(EXIT_FAILURE);}
//
//			if (max_score < hits_B->front()->bitscore) {
//				max_score  = hits_B->front()->bitscore;
//			}
//
//			vector<Match_hit * >::iterator hit_it;
//
//
//			double min_bitscore = max_score*(1-top_percent);
//
//			for (hit_it = hits_A->begin(); hit_it != hits_A->end(); hit_it++){
//
//				if  ( (*hit_it)->bitscore >= min_bitscore ) {
//
//					if ((*hit_it)->ncbi_tax_id == 0) {
//						cerr << "zero" << endl;
//						exit(EXIT_FAILURE);
//					}
//
//					Match_hit * hit = new Match_hit((*hit_it)->ncbi_tax_id, (*hit_it)->bitscore, (*hit_it)->evalue, (*hit_it)->positives);
//					hits_merged->push_back(hit);
//				}
//			}
//
//			for (hit_it = hits_B->begin(); hit_it != hits_B->end(); hit_it++){
//
//				if  ( (*hit_it)->bitscore >= min_bitscore ) {
//
//					if ((*hit_it)->ncbi_tax_id == 0) {
//						cerr << "zero" << endl;
//						exit(EXIT_FAILURE);
//					}
//
//					Match_hit * hit = new Match_hit((*hit_it)->ncbi_tax_id, (*hit_it)->bitscore, (*hit_it)->evalue, (*hit_it)->positives);
//					hits_merged->push_back(hit);
//				}
//			}
//		}


//		MatchingQuery * matching_query = new MatchingQuery();
//		matching_query->query=(*mq_A_it).second->query;
//
//		matching_query->hits=hits_merged;
//		matching_query->position_in_blast_file=-1;
//
//		//matching_queries->push_back(matching_query);
//		mq_merged->insert(pair<string,MatchingQuery * >(matching_query->query,matching_query));




	}

cerr << "a (blastx): " << a << endl;
cerr << "b (blastn): " << b << endl;

	return mq_merged;
}



int reduceTaxon(Blast_HSP * blast_hsp, tax_id * buffer_taxid, CARMA_BASE * basic){


	string assignment_rank_string;

	assert(blast_hsp);



	int query_length = abs(abs(blast_hsp->query_start) - abs(blast_hsp->query_end));

	int gaps = blast_hsp->gaps;
	int per_id = blast_hsp->percent_identities;
	int total_plus_percent = blast_hsp->positives - blast_hsp->identities;


//cout << "query_length: " << query_length << endl;
//cout << "per_id: " << per_id << endl;
//cout << "total_plus_percent: " << total_plus_percent << endl;

	if (query_length >= 550) {
		if (per_id >= 90) {

			if (gaps == 0) {

				if (total_plus_percent == 0) {
					assignment_rank_string = "species";
				} else  {
					assignment_rank_string = "genus";
				}

			} else {
				assignment_rank_string = "genus";
			}



		} else if (per_id >= 75) {

			if (total_plus_percent <= 2) {
				assignment_rank_string = "genus";
			} else  {
				assignment_rank_string = "family";
			}



		} else if (per_id >= 60 ) {
			assignment_rank_string = "class";
			//assignment_rank_string = "order";
		} else if (per_id >= 45 ) {
			assignment_rank_string = "phylum";
			//assignment_rank_string = "class";
		} else { // identities_int <= 40 !
			assignment_rank_string = "no rank";
		}


	} else if (query_length >= 350) {



		if (per_id >= 65) {

			if (total_plus_percent <= 6) {
				assignment_rank_string = "genus";
			} else { // positives-identities > 6
				assignment_rank_string = "family";
			}
		} else if (per_id >= 60) {
			assignment_rank_string = "family";
		} else if (per_id >= 55) {
			//assignment_rank_string = "family"; # original
			assignment_rank_string = "order";
		} else if (per_id >= 50) {
			assignment_rank_string = "class";
		} else if (per_id >= 35) {
			assignment_rank_string = "phylum";
		} else { // identities_int <= 40 !
			assignment_rank_string = "no rank";
		}


	} else if (query_length >= 150) {


		if (per_id >= 65) {

			if (total_plus_percent <= 6) {
				assignment_rank_string = "genus";
			} else {
				assignment_rank_string = "family";
			}
		} else if (per_id >= 60) {
			assignment_rank_string = "family";
		} else if (per_id >= 55) {
			assignment_rank_string = "order";
		} else if (per_id >= 50) {
			assignment_rank_string = "class";
		} else if (per_id >= 45) {
			assignment_rank_string = "phylum";
		} else { // identities_int <= 40 !
			assignment_rank_string = "no rank";
		}


	} else {

		if (per_id >= 65) {

			if (total_plus_percent <= 1) {
				assignment_rank_string = "genus";
			} else {
				assignment_rank_string = "family";
			}
		} else if (per_id >= 50) {
			assignment_rank_string = "phylum";
		} else { // identities_int <= 40 !
			assignment_rank_string = "no rank";
		}


	}


//		// values according SOrt-ITEMS paper
//		if (per_id >= 66) {
//			if (blast_hsp->positives - blast_hsp->identities == 0) {
//				assignment_rank_string = "species";
//			} else if (blast_hsp->positives - blast_hsp->identities <= 6) {
//				assignment_rank_string = "genus";
//			} else { // positives-identities > 6
//				assignment_rank_string = "family";
//			}
//		} else if (per_id >= 61) {
//			assignment_rank_string = "family";
//		} else if (per_id >= 56) {
//			assignment_rank_string = "order";
//		} else if (per_id >= 51) {
//			assignment_rank_string = "class";
//		} else if (per_id >= 41) {
//			assignment_rank_string = "phylum";
//		} else { // identities_int <= 40 !
//			assignment_rank_string = "no rank";
//		}




	map<string, rank_t>::iterator map_it;
	map_it =  basic->rank_to_id->find(assignment_rank_string);

	if (map_it == basic->rank_to_id->end()) {
		cerr << "rank: " << assignment_rank_string << " not found... ?!" << endl;
		exit(EXIT_FAILURE);
	}


	rank_t assignment_rank_id = (*map_it).second;


	//cout << "assignment_rank_id: "  << (int) assignment_rank_id << endl;



	int last_valid_node = reduceTaxonByRank(buffer_taxid, basic->ranks, assignment_rank_id, basic->NCBI_MAX);

	return last_valid_node;

//cout << "----------- " << endl;


}









void mysigint(int signum)
{
    printf("I caught the SIGINT signal!\n");
    return;
}

/* Our own SIGKILL handler */
void mysigkill(int signum)
{
    printf("I caught the SIGKILL signal!\n");
    return;
}

/* Our own SIGHUP handler */
void mysighup(int signum)
{
    printf("I caught the SIGHUP signal!\n");
    return;
}

/* Our own SIGTERM handler */
void mysigterm(int signum)
{
    printf("I caught the SIGTERM signal!\n");
    return;
}







MatchingQuery * CARMA_BASE::parseBlastEntry(ifstream &myfile, ifstream::pos_type streamposition,  map<string, tax_id> * NameToNCBI, bool use_hard_threshold, bool report_only_first_hsp) {
	bool use_fastacmd = false;
	MatchingQuery * matching_query = 0;
	Match_hit * blasthit = 0;
	Blast_HSP * blasthsp = 0;
//cerr << "------------------------------------------------------------------------ A " << endl;

	string * blast_database;

	//ifstream myfile (blastx_file);
	if (myfile.is_open())
	{

		myfile.clear(); // in order to unset EOF-flag which might have been set in a previous run...
		myfile.seekg(streamposition);
//cerr << "go: " << (int) streamposition << endl;
		int match_type = MATCH_TYPE_UNKNOWN;
		blast_database = blast_nr_database;


		int HSP_counter=0;
		int HIT_counter=0;

		bool validHSP=false;
		bool validHit=false;
		//int query_count = 0;
		bool hit_already_added;
		string line;


		if (myfile.eof()) {
			cerr << "EOF !?!"<< endl;
			exit(EXIT_FAILURE);
		}
		int linecount = -1;
		while (! myfile.eof() )	{
			linecount++;
			std::getline(myfile,line );
			if (linecount == 0) {
//cout << "R: "<< line << endl;
			}


			if ( (mycompare(line, 0, 5, "BLAST")== 0) || (mycompare(line, 1, 5, "BLAST")== 0) ) {

				if (mycompare(line, 0, 6, "BLASTX")== 0) {
					match_type = MATCH_TYPE_BLASTX;
					blast_database = blast_nr_database;
				} else if (mycompare(line, 0, 6, "BLASTP")== 0) {
					match_type = MATCH_TYPE_BLASTP;
					blast_database = blast_nr_database;
				} else if (mycompare(line, 0,6,"BLASTN")== 0) {
					match_type = MATCH_TYPE_BLASTN;
					blast_database = blast_nt_database;
				} else {
					cerr << line << endl;
					cerr << "blast type not yet supported" << endl;
					exit(EXIT_FAILURE);
				}

				// end of BLAST entry
				if (linecount > 0) {
					assert(matching_query);
					return matching_query;
				}

				continue;
			}

			// new read
			if (mycompare(line, 0,6, "Query=")== 0) {

				matching_query = new MatchingQuery();
				if (match_type == MATCH_TYPE_UNKNOWN) {
						cerr << "error: match_type == UNKNOWN" << endl;
						exit(EXIT_FAILURE);
				}
				matching_query->match_type=match_type;

				string query_string = line.substr(7);
				size_t sep = query_string.find_first_of(" 	");

				if (sep != string::npos) {
					query_string = query_string.substr(0, sep);
				}

				matching_query->query = query_string;
				//matching_queries->insert(pair<string,MatchingQuery * >(matching_query->query,matching_query));
//cerr << "query: " << matching_query->query << endl;
				validHit = false;
				validHSP = false;

				blasthit=0;
				blasthsp=0;


				matching_query->position_in_blast_file = myfile.tellg();

				//query_count++;
				//if (query_count % 1000 == 0) {
				//	break;
				//	cerr << "\rQuery: " << query_count << "   ";
				//}
				continue;
			}

			// start of a hit (parse NCBI tax_id)
			if ((line.length() > 1) && (line[0] == '>') ) {
				if (blasthit != 0 && hit_already_added == false) {
					delete blasthit;
				}

				if (report_only_first_hsp == true && HIT_counter > 0) {
					return matching_query;
				}

				blasthit = new Match_hit();
				hit_already_added = false;

				HIT_counter++;
				HSP_counter=0;

				blasthit->ncbi_tax_id = 0;
				string previous_line="";

				int linecount = 0;

				while (1) {
					linecount++;
					if (linecount > 1) {

						std::getline(myfile,line );

						if (myfile.eof()) {
							cerr << "error: EOF came unexpectedly "<< endl;
							exit(EXIT_FAILURE);
						}

					}

					if (mycompare(line, 0, 16, "          Length")==0) {

						blasthit->ncbi_tax_id = extractTaxonomyInformation(previous_line, NameToNCBI, use_fastacmd, match_type, blast_database, fastacmd_script);

						if (blasthit->ncbi_tax_id != 0) {
							blasthit->blast_id = extractSequenceID(previous_line);
						}


						if ((PRINT_UNPARSABLE_BLAST_DESCRIPTIONS == 1) && (blasthit->ncbi_tax_id == 0)) {
								cerr << "nothing found in: " << previous_line << endl;
						}
						break;
					}



					if (mycompare(line, 0, 5, "     ")==0) {

						size_t line_start = line.find_first_not_of(' ', 0);
						if (line_start==string::npos) {
							cerr << "error: no non-space found !?!" << endl;

							exit(EXIT_FAILURE);
						}
						previous_line += " " + line.substr(line_start);

					} else {
						// this is a new line, the old one can be parsed now:

						if (previous_line.length() > 0) {


							blasthit->ncbi_tax_id = extractTaxonomyInformation(previous_line, NameToNCBI, use_fastacmd, match_type, blast_database, fastacmd_script);

							if (blasthit->ncbi_tax_id != 0) {
								blasthit->blast_id = extractSequenceID(previous_line);
							}


							if ((PRINT_UNPARSABLE_BLAST_DESCRIPTIONS == 1) && (blasthit->ncbi_tax_id == 0)) {
								cerr << "nothing found in: " << previous_line << endl;
							}
							if ( blasthit->ncbi_tax_id != 0) {
								break;
							}
						}
						previous_line=line;
					}
				} // inner while loop to find tax_id.



				if (blasthit->ncbi_tax_id == 0) {
					validHit=false;
				} else {
					validHit=true;
				}
				continue;
			}

			if (!validHit) {
					continue;
			}

			// start of a HSP
			if (mycompare(line, 0, 8, " Score =")== 0 ) {

				if (report_only_first_hsp == true && HSP_counter > 0) {
					return matching_query;
				}

				blasthsp = 0;

				double bitscore = extractBLASTBitScore(line);

					// -------------------- E-value

				double evalue = extractBLASTEvalue(line);

				//double bitscore_threshold = matching_query->best_hit_bitscore*(1-top_percent);




			//cout << "--"<< endl;
//cout << matching_query->best_hit_bitscore << endl;
//cout << bitscore_threshold << endl;
				// check if I have to add the hit..
				if (hit_already_added == false) {

					if (matching_query->hits==0) {
						matching_query->hits = new vector<Match_hit * >;
					}

					// add hit
					if ( bitscore > matching_query->best_hit_bitscore ) {
						matching_query->best_hit_bitscore=bitscore;
						//matching_query->best_length = length_int;
						matching_query->hits->insert( matching_query->hits->begin(), blasthit);
					} else  {
						matching_query->hits->push_back(blasthit);
					}
					hit_already_added = true;
				}

//cout << bitscore << endl;



				validHSP=true;
				blasthsp = new Blast_HSP();
				blasthsp->bitscore=bitscore;
				blasthsp->evalue=evalue;

				if (blasthit->HSPs == 0) {
						blasthit->HSPs = new vector<Blast_HSP * >();
				}

				blasthit->HSPs->push_back(blasthsp);

				//if (bitscore > blasthit->best_hsp_bitscore) {
				//	blasthit->best_hsp_bitscore = bitscore;
				//}
				if ( bitscore > matching_query->best_hit_bitscore ) {
					matching_query->best_hit_bitscore=bitscore;
				}
				HSP_counter++;





				continue;
			}

			// continue HSP
			if ( (mycompare(line, 0, 13, " Identities =")== 0) && (validHSP==true) ) {
				// example: " Identities = 17/23 (73%), Positives = 18/23 (78%)"

				assert(blasthsp);

				blasthsp->identities         = (int) extractIdentities(line);
				blasthsp->percent_identities = (int) extractPercentIdentities(line);
				blasthsp->positives          = (int) extractPositives(line);
				blasthsp->gaps 				 = (int) extractGaps(line);

				blasthsp->alignment_length = (int) extractAlignmentLength(line);

				//cout << line << endl;
				//cout << "blasthsp->alignment_length: " << blasthsp->alignment_length << endl;

				//exit(0);

				continue;
			}

			if ( (mycompare(line, 0, 8, " Frame =")== 0) && (validHSP==true) ) {
				blasthsp->frame = extractFrame(line);
			}

			if ( (mycompare(line, 0, 9, " Strand =")== 0) && (validHSP==true) ) {

				if (mycompare(line, 0, 21, " Strand = Plus / Plus") == 0) {
					blasthsp->frame = 1;
				} else if (mycompare(line, 0, 22, " Strand = Plus / Minus") == 0) {
					blasthsp->frame = -1;
				} else {
						cerr << line << endl;
						exit(EXIT_FAILURE);
				}
			}

			// continue HSP
			if (mycompare(line, 0, 6, "Sbjct:")== 0 && (validHSP==true)) {

				string subject_sequence = extractSequence(line);
				int start = extractSequenceStartNumber(line);
				int end = extractSequenceEndNumber(line);

				if (blasthsp->subject_sequence == 0) {
					blasthsp->subject_sequence = new string(subject_sequence);
					blasthsp->subject_start = start;
					blasthsp->subject_end = end;
				} else {
					blasthsp->subject_sequence->append(subject_sequence);
					// subject should always be increasing, therfore I only update the end:
					blasthsp->subject_end = end;
				}

				continue;
			}

			// continue HSP
			if ( mycompare(line, 0, 6, "Query:")== 0 && (validHSP==true) ) {

				string query_sequence = extractSequence(line);

				//int frameshiftnum = getFrameShiftNum(line);
				//blasthsp->query_frameshifts += frameshiftnum;

				int start;
				int end;

				if (blasthsp->frame > 0) {
					start = extractSequenceStartNumber(line);
					end = extractSequenceEndNumber(line);
				} else {
					end = extractSequenceStartNumber(line);
					start = extractSequenceEndNumber(line);
				}

				if (blasthsp->query_sequence == 0) {
					blasthsp->query_sequence = new string(query_sequence);
					blasthsp->query_start = start;
					blasthsp->query_end = end;
				} else {
					blasthsp->query_sequence->append(query_sequence);

					if (blasthsp->frame < 0) {
						blasthsp->query_start = start;
					} else {
						blasthsp->query_end = end;
					}
				}



				continue;
			}

		} // end of while
//cerr << "------------------- B " << endl;

		return matching_query;

	} else {

		cerr << "error file is not open" << endl;
		exit(EXIT_FAILURE);
		//return 0;
	}

	cerr << "error: should not have reached this point" << endl;
	exit(EXIT_FAILURE);
	return 0;

}





void CARMA_BASE::classifyBlastXOnly(string * blast_file, string * blast_egts, bool use_hard_threshold, enum Algo_Type algo_type, bool local_temp) {




	ifstream::pos_type blast_stream_position;

	MatchingQuery * matching_query_blast;

	//string * blast_file = 0;


	assert(blast_file);


	//map<string, ifstream::pos_type > * blast_index = 0;

	assert(this->blast_index);
	FILE * temporary_blastegts_file;

	fstream blast_egt_file;
	string blast_egt_temp_filename;

	if (blast_egts != 0) {


		if (local_temp) {
			temporary_blastegts_file = tmpfile();
			int temporary_blastegts_file_no=fileno(temporary_blastegts_file);
			blast_egt_temp_filename = "/proc/self/fd/";
			blast_egt_temp_filename.append(int2str(temporary_blastegts_file_no));
			//blast_egts_temp = new string (temp_blastegts_file_string);

		} else  {


			blast_egt_temp_filename = *(blast_egts);
			blast_egt_temp_filename.append(".part");

			if (FileExists((char *) blast_egt_temp_filename.c_str())) {
				remove(blast_egt_temp_filename.c_str());
			}
		}

		blast_egt_file.open(blast_egt_temp_filename.c_str(), ios::out);

		if (! blast_egt_file.is_open()) {
			cerr << "Error opening file: " << blast_egt_temp_filename  << endl;
			exit(1);
		}


	}
	bool already_printed_query;

	ifstream blast_stream (blast_file->c_str());

	if (not blast_stream.is_open()) {
			cerr << "error reading file: " << *(blast_file) << endl;
			exit(EXIT_FAILURE);
	}

//int top_003 =  0;
//int top_01 =  0;
//int top_03 =  0;
//int top_06 =  0;
//int top_1 =  0;

	map<string, ifstream::pos_type >::iterator blast_index_iterator;
//cerr << "SIZE: " << blastx_index->size() << endl;
	for (blast_index_iterator = blast_index->begin(); blast_index_iterator != blast_index->end(); blast_index_iterator++ ) {

		string query = blast_index_iterator->first;

		blast_stream_position = blast_index_iterator->second;


		matching_query_blast=0;

		ClassificationResult * classification_result = 0;


		blast_stream_position = blast_index_iterator->second;

		//if (blast_stream_position != -1) {
		if (blast_stream_position != std::ifstream::pos_type(-1)) {


			if (algo_type != CARMA3) {
				matching_query_blast = parseBlastEntry(blast_stream, blast_stream_position, NameToNCBI, use_hard_threshold, false);

				assert(matching_query_blast);

				classification_result = filterMatchingQuery(matching_query_blast, algo_type);
			}

			//if (matching_query_blast->hits != 0) {
				//sort(matching_query_blast->hits->begin(), matching_query_blast->hits->end(), compare_by_percent_identity);
			//}

			//cout << "matching_query_blast->hits contains:";
			//  for ( it=matching_query_blast->hits->begin() ; it < matching_query_blast->hits->end(); it++ ) {
			//  	Match_hit * testhit = *it;
			//    cout << "PI: " << testhit->HSPs->front()->percent_identities << endl;
			//  }



			if (classification_result != 0 ) {
					printClassificationResult(classification_result);
					delete classification_result;
					delete matching_query_blast;
					continue;
			}


		} else {

			continue;
		}


		Match_hit * mh;
		Blast_HSP * hsp;
		string evalue_string;

		if (algo_type != CARMA3) {
			assert(matching_query_blast);
			assert(matching_query_blast->hits);
			mh = getFirstUnmaskedHit(matching_query_blast->hits);
			assert(mh);
			assert(mh->HSPs);

			hsp = mh->HSPs->front();
			evalue_string = double2str(hsp->evalue);
		}
		switch (algo_type) {

			case CARMA3:



				this->CARMA_top_percent = 0;

				already_printed_query= false;

				matching_query_blast = parseBlastEntry(blast_stream, blast_stream_position, NameToNCBI, use_hard_threshold, false);
//cerr << "have: " <<  matching_query_blast->hits->front()->ncbi_tax_id << endl;
				assert(matching_query_blast);

				// top-percent: 0.1, 0.3, 0.6, 1
				while (1) {


//					if (matching_query_blast != 0) {
//						delete matching_query_blast;
//						matching_query_blast = 0;
//					}

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

					//cout << "#CARMA_top_percent: "<< CARMA_top_percent << endl;


					if (matching_query_blast->hits == 0){
						evalue_string = "0";
						mh=0;
						hsp=0;
					} else {
						mh = getFirstUnmaskedHit(matching_query_blast->hits);
						if (mh != 0) {
							assert(mh->HSPs);

							hsp = mh->HSPs->front();
							evalue_string = double2str(hsp->evalue);
						} else {
							evalue_string = "0";
							mh=0;
							hsp=0;
						}
					}



					classification_result = filterMatchingQuery(matching_query_blast, algo_type);

					if (classification_result == 0) {

						classification_result = computeClassification(matching_query_blast, true);

						if (classification_result->error_no_worse_hit == false) {
								//cout << "#3 / good "<< CARMA_top_percent << endl;
//								if (this->CARMA_top_percent == 0.03) {
//									top_003++;
//								}else if (this->CARMA_top_percent == 0.1) {
//									top_01++;
//								} else if (this->CARMA_top_percent == 0.3) {
//									top_03++;
//								} else if (this->CARMA_top_percent == 0.6) {
//									top_06++;
//								} else {
//									top_1++;
//								}

								break;
						}

						if (matching_query_blast->top_percent_exhausted) {
							//cerr << "XXXXXXXXXXXXXXXXXXXXXXX" << endl;
							break;
						}
						//cout << "#3 / retry next top-percent" << endl;
					}

					if (already_printed_query == false && blast_egts != 0 && hsp != 0) {

						assert(hsp->query_sequence);
						already_printed_query = true;
						blast_egt_file << "> " << query <<  endl;
						blast_egt_file << *(hsp->query_sequence) << endl;
					}



					if (this->CARMA_top_percent >= 1) {
							break;
					}
				}

//cerr << "stop" << endl;
			break;

			case LCA:
				{
				tax_id lca_result = computeLCA(matching_query_blast->hits);
				string property = "LCA";
				classification_result = new ClassificationResult(matching_query_blast->query, property, lca_result, 0, evalue_string);
				}
			break;

			case BESTBLAST:
				{
				string property = "best_hit";
				classification_result = new ClassificationResult(matching_query_blast->query, property, getFirstUnmaskedHit(matching_query_blast->hits)->ncbi_tax_id, 0, evalue_string);
				printClassificationResult(classification_result);
				delete classification_result;
				}
			break;
		}



		if (classification_result == 0) {

			string prop = "carma_";
			prop.append(MatchTypeToString(matching_query_blast->match_type));

			classification_result = new ClassificationResult(query, prop, 0, 0, evalue_string);
		}

//cerr << "top_003: " << top_003 << endl;
//cerr << "top_01: " << top_01 << endl;
//cerr << "top_03: " << top_03 << endl;
//cerr << "top_06: " << top_06 << endl;
//cerr << "top_1: " << top_1 << endl;

		printClassificationResult(classification_result);

		delete matching_query_blast;
		delete classification_result;


	}  // end for-loop

//cerr << "top_003: " << top_003 << endl;
//cerr << "top_01: " << top_01 << endl;
//cerr << "top_03: " << top_03 << endl;
//cerr << "top_06: " << top_06 << endl;
//cerr << "top_1: " << top_1 << endl;


	blast_stream.close();

	if (blast_egts != 0) {



		if (local_temp) {
			fflush(temporary_blastegts_file);
			bool copy_ret = copyFile( blast_egt_temp_filename.c_str(), blast_egts->c_str());
			if (copy_ret == false) {
				cerr << "(classifyBlastXOnly) error copying file from " <<  blast_egt_temp_filename << " to " << *blast_egts << endl;
				exit(1);
			}
		} else {
			blast_egt_file.close();
			int rename_ret = rename ( blast_egt_temp_filename.c_str(), blast_egts->c_str());
			if (rename_ret != 0) {
					cerr << "error on renameing file" << endl;
			}
		}

	}


}



ReciprocResults * CARMA_BASE::parseReciprocBLASTNP(tax_id best_hit_tax_id, rank_t * ranks, char * reciprok_blast_results_file_name, int match_type){
	string file = reciprok_blast_results_file_name;

	bool debug_mode = false;
//string cmd = "cat ";
//cmd.append(file);
//system(cmd.c_str());
	//exit(0);

	ReciprocResults * rResults = new ReciprocResults();
	rResults->match_type = match_type;



	tax_id * best_hit_taxonomy=0;
	tax_id * other_hit_taxonomy=0;

	getTaxa(best_hit_tax_id, &best_hit_taxonomy, this->parent_taxid, NCBI_MAX);

	BLAST_M9_Parser * blastparser = new BLAST_M9_Parser(file, zcat_bin);


	tax_id current_taxid = 0;
	rank_t lct_rank = 0;

	vector<vector<string> * > * blast_hits = 0;
	string current_sbj="";
	string query_id;
	blastparser->getNextBLASTHits(&blast_hits, query_id);

	//double max_score=DBL_MIN;
	double min_score=DBL_MAX;

	for (int i =0 ; i < blast_hits->size() ; ++i) {
		string subject = blast_hits->at(i)->at(1);
		//cerr << "q: "<< blast_hits->at(i)->at(0) << endl;
		//cerr << "s: "<< blast_hits->at(i)->at(1) << endl;


		if (subject.compare(current_sbj) != 0) {
//cerr <<"i: "<< i << endl;





			double evalue = string2double(blast_hits->at(i)->at(10));
			double bitscore = string2double(blast_hits->at(i)->at(11));
//cerr << "bitscore: "<<  bitscore << endl;
//cerr << "rResults->query_bitscore: " << rResults->query_bitscore << endl;
//cerr << "rResults->self_bitscore: " << rResults->self_bitscore << endl;
			if (bitscore < min_score) {
				min_score = bitscore;
			}
//			if (bitscore > max_score) {
//				max_score = bitscore;
//			}
//cerr << bitscore << " " << subject << endl;
			if (subject.compare("query")==0) {
				rResults->query_bitscore = bitscore;
				if (debug_mode) {
						cout << "#\t" << "query" << "\t" << 8 << "\t" << bitscore << endl;
				}
			} else if (subject.compare("self")==0) {
				rResults->self_bitscore= bitscore;
				if (debug_mode) {
						cout << "#\t" << "self" << "\t" << 0 << "\t" << bitscore << endl;
				}
			} else {



				if (bitscore >= this->carma_blastn_bitscore) {

					current_taxid = (tax_id)str2int(subject) ;
					getTaxa(current_taxid, &other_hit_taxonomy, this->parent_taxid, NCBI_MAX);



					tax_id lct = getLowestCommonTaxID(best_hit_taxonomy, other_hit_taxonomy, NCBI_MAX);
//	cerr << "last common tax_id: " << lct << endl;
					lct_rank = getLowestRank(lct, parent_taxid, ranks, NCBI_MAX);

					if (debug_mode) {
						int simpleRank = mapRankToSimpleRank(lct_rank);
						cout << "#\t" << (int)current_taxid << "\t" << simpleRank << "\t" << bitscore << endl;
					}

//	cerr << "lct_rank: " << (int)lct_rank << endl;
					pair<double, rank_t> * ppp = new pair<double, rank_t>(bitscore, lct_rank);

					if (rResults->other_hits == 0) {
						rResults->other_hits = new vector<pair<double, rank_t> * >();
					}

					rResults->other_hits->push_back(ppp);

				}

			}




		}
		current_sbj = subject;
	}
	//cerr << "number of hits: "<< blast_hits->size() << endl;
//exit(0);






	if (best_hit_taxonomy != 0 ){
		delete [] best_hit_taxonomy;
	}
	if (other_hit_taxonomy != 0 ){
		delete [] other_hit_taxonomy;
	}

	//if (blast_hits->size() < 250) {
	//if (true) {
		if (rResults->query_bitscore==-1) {
			//cout << "warning: rResults->query_bitscore: " << rResults->query_bitscore << endl;

			rResults->query_bitscore=min_score;
			//exit(1);
		}

		if (rResults->self_bitscore == -1) {
			cerr << "rResults->self_bitscore: " << rResults->self_bitscore << endl;
			exit(1);
		}
	//} else {
//		if (rResults->query_bitscore==-1) {
//			rResults->query_bitscore = min_score;
//
//		}
//
//		if (rResults->self_bitscore == -1) {
//			rResults->self_bitscore = max_score;
//		}
//	}


	delete blast_hits;
	delete blastparser;


	return rResults;
	//exit(1);
}

ReciprocResults * CARMA_BASE::parseReciprocBLASTP(tax_id best_hit_tax_id, rank_t * ranks, char * reciprok_blast_results_file_name, int match_type){

	tax_id * best_hit_taxonomy=0;
	tax_id * other_hit_taxonomy=0;

	ReciprocResults * rResults = new ReciprocResults();
	rResults->match_type = match_type;

	getTaxa(best_hit_tax_id, &best_hit_taxonomy, this->parent_taxid, NCBI_MAX);

//	int i = 0;
//	while (best_hit_taxonomy[i] != -1) {
//		cout << i << "  "<< best_hit_taxonomy[i] << endl;
//		i++;
//	}

	//exit(0);
	int hit_type = 0; // 0:other 1: self 2:query
	bool newhit=false;
	tax_id current_taxid = 0;
	rank_t lct_rank = 0;
	string blast_file = reciprok_blast_results_file_name;

	ifstream myfile (blast_file.c_str());
	if (myfile.is_open())
	{

		string line;

		while (! myfile.eof() )	{
			std::getline(myfile,line );
			//cerr << "L: " << line << endl;

			if (line.length() == 0) {
				continue;
			}

			if (line[0] == '>') {
				newhit=true;
				//cout << "L: " << line << endl;
				//cout << line.substr(1) << endl;

				if (mycompare(line, ">query") == 0) {
					//cout << "query" << endl;
					hit_type = 2;
				} else if (mycompare(line, ">self") == 0) {
					//cout << "self" << endl;
					hit_type = 1;
				} else {
					hit_type = 0;
					string bla = line.substr(1);
					current_taxid = (tax_id)str2int(bla) ;

//cerr << "line: " <<  line << endl;
//cerr << "current_taxid: " <<  (int)current_taxid << endl;

//					if (current_taxid == 0) {
//						cout << line << endl;
//						exit(EXIT_FAILURE);
//					}

					//cout << "compare: " << best_hit_tax_id << " and "<<current_taxid << endl;

					getTaxa(current_taxid, &other_hit_taxonomy, this->parent_taxid, NCBI_MAX);

					tax_id lct = getLowestCommonTaxID(best_hit_taxonomy, other_hit_taxonomy, NCBI_MAX);

					//cout << "last common tax_id: " << lct << endl;
					lct_rank = getLowestRank(lct, parent_taxid, ranks, NCBI_MAX);
					//cout << "last common lct_rank " << (int)lct_rank << endl;

				}

				//exit(0);
				continue;
			}


			if (newhit == true && mycompare(line, 0, 8, " Score =")==0) {
			//if (newhit == true && mycompare(line, 0, 13, " Identities =")==0) {

				double bitscore = extractBLASTBitScore(line);
				//double bitscore = extractPercentIdentities(line);

				//cout << line << endl;
				//cout << bitscore << endl;

				if (hit_type == 0) {

					//cout << bitscore << " - " << (int) lct_rank << endl;
					double bit_threshold = 0;

					if (match_type == MATCH_TYPE_BLASTX) {
						bit_threshold = this->carma_blastx_bitscore;
					} else if (match_type == MATCH_TYPE_BLASTN) {
						bit_threshold = this->carma_blastn_bitscore;
					} else {
						cerr << "blast_type not supported" << endl;
						exit(1);
					}


					if (bitscore >= bit_threshold) {
						pair<double, rank_t> * ppp = new pair<double, rank_t>(bitscore, lct_rank);

						if (rResults->other_hits == 0) {
							rResults->other_hits = new vector<pair<double, rank_t> * >();
						}

						rResults->other_hits->push_back(ppp);
					}
				} else if (hit_type == 1)  {
						rResults->self_bitscore=bitscore;
				} else if (hit_type == 2)  {
						rResults->query_bitscore=bitscore;
				}

				newhit = false;
				continue;
			}

			//cout << line << endl;
		}


	} else {
		cerr << "Error: Unable to open file " << blast_file << endl;
		exit(1);
	}

	myfile.close();

	if (best_hit_taxonomy != 0 ){
		delete [] best_hit_taxonomy;
	}
	if (other_hit_taxonomy != 0 ){
		delete [] other_hit_taxonomy;
	}


	return rResults;
}






tax_id CARMA_BASE::computeLCA(vector<Match_hit * > * hits) {

	assert(hits);
	assert(getFirstUnmaskedHit(hits)->HSPs);

	//double min_score = (1-lca_top_percent) * hits->front()->HSPs->front()->bitscore;

	tax_id one = 1;
	PhyloNode * tree = new PhyloNode(one);
	PhyloNode * treePointer;


	tax_id * buffer_taxid=0;

	vector<Match_hit * >::iterator hit_it;

	// filter "Other" and "Unclassified" and Hits below percent threshold
	for (hit_it = hits->begin(); hit_it != hits->end(); hit_it++){

		Match_hit * match_hit = *hit_it;



		getTaxa(match_hit->ncbi_tax_id, &buffer_taxid, this->parent_taxid, NCBI_MAX);

		//# skip "Other" and "Unclassified"
		if (buffer_taxid[1]== 28384 || buffer_taxid[1] == 12908) {
			continue;
		}

		treePointer = tree;

		int node_pos=1;
		tax_id ncbi_tax_id;
		while ((ncbi_tax_id = buffer_taxid[node_pos]) != -1) {

			treePointer = treePointer->getNode(ncbi_tax_id);
			node_pos++;
		}
		treePointer->addHit(*hit_it);

	}



	if (tree->children!=0) {



		tax_id * lca = tree->getLCA();
		int node_pos = 0;
		while (lca[node_pos] != -1) {
			//cout << "lca["<< node_pos <<"] =" << lca[node_pos] << endl;
			node_pos++;

			if (node_pos >= MAX_TREE_DEPTH) {
				cerr << "error: MAX_TREE_DEPTH" << endl;
				exit(EXIT_FAILURE);
			}

		}

		tax_id classification_result = 0;
		if (node_pos > 0) {
			//cout << "LCA: " << lca[node_pos-1] << endl;
			classification_result = lca[node_pos-1];
		}

		//if ( (*mq_it)->hits->size() < 6) {
//						cout << (*mq_it)->query << "\tpfam1\tgoterms\t" << classification_result << endl;
//					} else {
//						cout << (*mq_it)->query << "\tpfam2\tgoterms\t" << classification_result << endl;
//					}



		delete [] lca;
		tree->setHitsZero();
		delete tree;
		return classification_result;
	} else {
			return -1;
	}


}

CARMA_BASE::~CARMA_BASE(){

	if (rank_to_id != 0) {
		delete rank_to_id;
	}

	if (NameToNCBI != 0) {
		delete NameToNCBI;
	}

	if (taxid_to_name != 0) {
		delete taxid_to_name;
	}

	if (pfam_to_GOid != 0) {
		delete pfam_to_GOid;
	}

	if (parent_taxid != 0) {
		delete [] parent_taxid;
	}

	if ( ranks != 0) {
		delete [] ranks;
	}

	if (nucleobase_ascii2num != 0) {
		delete [] nucleobase_ascii2num;
	}

	if (buffer_taxid != 0) {
		delete [] buffer_taxid;
	}

	if (path != 0) {
		delete [] path;
	}

	if (parameters != 0) {
		delete parameters;
	}

	if (blast_index != 0) {
		delete blast_index;
	}



//------------------------

	if (nodes_dmp != 0) {
		delete nodes_dmp;
	}

	if (merged_dmp != 0) {
		delete merged_dmp;
	}

	if (names_dmp != 0) {
		delete names_dmp;
	}


	if (blast_nr_database != 0) {
		delete blast_nr_database;
	}

	if (blast_nt_database != 0) {
		delete blast_nt_database;
	}

	if (pfamId2TaxId_file != 0) {
		delete pfamId2TaxId_file;
	}

	if (blastall_script != 0) {
		delete blastall_script;
	}

	if (fastacmd_script != 0) {
		delete fastacmd_script;
	}

	if (formatdb_script != 0) {
		delete formatdb_script;
	}

	if (hmmfetch_bin != 0) {
		delete hmmfetch_bin;
	}

	if (hmmalign_bin != 0) {
		delete hmmalign_bin;
	}

	if (hmmscan_bin != 0) {
		delete hmmscan_bin;
	}

	if (pfam_A_hmm_file != 0) {
		delete pfam_A_hmm_file;
	}

	if (pfam_fasta_dir != 0) {
		delete pfam_fasta_dir;
	}

	if (blosum_file != 0) {
		delete blosum_file;
	}

	if (pfamA_txt_file != 0) {
		delete pfamA_txt_file;
	}

	if (gene_ontology_txt_file != 0) {
		delete gene_ontology_txt_file;
	}


	if (blast_rdp_database != 0) {
		delete blast_rdp_database;
	}

	if (gzip_bin != 0) {
		delete gzip_bin;
	}

	if (zcat_bin != 0) {
		delete zcat_bin;
	}

	//if ( != 0) {
	//	delete ;
	//}

}

CARMA_BASE::CARMA_BASE(string * config_file, string * config_overlay_string){

	this->nodes_dmp = 0;
	this->merged_dmp = 0;
	this->names_dmp = 0;
	this->blast_nr_database = 0;
	this->blast_nt_database = 0;
	this->pfamId2TaxId_file = 0;
	this->blastall_script = 0;
	this->fastacmd_script = 0;
	this->formatdb_script = 0;
	this->hmmfetch_bin = 0;
	this->hmmalign_bin = 0;
	this->hmmscan_bin = 0;
	this->pfam_A_hmm_file = 0;
	this->pfam_fasta_dir = 0;
	this->blosum_file = 0;
	this->pfamA_txt_file = 0;
	this->gene_ontology_txt_file = 0;

	this->blast_rdp_database = 0;
	this->NameToNCBI = 0;
	this->parameters = 0;
	this->blast_index=0;

	this->rank_to_id=0;
	//this->blast_egts=0;
	this->parent_taxid = 0;

	this->buffer_taxid=0;
	this->taxid_to_name=0;
	this->ranks = 0;
	this->gzip_bin = 0;
	this->zcat_bin = 0;
	//this->blast_file = 0;


	//this->dna_fasta_file = 0;

	this->pfamid2taxid_map=0;
	this->pfam_to_GOid=0;

	this->blosum_matrix = 0;

	this->filtered_taxa=0;



	this->blast_egts_map=0;

	this->readConfiguration(config_file);
	// overwrite using config_overlay_string

	
	if (config_overlay_string != 0 ) {
		vector<string > * config_overlay_vector = parse_column_data_line(*config_overlay_string, ',');
		vector<string > * config_element = 0;
		map<string, string>::iterator map_it;
		for (int i =0; i < config_overlay_vector->size(); i++) {
				//cerr << "got: " << config_overlay_vector->at(i) << endl;
				config_element = parse_column_data_line(config_overlay_vector->at(i), '=');
				if (config_element->size()==2) {
					cerr << "additional command line option: " << config_element->at(0) << " = " << config_element->at(1) << endl;
					map_it = parameters->find(config_element->at(0));
					if (map_it == parameters->end()) {
						cerr << "error: parameter \"" << config_element->at(0) << "\" not known." << endl;
						exit(1);
					}
					(*parameters)[config_element->at(0)] = config_element->at(1);

				} else {
					cerr << "error: config_overlay parsing: " << config_overlay_vector->at(i)  << endl;
					exit(1);
				}
				delete config_element;
				config_element = 0;
		}
	}
	
	this->setBaseVariables();

	this->path = new char[this->max_blast_description_length];



	nucleobase_ascii2num = new int [256];
	for (int i=0; i< 256; i++) {nucleobase_ascii2num[i] = -1;}
	nucleobase_ascii2num['A']=0;
	nucleobase_ascii2num['T']=1;
	nucleobase_ascii2num['G']=2;
	nucleobase_ascii2num['C']=3;
	nucleobase_ascii2num['a']=0;
	nucleobase_ascii2num['t']=1;
	nucleobase_ascii2num['g']=2;
	nucleobase_ascii2num['c']=3;
	nucleobase_ascii2num['U']=1;
	nucleobase_ascii2num['u']=1;

	fastacmd_class = new Fastacmd_Class();
}




ReciprocResults * CARMA_BASE::doReciprocalStuffOnCompleteBlastDB(MatchingQuery * matching_query) {

	vector<Match_hit * >::iterator hit_it;
	vector<Match_hit * > * hits=matching_query->hits;

	char besthit_file_name [L_tmpnam];
	char reciprok_blast_results_file_name [L_tmpnam];

	char * tmpnam_res3 = tmpnam ( besthit_file_name );
	char * tmpnam_res4 = tmpnam ( reciprok_blast_results_file_name );

	int match_type = matching_query->match_type;

	// ------------------------------
	ofstream file_besthit (besthit_file_name);

	Match_hit * match_hit = getFirstUnmaskedHit(matching_query->hits);
	Blast_HSP * blast_hsp = match_hit->HSPs->front();
	tax_id best_hit_taxid = match_hit->ncbi_tax_id;

	file_besthit << "> " << match_hit->ncbi_tax_id << endl;
	file_besthit << *(blast_hsp->subject_sequence) << endl;
//cout << *(blast_hsp->subject_sequence) << endl;
	file_besthit.flush();
	file_besthit.close();
	// ------------------------------


	// We blasted the fragment against NR(or NT), now we take the best hit and blast it again against NR(or NT):

	if (match_type == MATCH_TYPE_BLASTX || match_type == MATCH_TYPE_BLASTX) {

		//string blastp_command  = "blastall -p blastp -d ";
		string blastp_command = *blastall_script;
		blastp_command.append(" -p blastp -d ");

		blastp_command.append(*blast_nr_database);
		blastp_command.append(" -F F -C 0 -i "); // -F F : turn off filtering - can filter away complete query and then crashes...
												// -C 0 : no composition-based statistics, scores higher identity lower, stupid...
		blastp_command.append(besthit_file_name);
		blastp_command.append(" -o ");
		blastp_command.append(reciprok_blast_results_file_name);

		if (PRINT_EXTERNAL_CALLS == 1) {
			cerr << "call: " << blastp_command << endl;
		}
		int ret_blastp = system(blastp_command.c_str());
		//int ret_blastp = system("blastall -p blastp -d /ramdisk/small_blast_db -i /ramdisk/best.fas -o /ramdisk/test.blastp");
		if (ret_blastp != 0) {
			cerr << "error: blastall -p blastp ..." << endl;
			cerr << matching_query->query << endl;

//			string copy_command = "cp ";
//			copy_command.append(hits_fasta_file_name);
//			copy_command.append(" ");
//			copy_command.append(besthit_file_name);
//			copy_command.append(" /vol/cluster-data-nocache/wgerlach/error/");
//			int system_res = system(copy_command.c_str());
//
//			cerr << copy_command << endl;

			exit(EXIT_FAILURE);
		}
	} else if (match_type == MATCH_TYPE_BLASTN) {
		//string blastn_command  = "blastall -p blastn -d ";
		string blastn_command = *blastall_script;
		blastn_command.append(" -p blastn -d ");

		blastn_command.append(*blast_nt_database);
		blastn_command.append(" -F F -C 0 -i "); // -F F : turn off filtering - can filter away complete query and then crashes...
												// -C 0 : no composition-based statistics, scores higher identity lower, stupid...
		blastn_command.append(besthit_file_name);
		blastn_command.append(" -o ");
		blastn_command.append(reciprok_blast_results_file_name);

		//cout << blastn_command << endl;
		if (PRINT_EXTERNAL_CALLS == 1) {
			cerr << "call: " << blastn_command << endl;
		}
		int ret_blastn = system(blastn_command.c_str());
		//int ret_blastp = system("blastall -p blastp -d /ramdisk/small_blast_db -i /ramdisk/best.fas -o /ramdisk/test.blastp");
		if (ret_blastn != 0) {
			cerr << "error: blastall -p blastp ..." << endl;
			cerr << matching_query->query << endl;

//			string copy_command = "cp ";
//			copy_command.append(hits_fasta_file_name);
//			copy_command.append(" ");
//			copy_command.append(besthit_file_name);
//			copy_command.append(" /vol/cluster-data-nocache/wgerlach/error/");
//			int system_res = system(copy_command.c_str());
//
//			cerr << copy_command << endl;

			exit(EXIT_FAILURE);
		}
	} else {
			cerr << "unknown BLAST type" << endl;
			exit(1);
	}

	// parse result file:
	ifstream reciprocal_result_stream (reciprok_blast_results_file_name);
	if (! reciprocal_result_stream.is_open()) {
		cerr << "error: file not open: " << reciprok_blast_results_file_name << endl;
		exit(EXIT_FAILURE);
	}
	MatchingQuery * reciprocal_results = parseBlastEntry(reciprocal_result_stream, 10,  NameToNCBI, true, false);
	reciprocal_result_stream.close();



	deleteFile(besthit_file_name);
	deleteFile(reciprok_blast_results_file_name);

	if (reciprocal_results == 0 ) {
		cerr << "reciprocal_results == 0 " << endl;
		exit(EXIT_FAILURE);
	}

	if (reciprocal_results->hits == 0 ) {
		cerr << "reciprocal_results->hits == 0 " << endl;
		exit(EXIT_FAILURE);
	}

	if (reciprocal_results->hits->size() == 0 ) {
		cerr << "reciprocal_results->hits->size() == 0 " << endl;
		exit(EXIT_FAILURE);
	}



// ------------------------------------------------------------------------------------------------------
// copy content of reciprocal_results to recRes:


	ReciprocResults * recRes = new ReciprocResults();
	recRes->match_type = matching_query->match_type;
	recRes->other_hits = new vector<pair<double, rank_t> * >();
	recRes->best_hit_taxid = match_hit->ncbi_tax_id;
	recRes->self_bitscore = 0;
	//recRes->self_bitscore = blast_hsp->bitscore;
	//recRes->query_bitscore = blast_hsp->bitscore;


	tax_id * best_hit_taxonomy=0;
	tax_id * other_hit_taxonomy=0;


	getTaxa(match_hit->ncbi_tax_id, &best_hit_taxonomy, this->parent_taxid, NCBI_MAX);



	for (hit_it = reciprocal_results->hits->begin(); hit_it != reciprocal_results->hits->end(); hit_it++) {
			Match_hit * reciprocal_hit = *hit_it;
			Blast_HSP * reciprocal_hsp = reciprocal_hit->HSPs->front();


			if (reciprocal_hsp->bitscore > recRes->self_bitscore) {
				recRes->self_bitscore = reciprocal_hsp->bitscore;
			}
			//cout << "bitscore: " << reciprocal_hsp->bitscore << endl;
			//cout << "tax_id: " << (int) reciprocal_hit->ncbi_tax_id << endl;

			getTaxa(reciprocal_hit->ncbi_tax_id, &other_hit_taxonomy, this->parent_taxid, NCBI_MAX);

			tax_id lct = getLowestCommonTaxID(best_hit_taxonomy, other_hit_taxonomy, NCBI_MAX);

			//cout << "last common tax_id: " << lct << endl;
			rank_t lct_rank = getLowestRank(lct, parent_taxid, ranks, NCBI_MAX);
			//cout << "rank: " << (int) lct_rank << endl;

			pair<double, rank_t> * ppp = new pair<double, rank_t>(reciprocal_hsp->bitscore, lct_rank);

			recRes->other_hits->push_back(ppp);

	}

	delete reciprocal_results;
// ------------------------------------------------------------------------------------------------------

	// We just obtained the distance of the best hit to all other hits, but the distance to the metagenomic read is still unknown:

	if (match_type == MATCH_TYPE_BLASTX) {
		// get GI of best hit so I can blast translated query gainst A in NR
		// get score A-Q


		// --- translated read:
		char translated_read_file_name [L_tmpnam];
		char * tmpnam_translated_read = tmpnam ( translated_read_file_name );
		//cout << "name: " << *match_hit->blast_id << endl;
		//cout << "translated read (query): " << *blast_hsp->query_sequence << endl;

		ofstream translated_read_file_stream(translated_read_file_name);
		translated_read_file_stream << "> translated read" << endl;
		translated_read_file_stream << *blast_hsp->query_sequence << endl;
		translated_read_file_stream.close();



		// --- gi list:
		char gilist_file_name [L_tmpnam];
		char * tmpnam_gilist = tmpnam ( gilist_file_name );

		int gi = extractGiFromBlastDB(match_hit->blast_id, NameToNCBI, match_type, blast_nr_database, fastacmd_script);

		if (gi == 0) {
			delete recRes;
			return 0;
		}

		//cout << "gi: " << gi << endl;
		ofstream gilist_file_stream (gilist_file_name);
		gilist_file_stream << gi << endl;
		gilist_file_stream.close();




		// --- result file
		char blastp_excl_result_file_name [L_tmpnam];
		char * tmpnam_blastp_excl_result = tmpnam ( blastp_excl_result_file_name );




		// -----------------------------
// example: blastall -p blastp -d /vol/cluster-data/wgerlach/NR25_genus/nr -F F -C 0 -i test.fas -l gi.list


		//string blastp_exclusive_command  = "blastall -p blastp -d ";
		string blastp_exclusive_command  = *blastall_script;
		blastp_exclusive_command.append(" -p blastp -d ");

		blastp_exclusive_command.append(*blast_nr_database);
		blastp_exclusive_command.append(" -F F -C 0 -i "); // -F F : turn off filtering - can filter away complete query and then crashes...
												// -C 0 : no composition-based statistics, scores higher identity lower, stupid...
		blastp_exclusive_command.append(translated_read_file_name);
		blastp_exclusive_command.append(" -o ");
		blastp_exclusive_command.append(blastp_excl_result_file_name);
		blastp_exclusive_command.append(" -l ");
		blastp_exclusive_command.append(gilist_file_name);

		if (PRINT_EXTERNAL_CALLS == 1) {
			cerr << "call: " << blastp_exclusive_command << endl;
		}
		int ret_blastp_ex = system(blastp_exclusive_command.c_str());
		//int ret_blastp = system("blastall -p blastp -d /ramdisk/small_blast_db -i /ramdisk/best.fas -o /ramdisk/test.blastp");
		if (ret_blastp_ex != 0) {
			cerr << "error: blastall -p blastp ..." << endl;
			cerr << matching_query->query << endl;

//			string copy_command = "cp ";
//			copy_command.append(hits_fasta_file_name);
//			copy_command.append(" ");
//			copy_command.append(besthit_file_name);
//			copy_command.append(" /vol/cluster-data-nocache/wgerlach/error/");
//			int system_res = system(copy_command.c_str());
//
//			cerr << copy_command << endl;

			exit(EXIT_FAILURE);
		}

		// parse result file:
		ifstream blastp_excl_result_stream (blastp_excl_result_file_name);
		if (! blastp_excl_result_stream.is_open()) {
			cerr << "error: file not open: " << blastp_excl_result_file_name << endl;
			exit(EXIT_FAILURE);
		}
		MatchingQuery * blastp_excl_results = parseBlastEntry(blastp_excl_result_stream, 10,  NameToNCBI, true, true);
		blastp_excl_result_stream.close();

		if (blastp_excl_results == 0 ) {
				cerr << "blastp_excl_results == 0 " << endl;
				exit(EXIT_FAILURE);
		}

		if (blastp_excl_results->hits == 0 ) {
				cerr << matching_query->query << endl;
				cerr << "blastp_excl_results->hits == 0 " << endl;
				return 0;
				exit(EXIT_FAILURE);
		}

		if (blastp_excl_results->hits->size() == 0 ) {
				cerr << "blastp_excl_results->hits->size() == 0 " << endl;
				exit(EXIT_FAILURE);
		}

		// we consider only first Hit/HSP (theres is only one hit, but theoretically more HSPs)
		Match_hit * blastp_excl_hit = getFirstUnmaskedHit(blastp_excl_results->hits);
		Blast_HSP * blastp_excl_hsp = blastp_excl_hit->HSPs->front();

		//cout << "_bitscore: " << blastp_excl_hsp->bitscore << endl;
		//cout << "_bitscore: " << (int) blastp_excl_hit->ncbi_tax_id << endl;
		recRes->query_bitscore = blastp_excl_hsp->bitscore;


		deleteFile(translated_read_file_name);
		deleteFile(gilist_file_name);
		deleteFile(blastp_excl_result_file_name);


	} else if (match_type == MATCH_TYPE_BLASTN) {

		recRes->query_bitscore = blast_hsp->bitscore;

		// do something !!!!!!!!
	} else {
			cerr << "not supported?!" << endl;
			exit(EXIT_FAILURE);
	}


	if (best_hit_taxonomy != 0) {
		delete [] best_hit_taxonomy;
	}

	if (other_hit_taxonomy != 0) {
		delete [] other_hit_taxonomy;
	}


	return recRes;
}





//enum Match_Type { BLASTX, BLASTN, BLASTP, TBLASTX, TBLASTN, PFAM_HMM};
string CARMA_BASE::MatchTypeToString(int mt) {
	string result_string = "unknown";
	switch (mt) {
		case MATCH_TYPE_BLASTP:
			result_string = "BLASTP";
			break;
		case MATCH_TYPE_BLASTX:
			result_string = "BLASTX";
			break;
		case MATCH_TYPE_BLASTN:
			result_string = "BLASTN";
			break;
		case MATCH_TYPE_PFAM_HMM:
			result_string = "PFAM_HMM";
			break;
		case MATCH_TYPE_RDP:
			result_string = "RDP";
			break;
	}

	return result_string;
}


ReciprocResults * CARMA_BASE::doReciprocalStuff(MatchingQuery * matching_query, int match_type) {


	vector<Match_hit * >::iterator hit_it;
	vector<Match_hit * > * hits=matching_query->hits;

	char hits_fasta_file_name [L_tmpnam];
	char small_blast_db_file_name [L_tmpnam];
	char besthit_file_name [L_tmpnam];
	char reciprok_blast_results_file_name [L_tmpnam];


	char * tmpnam_res1 = tmpnam ( hits_fasta_file_name );
	char * tmpnam_res2 = tmpnam ( small_blast_db_file_name );
	char * tmpnam_res3 = tmpnam ( besthit_file_name );
	char * tmpnam_res4 = tmpnam ( reciprok_blast_results_file_name );


	string rm_command = "rm -rf ";
	rm_command.append(hits_fasta_file_name);
	rm_command.append(" ");
	rm_command.append(small_blast_db_file_name);
	rm_command.append(".* ");
	rm_command.append(besthit_file_name);
	rm_command.append(" ");
	rm_command.append(reciprok_blast_results_file_name);


	ofstream hits_fasta_file_stream (hits_fasta_file_name);
	ofstream file_besthit (besthit_file_name);

	assert(matching_query);
	assert(matching_query->hits);

	Match_hit * best_hit = getFirstUnmaskedHit(matching_query->hits);

	if (best_hit == 0) {
		cerr << "error: best_hit == 0" << endl;
	}

	assert(best_hit->HSPs);

	Blast_HSP * best_hsp = best_hit->HSPs->front();





	hits_fasta_file_stream << "> query" << endl;
	if (best_hsp->query_sequence!= 0) {
	//assert(best_hsp->query_sequence);

		hits_fasta_file_stream << *(best_hsp->query_sequence) << endl;
	} else if (matching_query->query_sequence != 0) {
		hits_fasta_file_stream << *(matching_query->query_sequence) << endl;
	} else {
		cerr << "error: no query sequence found" << endl;
		exit(1);
	}

	//int frameshiftNum = matching_query->hits->front()->HSPs->front()->query_frameshifts;

//cout << "> query" << endl;
//cout << *(best_hit->HSPs->front()->query_sequence) << endl;


//cout <<  *(matching_query->hits->front()->HSPs->front()->query_sequence)  << endl;

	if (false) {
		tax_id one = 1;
		PhyloNode * tree = new PhyloNode(one);
	}
//				PhyloNode * treePointer;

	int positives_total=0;


	tax_id best_hit_taxid = 0;



	int best_query_start;
	int best_query_end;



	// iterate over all hits of the metagenomic read.

	int written_sequences = 1; // 1 for query
	bool found_best=false;
	bool is_best=false;
	int i=-1;
	for (hit_it = hits->begin(); hit_it != hits->end(); hit_it++){
		i++;
		Match_hit * match_hit = *hit_it;

		if (match_hit->tail) {
			break;
		}

		if (match_hit->masked) {
			continue;
		}

		Blast_HSP * blast_hsp = match_hit->HSPs->front();

		while (1) {
			//positives_total += match_hit->positives;

			//cout << match_hit->bitscore << endl;
			//cout << "blast_id: " << *(match_hit->blast_id) << endl;

			//FILE *fp;
			int status;

			if (found_best==false) {
				is_best=true;
				found_best=true;
			}


			if (is_best) {
				best_hit_taxid = match_hit->ncbi_tax_id;
				file_besthit << "> " << match_hit->ncbi_tax_id << endl;
				file_besthit << *(blast_hsp->subject_sequence) << endl;
				file_besthit.flush();

//cout << "This is A, the new query:" << endl;
//cout << "> " << match_hit->ncbi_tax_id << endl;
//cout << *(blast_hsp->subject_sequence) << endl;



				//cout << "BEST: " << endl;
				best_query_start = blast_hsp->query_start;
				best_query_end = blast_hsp->query_end;
			}

			string extracted_sequence;
			if (match_type == MATCH_TYPE_BLASTX) {
				extracted_sequence = extractSequenceFromBlastDB(match_hit, blast_hsp, match_type, best_query_start, best_query_end, blast_nr_database, blast_nt_database, path, max_blast_description_length, this->fastacmd_script);
			} else if (match_type == MATCH_TYPE_BLASTN || match_type == MATCH_TYPE_BLASTP) {
				extracted_sequence = *(blast_hsp->subject_sequence);
			} else {
				cerr << "error: match_type"<< endl;
				exit(1);
			}


//cout  <<  "got: " << extracted_sequence << endl;

//cout  <<  "q: " << matching_query->query << endl;
			//if (extracted_sequence.compare("") != 0) {
			if (extracted_sequence.length() >= 10) {
				written_sequences++;

				if (is_best) {
					hits_fasta_file_stream << "> self" << endl;
					//cout << "> self" << endl;
				} else {
					hits_fasta_file_stream << "> " << match_hit->ncbi_tax_id << endl;
					//cout << "> " << match_hit->ncbi_tax_id << endl;
				}
				hits_fasta_file_stream << extracted_sequence << endl;
				//cout << "extracted_sequence: " << extracted_sequence << endl;
			}
//exit(1);
			hits_fasta_file_stream.flush();
//cout << "f4" << endl;
//cerr << "call close... " << endl;
//fflush(fp);





//cerr << "feof: " <<  feof(fp) << endl;
//cerr << "ferror: " <<  ferror(fp) << endl;




//cerr << "call closed. " << endl;
			is_best=false;
			break;
		} //end while (1)
	} // end iterate over hits
//cout << "g" << endl;
	hits_fasta_file_stream.close();
	file_besthit.close();


//cout << "#written_sequences: " << written_sequences << endl;

	//formatdb  -i test.fas -p T -n small_blast_db
	//blastall -p blastp -d ./small_blast_db -i test.fas

	//system("rm -f small_blast_db.*");
//system("cat test.fas");
	string formatdb_command = *formatdb_script;
	formatdb_command.append(" -l /dev/null -i ");
	formatdb_command.append(hits_fasta_file_name);
	if (match_type == MATCH_TYPE_BLASTX || match_type == MATCH_TYPE_BLASTP) {
		formatdb_command.append(" -p T -n ");
	} else if (match_type == MATCH_TYPE_BLASTN) {
		formatdb_command.append(" -p F -n ");
	} else {
			cerr << "unknown BLAST type" << endl;
			exit(1);
	}
	formatdb_command.append(small_blast_db_file_name);

	//int ret_formatdb = system("formatdb  -i /ramdisk/test.fas -p T -n /ramdisk/small_blast_db");
	if (PRINT_EXTERNAL_CALLS == 1) {
		cerr << "call: " << formatdb_command << endl;
	}
	int ret_formatdb = system(formatdb_command.c_str());

	if (ret_formatdb != 0) {
		cerr << "error: formatdb exited abnormally" << endl;
		cerr << "query: " << matching_query->query << endl;

		//string error_copy = "cp ";
		//error_copy.append(hits_fasta_file_name);
		//error_copy.append(" /homes/wgerlach/cluster/error/");
		//cerr << error_copy << endl;
		//int system_out = system(error_copy.c_str());

		exit(EXIT_FAILURE);
	}
//cerr << "a" << endl;
	if (match_type == MATCH_TYPE_BLASTX || match_type == MATCH_TYPE_BLASTP) {

		//string blastp_command  = "blastall -p blastp -d ";
		string blastp_command  = *blastall_script;
		blastp_command.append(" -p blastp -d ");

		blastp_command.append(small_blast_db_file_name);

		if (match_type == MATCH_TYPE_BLASTP) {
			blastp_command.append(" -m 9");
		}

		blastp_command.append(" -F F -C 0 -i "); // -F F : turn off filtering - can filter away complete query and then crashes...
												// -C 0 : no composition-based statistics, scores higher identity lower, stupid...
		blastp_command.append(besthit_file_name);
		blastp_command.append(" -o ");
		blastp_command.append(reciprok_blast_results_file_name);

		if (PRINT_EXTERNAL_CALLS == 1) {
			cerr << "call: " << blastp_command << endl;
		}
		int ret_blastp = system(blastp_command.c_str());
		//int ret_blastp = system("blastall -p blastp -d /ramdisk/small_blast_db -i /ramdisk/best.fas -o /ramdisk/test.blastp");
		if (ret_blastp != 0) {
			cerr << "error: blastall -p blastp ..." << endl;
			cerr << matching_query->query << endl;

//			string copy_command = "cp ";
//			copy_command.append(hits_fasta_file_name);
//			copy_command.append(" ");
//			copy_command.append(besthit_file_name);
//			copy_command.append(" /vol/cluster-data-nocache/wgerlach/error/");
//			int system_res = system(copy_command.c_str());
//
//			cerr << copy_command << endl;

			exit(EXIT_FAILURE);
		}
	} else if (match_type == MATCH_TYPE_BLASTN) {
		//string blastn_command  = "blastall -p blastn -d ";
		string blastn_command  = *blastall_script;
		blastn_command.append(" -p blastn -d ");

		blastn_command.append(small_blast_db_file_name);
		blastn_command.append(" -m 9 -b 300 -F F -C 0 -i "); // -F F : turn off filtering - can filter away complete query and then crashes...
												// -C 0 : no composition-based statistics, scores higher identity lower, stupid...
		blastn_command.append(besthit_file_name);
		blastn_command.append(" -o ");
		blastn_command.append(reciprok_blast_results_file_name);

		//cout << blastn_command << endl;
		if (PRINT_EXTERNAL_CALLS == 1) {
			cerr << "call: " << blastn_command << endl;
		}
		int ret_blastn = system(blastn_command.c_str());
		//int ret_blastp = system("blastall -p blastp -d /ramdisk/small_blast_db -i /ramdisk/best.fas -o /ramdisk/test.blastp");
		if (ret_blastn != 0) {
			cerr << "error: blastall -p blastp ..." << endl;
			cerr << matching_query->query << endl;

			//string copy_command = "cp ";
			//copy_command.append(hits_fasta_file_name);
			//copy_command.append(" ");
			//copy_command.append(besthit_file_name);
			//copy_command.append(" /vol/cluster-data-nocache/wgerlach/error/");
			//int system_res = system(copy_command.c_str());

			//cerr << copy_command << endl;

			exit(EXIT_FAILURE);
		}
	} else {
			cerr << "unknown BLAST type" << endl;
			exit(1);
	}
//cout << "g2" << endl;


	ReciprocResults * recRes;

	if (match_type == MATCH_TYPE_BLASTX) {
		recRes = parseReciprocBLASTP(best_hit->ncbi_tax_id, this->ranks, reciprok_blast_results_file_name, matching_query->match_type);
	} else if (match_type == MATCH_TYPE_BLASTN || match_type == MATCH_TYPE_BLASTP) {
		recRes = parseReciprocBLASTNP(best_hit->ncbi_tax_id, this->ranks, reciprok_blast_results_file_name, matching_query->match_type);
	} else {
		cerr << "match_type unknown" << endl;
		exit(1);
	}


//	if (recRes != 0 && recRes->other_hits != 0) {
//		cout << "#hits obtained in reciprocal search (without self & query): " << recRes->other_hits->size() << endl;
//	} else {
//		cout << "#hits obtained in reciprocal search (without self & query): zero" << endl;
//	}

	// TODO only for blastx !
	//if (frameshiftNum != 0) {
		//recRes->query_bitscore = recRes->query_bitscore - (frameshiftNum*1);
		//recRes->query_bitscore = recRes->query_bitscore - 1;
		//cout << "#frameshift penalty: " << frameshiftNum << endl;
	//}

	recRes->best_hit_taxid = best_hit_taxid;


	int system_res = system(rm_command.c_str());

	return recRes;
}


ClassificationResult * CARMA_BASE::coreClassificationAlgorithm(string query_name, string evalue_str, Blast_HSP * blast_hsp, ReciprocResults * recRes, bool doReduceTaxon){

	bool debug_mode = false;

	ClassificationResult * cr;




	//string evalue_str = "-1";

	//if (match_hit != 0) {
	//	evalue_str = double2str(match_hit->HSPs->front()->evalue);
	//}



	if ((recRes == 0) || (recRes->other_hits == 0) || (recRes->query_bitscore == -1) || (recRes->self_bitscore == -1)) {
			//Match_hit * match_hit = matching_query->hits->front();

//			cerr << "recRes->other_hits: " << recRes->other_hits << endl;
//			cerr << "recRes->query_bitscore: " << recRes->query_bitscore << endl;
//			cerr << "recRes->self_bitscore: " << recRes->self_bitscore << endl;


		//	if (do_reduce_taxon) {
			//	int last_valid_node = reduceTaxon(match_hit->HSPs->front(), buffer_taxid, this);
				//cout << matching_query->query << "\tw2\tgoterms\t" << buffer_taxid[last_valid_node] << endl;
				//string property = "no_reziprok";
				//cr = new ClassificationResult(query_name, property, buffer_taxid[last_valid_node], 0, evalue_str);
			//} else {
				//cout << matching_query->query << "\tw2\tgoterms\t" << match_hit->ncbi_tax_id << endl;
				return 0;
				string none_evalue = "-";
				string property = "no_reziprok";
				//cr = new ClassificationResult(query_name, property, match_hit->ncbi_tax_id, 0, evalue_str);

				cr = new ClassificationResult(query_name, property, 0, 0, none_evalue);


			//}

//			if (recRes != 0) {
//				delete recRes;
//			}
				return cr;
	}


	vector<pair<double, rank_t> * >:: iterator recIt;

		// if hits are too few, discard prediction! Reason: different query reading frames !?!
//			int numberOfGoodHits = 0;
//			for (recIt = recRes->other_hits->begin(); recIt != recRes->other_hits->end(); recIt++) {
//
//				double hit_score = (*recIt)->first;
//				if (hit_score > 35) {
//					numberOfGoodHits++;
//				}
//			}

		//cout << "#written_sequences: " << written_sequences << endl;
		//cout << "numberOfGoodHits: " << numberOfGoodHits << endl;

//			if (numberOfGoodHits <= written_sequences*0.33) {
//				string property = "bad_reziprok_reading_frame";
//				ClassificationResult * cr = new ClassificationResult(matching_query->query, property, 0, 0);
//				delete recRes;
//				int system_res = system(rm_command.c_str());
//
//				return cr;
//			}


//cout << "i" << endl;
//					double gradients[8];
	double best_scores[8];
	double worst_scores[8];

	for (int i = 0; i< 8; i++) {
		best_scores[i]=0;
		worst_scores[i]=0;
	}
	best_scores[0]=recRes->self_bitscore;
	worst_scores[0]=recRes->self_bitscore;

	if (recRes->self_bitscore <= 0) {
		cerr << "recRes->self_bitscore: " << recRes->self_bitscore << endl;
		exit(1);
	}

//					int simpleRankQuery = 0;

//					double gradient_sum = 0;
//					int gradients_count = 0;
	if (debug_mode) {
		cout << "#recRes->self_bitscore: " << recRes->self_bitscore << endl;
		cout << "#recRes->query_bitscore: " << recRes->query_bitscore << endl;
	}


	//double query_multiplier =  1.03;  // example: 1.032;
	if (this->query_multiplier != 1) {
		double mul_result = recRes->query_bitscore * this->query_multiplier;
		// increase query bitscore maximal to self bitscore:
		recRes->query_bitscore = (mul_result<recRes->self_bitscore)?mul_result:recRes->self_bitscore;
		if (debug_mode) {
			cout << "#multiplied query_bitscore: " << recRes->query_bitscore << endl;
		}
	}


	// -------------------------------------------------------------------------------------------------------------------------
	for (recIt = recRes->other_hits->begin(); recIt != recRes->other_hits->end(); recIt++) {
		double hit_score = (*recIt)->first;


		int simpleRank = mapRankToSimpleRank((*recIt)->second);

	//cout << "# "<< hit_score << " " << simpleRank << " " << (int)(*recIt)->second << endl;

		if (simpleRank > 0) {

			if (hit_score > best_scores[simpleRank]) {
				best_scores[simpleRank] = hit_score;
				//gradients[simpleRank] = gradient;
			}



			if (worst_scores[simpleRank] == 0 || hit_score < worst_scores[simpleRank]) {
				worst_scores[simpleRank] = hit_score;
			}


//							cout << endl;
//							cout << "hit_score: " << hit_score << endl;
//							cout << "recRes->self_bitscore: " << recRes->self_bitscore << endl;
//							cout << "recRes->query_bitscore: " << recRes->query_bitscore << endl;
//							cout << "simpleRank: " << simpleRank << endl;
//							cout << "gradient: " << gradient << endl;
		}



		//cout << "have: " << (*recIt)->first << endl;

	}
//exit(0);


//					for (int i = 0; i< 8; i++) {
//							cout << i << "  " << best_scores[i] << endl;
//					}
	if (debug_mode) {
					cout << "#";
					for (int i = 0; i< 8; i++) {
							cout << " " << best_scores[i];
					}
					cout << endl;

					cout << "#";
					for (int i = 0; i< 8; i++) {
							cout << " " << worst_scores[i];
					}
					cout << endl;
	}

	// make it decreasing:
// var1 works!
	double max_score=0;
	for (int i = 7; i> 0; i--) {
			if (best_scores[i] > max_score) {
				max_score = best_scores[i];
			} else if (best_scores[i] < max_score) {
				best_scores[i] = 0;
			}
	}

// variant2
//	int max_score=0;
//	for (int i = 7; i> 0; i--) {
//			if (best_scores[i] != 0) {
//				if (best_scores[i] > max_score) {
//					max_score = best_scores[i];
//				} else if (best_scores[i] < max_score) {
//					best_scores[i] = max_score;
//				}
//			}
//	}


// var3
//	int max_score=0;
//	for (int i = 7; i> 0; i--) {
//			if (best_scores[i] != 0) {
//				if (best_scores[i] > max_score) {
//					max_score = best_scores[i];
//				} else if (best_scores[i] < max_score) {
//					best_scores[i] = 0;
//					worst_scores[i] = 0;
//				}
//			}
//	}


	if (max_score > best_scores[0]) {
			cerr << "warning: max_score is bigger than the best self_score... (" << query_name << ")" << endl;
	}

	if (debug_mode) {
				cout << "#";
				for (int i = 0; i< 8; i++) {
						cout << " " << best_scores[i];
				}
				cout << endl;
	}

			//		for (int i = 0; i< 8; i++) {
			//			cout << i << "  " << best_scores[i] << endl;
			//		}
	// interpolate:
	//int max_score=0;
	for (int i = 1; i<7 ; i++) {
		if (best_scores[i] == 0) {
			//find bigger value on the left (in direction species):
			int left=i-1;
			while (best_scores[left] == 0 && left >= 0) {
					left--;
			}
			if (best_scores[left] == 0) {
					cerr << "best_scores[left] == 0" << endl;
					exit(EXIT_FAILURE);
			}

			//find bigger value on the right (in direction unknown/superkingdom):
			int right=i+1;
			while (best_scores[right] == 0 && right < 7) {
					right++;
			}
			if (best_scores[right] == 0 || i==right) {
					break;
			}

			if (best_scores[right]==best_scores[left]) {
				best_scores[i]=best_scores[right];
			} else {
				//cout << "right: " << right << " : " << best_scores[right]<< endl;
				double gradient = (best_scores[left]-best_scores[right]) / (double)(left-right);
				best_scores[i]=best_scores[left]+(i-left)*gradient;
			}

			// test:
			if (worst_scores[i]==0) {
				worst_scores[i]=best_scores[i];
			}
		}

	}

//cerr << "d" << endl;
//					if (match_type == BLASTN) {
//						cout << "BLASTN" << endl;
//					} else {
//						cout << "BLASTX" << endl;
//					}

//					for (int i = 0; i< 8; i++) {
//							cout << i << "  " << best_scores[i] << endl;
//					}

	if (debug_mode) {
				cout << "#";
				for (int i = 0; i< 8; i++) {
						cout << " " << best_scores[i];
				}
				cout << endl;
	}


	// search best fitting rank for metagenomic fragment, start with "unknown/superkingdom"
	int simpleRank = 7;
	while ( (recRes->query_bitscore > best_scores[simpleRank-1] ) && simpleRank > 0) {

		simpleRank--; // ok, try next taxonomic rank

		//if (simpleRank < 0) {
		//		cerr << recRes->query_bitscore << endl;
		//		cerr << "recRes->query_bitscore >= best_scores[simpleRank] and i<0 !?" << endl;
		//		exit(EXIT_FAILURE);
		//}
	}
	//simpleRank++;
//cerr << "simpleRank: " << simpleRank << endl;
	// try to improve prediction, by jumping in next interval:
	if (simpleRank > 0) {

		if ( recRes->query_bitscore > best_scores[simpleRank]  ) { // when we are better than the current interval.. (special case: unknown-level)

//				if (recRes->query_bitscore == best_scores[simpleRank-1]) {
//					simpleRank--;
//				} else if ((recRes->query_bitscore >= worst_scores[simpleRank-1]) && (worst_scores[simpleRank-1] > 0)) { // > 0 not needed anymore?
//					// if score is better than other at this rank and at the next rank we would be in the interval, take the next rank.
//					simpleRank--;
//				}


				//if ( recRes->query_bitscore >= worst_scores[simpleRank-1] && (worst_scores[simpleRank-1] > 0) ) { // > 0 not needed anymore?
				//cerr << "recRes->query_bitscore: " << recRes->query_bitscore << endl;
				//cerr << "worst_scores[simpleRank-1]: " << worst_scores[simpleRank-1] << endl;
				//cerr << "simpleRank: " << simpleRank << endl;
				if ( recRes->query_bitscore >= worst_scores[simpleRank-1] ) { // > 0 not needed anymore?
				//if ( recRes->query_bitscore >= worst_scores[simpleRank-1] ) { // > 0 not needed anymore?
					// if score is better than other at this rank and at the next rank we would be in the interval, take the next rank.
					simpleRank--;
				}
				//cerr << "simpleRank: " << simpleRank << endl;

		}

	}
//cerr << "simpleRank: " << simpleRank << endl;

	bool error_no_worse_hit = false;

	if (best_scores[simpleRank] == 0 && simpleRank < 7) {
			//simpleRank = 7;
//			simpleRank--; // to achive LCA effect

//			return 0;

			error_no_worse_hit = true;

			//string property = "no_limit";
			//string none_evalue = "-";
			//cr = new ClassificationResult(query_name, property, 0, 0, none_evalue);
			//cr->error_no_worse_hit = true;
			//return cr;

	}

	//cout << "simpleRank: " << simpleRank << endl;

	rank_t computed_rank;
	//int simpleRank = (int) ceil(rank_double);

	if (simpleRank <= 0) {
		computed_rank = 27; // species
		//cout << "SPECIES" << endl;
	} else if (simpleRank == 1) {
		computed_rank = 23;
	} else if (simpleRank == 2) {
		computed_rank = 19;
	} else if (simpleRank == 3) {
		computed_rank = 14;
	} else if (simpleRank == 4) {
		computed_rank = 8;
	} else if (simpleRank == 5) {
		computed_rank = 5;
	} else if (simpleRank == 6) {
		computed_rank = 1; // superkingdom
	} else {
		computed_rank = 0; // unknown
	}
//cerr << "e" << endl;
//						if (computed_rank > 0 ) {
//							computed_rank--;
//						}
	//cout << "computed_rank: " << (int) computed_rank << endl;
	//cout << "best_hit_taxid: " << (int) best_hit_taxid << endl;

//			if (matching_query->query.compare("134324_218547440") == 0) {
//
//				cerr << rm_command << endl;
//				exit(1);
//			}
//cerr << "recRes->best_hit_taxid: " << recRes->best_hit_taxid<< endl;
	getTaxa(recRes->best_hit_taxid, &buffer_taxid, this->parent_taxid, NCBI_MAX);
	int last_valid_node = reduceTaxonByRank(buffer_taxid, this->ranks, computed_rank, NCBI_MAX);
	//cout << "last_valid_node: "<< last_valid_node << endl;



	if (last_valid_node == -1) {
		//cout << matching_query->query << "\twolfgang\tgoterms\t" << 0 << endl;

		//cout << "no_valid: " << (int) best_hit_taxid << endl;
		//cout << "wished rank: " << (int) simpleRank << endl;
		string property = "no_valid_rank";
		//string none_evalue = "-";
		cr = new ClassificationResult(query_name, property, 0, 0, evalue_str);

	} else {

		tax_id classi_result;
		if (last_valid_node >= 0) {
			classi_result = buffer_taxid[last_valid_node];
		} else {
			classi_result = 0;
		}

		if (classi_result > NCBI_MAX) {
			cerr << "1 classi_result > NCBI_MAX" << endl;
			exit(1);
		}

		if (classi_result > 0 && error_no_worse_hit && (blast_hsp != 0)) {  // TODO blast_hit check necessary?

			if (debug_mode) {
				cout << "#previous: " << (int)classi_result << endl;
			}


			//assert(match_hit->HSPs);

			// make more conservative prediction:
			//Blast_HSP * blast_hsp = match_hit->HSPs->front();
			//getTaxa(classi_result, &buffer_taxid, this->parent_taxid, NCBI_MAX);

			last_valid_node = reduceTaxon(blast_hsp, buffer_taxid , this);

			if (last_valid_node >= 0) {
				classi_result = buffer_taxid[last_valid_node];
			} else {
				classi_result = 0;
			}

			if (classi_result > NCBI_MAX) {
					cerr << "2 classi_result > NCBI_MAX" << endl;
					exit(1);
			}
			if (debug_mode) {
				cout << "#later: " << (int)classi_result << endl;
			}
		}

		//cout << matching_query->query << "\tcarma_"<< ((match_type==BLASTN)?"blastn":"blastx") <<"\tgoterms\t" << buffer_taxid[last_valid_node]  << endl;
		int match_type = recRes->match_type;

		string property = "carma_";
		property.append( MatchTypeToString(match_type) );


		cr = new ClassificationResult(query_name, property, classi_result, simpleRank, evalue_str);
		if (error_no_worse_hit) {
				cr->error_no_worse_hit = true;
		}

		//if (debug_mode) {
				//cout << rm_command << endl;
				//cerr << "#blubb, debug mode, exit now..."<< endl;
				//exit(0);
		//}

		//int system_res = system(rm_command.c_str());
		//return cr;
	}


		//cout << "gradients_count: " << gradients_count << endl;
		//cout << "gradient_sum: " <<  gradient_sum << endl;

		//tax_id lca_result = computeLCA(hits, 0.1, basic);

		//cout << (*mq_it).second->query << "\tlca\tgoterms\t" <<  lca_result << endl;





	//delete recRes;
	return cr;
	//cout << "closest_bitscore " << closest_bitscore << endl;
	//cout << "rank_of_closest_bitscore " << (int)rank_of_closest_bitscore << endl;


}


// returns 0 if ok
ClassificationResult * CARMA_BASE::filterMatchingQuery(MatchingQuery * matching_query, enum Algo_Type algo_type) {
	//cerr << "_have: " <<  matching_query->hits->front()->ncbi_tax_id << endl;
	bool do_reduce_taxon_for_all_hits = false;
//cerr << "xxx1" << endl;
	if (matching_query == 0) {
		string property = "no_hits1";
		string none_evalue = "-";
		ClassificationResult * cr = new ClassificationResult(matching_query->query, property, 0, 0, none_evalue);

		return cr;
	}

//cerr << "xxx" << endl;

	vector<Match_hit * > * hits=matching_query->hits;

	if (hits == 0 || hits->size()==0) {
		string property = "no_hits";
		string none_evalue = "-";
		ClassificationResult * cr = new ClassificationResult(matching_query->query, property, 0, 0, none_evalue);

		return cr;
	}

	//assert(getFirstUnmaskedHit(hits)->HSPs);




	vector<Match_hit * >::iterator hit_it;
	vector<Blast_HSP * >::iterator hsp_it;

	for (hit_it = hits->begin(); hit_it != hits->end(); ++hit_it){
		Match_hit * match_hit = *hit_it;

		// reset masking:
		match_hit->masked=false;
		match_hit->tail=false;

		if (match_hit->HSPs == 0) {
			cerr << "no HPSs !" << endl;
			exit(EXIT_FAILURE);
		}

		sort( match_hit->HSPs->begin(), match_hit->HSPs->end(), compare_hsp_by_bitscore );

	}


//cerr << "a "<<  matching_query->hits->front()->masked << endl;

	//string evalue_str = double2str(getFirstUnmaskedHit(hits)->HSPs->front()->evalue);


	double dMax=numeric_limits<double>::max();

	double evalue;
	double bitscore;
	int alignment_length;

	int match_type = matching_query->match_type;
	if (match_type == MATCH_TYPE_BLASTX) {
		evalue = carma_blastx_evalue;
		bitscore = carma_blastx_bitscore;
		alignment_length = carma_blastx_alignment_length;
	} else if (match_type == MATCH_TYPE_BLASTN) {
		evalue = carma_blastp_evalue;
		bitscore = carma_blastp_bitscore;
		alignment_length = carma_blastp_alignment_length;
	} else if (match_type == MATCH_TYPE_BLASTN) {
		evalue = carma_blastn_evalue;
		bitscore = carma_blastn_bitscore;
		alignment_length = carma_blastn_alignment_length;
	} 	else {
		evalue = dMax;
		bitscore = 0;
		alignment_length = 0;
	}


	double my_top_percent = 1;

	switch (algo_type) {
		case CARMA3:
			my_top_percent = CARMA_top_percent;
			break;
		case LCA:
			my_top_percent == LCA_top_percent;
			break;
	}


	if (evalue != dMax) {
		filterByEvalue(matching_query, evalue);
	}

	if (bitscore != 0) {
		filterByBitscore(matching_query, bitscore);
	}

	if (alignment_length != 0) {
		filterByAlignmentLength(matching_query, alignment_length);
	}


	if (my_top_percent != 1) {
		filterTopPercent(matching_query, my_top_percent);
	}


	filterOtherAndUnclassified(matching_query);



	if (do_reduce_taxon_for_all_hits) {
		filterReduceTaxa(matching_query);
	}



	int unmasked=0;
	for (hit_it = hits->begin(); hit_it != hits->end(); ++hit_it){
			Match_hit * match_hit = *hit_it;


			if (match_hit->tail)  {
				break;
			}

			if (match_hit->masked==false) {
				unmasked++;
			}
	}

	if (hits == 0 || unmasked == 0) {
		string property = "no_hits2";
		string none_evalue = "-";
		ClassificationResult * cr = new ClassificationResult(matching_query->query, property, 0, 0, none_evalue);

		return cr;
	}


	return 0;
}

// This is something SOrt-ITEMS does, CARMA does not do this. This was only for experimental purposes.
void CARMA_BASE::filterReduceTaxa(MatchingQuery * matching_query) {



	assert(matching_query);
	assert(matching_query->hits);

	vector<Match_hit * > * hits=matching_query->hits;



	vector<Match_hit * >::iterator hit_it;

	for (hit_it = hits->begin(); hit_it != hits->end(); ++hit_it){
		Match_hit * match_hit = *hit_it;


		if (match_hit->masked) {
			continue;
		}

		if (match_hit->HSPs == 0) {
			cerr << "no HPSs !" << endl;
			exit(EXIT_FAILURE);
		}




		getTaxa(match_hit->ncbi_tax_id, &buffer_taxid, this->parent_taxid, NCBI_MAX);
		int last_valid_node = reduceTaxon(match_hit->HSPs->front(), buffer_taxid, this);

		if (last_valid_node == -1) {
			match_hit->ncbi_tax_id = 0;

		} else {
			match_hit->ncbi_tax_id = buffer_taxid[last_valid_node];
		}


		if (match_hit->ncbi_tax_id > this->NCBI_MAX) {

			cerr << "buffer_taxid[last_valid_node] > this->NCBI_MAX" << endl;
			cerr << "ncbi_tax_id:" << (int)match_hit->ncbi_tax_id << endl;
			cerr << "last_valid_node:" << last_valid_node << endl;
			exit(1);
		}


	}


}

void CARMA_BASE::filterByEvalue(MatchingQuery * matching_query, double evalue) {



	assert(matching_query);
	assert(matching_query->hits);

	vector<Match_hit * > * hits=matching_query->hits;



	// --------------------------------------------------------------------------------------------------
	// filter all hits above global evalue

	vector<Match_hit * >::iterator hit_it;

	for (hit_it = hits->begin(); hit_it != hits->end(); ++hit_it){
		Match_hit * match_hit = *hit_it;

		if (match_hit->HSPs == 0) {
			cerr << "no HPSs !" << endl;
			exit(EXIT_FAILURE);
		}

		if (match_hit->HSPs->front()->evalue > evalue) {
			//delete match_hit;
			match_hit->masked=true;
			match_hit->tail=true;
		}


	}





}



void CARMA_BASE::filterByBitscore(MatchingQuery * matching_query, double bitscore) {

	assert(matching_query);
	assert(matching_query->hits);

	vector<Match_hit * > * hits=matching_query->hits;



	// --------------------------------------------------------------------------------------------------
	// filter all hits below global bitscore

	vector<Match_hit * >::iterator hit_it;

	for (hit_it = hits->begin(); hit_it != hits->end();  ++hit_it){
		Match_hit * match_hit = *hit_it;

		assert(match_hit->HSPs);
		if (match_hit->HSPs->size() == 0) {
			cerr << "match_hit->HSPs->size() == 0"  << endl;
			exit(1);
		}

		if (match_hit->HSPs->front()->bitscore < bitscore) {
			//delete match_hit;
			match_hit->masked=true;
			match_hit->tail=true;
		}


	}


}



void CARMA_BASE::filterByAlignmentLength(MatchingQuery * matching_query, int alignment_length) {

	assert(matching_query);
	assert(matching_query->hits);

	vector<Match_hit * > * hits=matching_query->hits;

	// --------------------------------------------------------------------------------------------------
	// filter all hits with short aligenment length

	vector<Match_hit * >::iterator hit_it;

	for (hit_it = hits->begin(); hit_it != hits->end(); ++hit_it){
		Match_hit * match_hit = *hit_it;

		if (match_hit->HSPs == 0) {
			cerr << "no HPSs !" << endl;
			exit(EXIT_FAILURE);
		}

		if (match_hit->HSPs->front()->alignment_length < alignment_length) {
			//delete match_hit;
			match_hit->masked=true;
		}


	}



}

void CARMA_BASE::filterTopPercent(MatchingQuery * matching_query, double _top_percent) {

	assert(matching_query);
	assert(matching_query->hits);

	vector<Match_hit * > * hits=matching_query->hits;


	// --------------------------------------------------------------------------------------------------
	// filter Hits below percent threshold

	double top_percent_min_score = (1-_top_percent) * matching_query->best_hit_bitscore;
	vector<Match_hit * >::iterator hit_it;


	bool breaked=false;
	for (hit_it = hits->begin(); hit_it != hits->end();  ++hit_it) {


		Match_hit * match_hit = *hit_it;

		assert(match_hit->HSPs);

		if ( match_hit->HSPs->front()->bitscore <  top_percent_min_score ) {
			//delete match_hit;
			match_hit->masked=true;
			match_hit->tail=true;
			breaked = true;
			break;
		}


	}
	if (! breaked) {
		matching_query->top_percent_exhausted=true;
	}

}


void CARMA_BASE::filterOtherAndUnclassified(MatchingQuery * matching_query) {

	assert(matching_query);

	assert(matching_query->hits);



	vector<Match_hit * > * hits=matching_query->hits;



	vector<Match_hit * >::iterator hit_it;

	int archaea = 0;
	int bacteria = 0;
	int eukaryota = 0;
	int viroids = 0;
	int viruses = 0;
	int somethingelse = 0;


	for (hit_it = hits->begin(); hit_it != hits->end();  ++hit_it) {


		Match_hit * match_hit = *hit_it;

		if (match_hit->tail) {
			break;
		}

		if (match_hit->masked) {
			continue;
		}


		if (match_hit->ncbi_tax_id == 0) {
			match_hit->masked = true;
			continue;
		}

		assert(match_hit->HSPs);



	//	bool have_deleted_hit = false;

		getTaxa(match_hit->ncbi_tax_id, &buffer_taxid, this->parent_taxid, NCBI_MAX);

		if ( abs( (int) buffer_taxid[0] ) != 1) {

				match_hit->masked = true; // can happen for deleted nodes (delnode.dmp)
				continue;
				//cerr << "abs( (int) buffer_taxid[0] ) != 1)" << endl;
				//cerr << "buffer_taxid[0] = " << (int) buffer_taxid[0]<< endl;
				//cerr << "match_hit->ncbi_tax_id: " << (int) match_hit->ncbi_tax_id<< endl;
				//exit(1);
		}




//		switch ((int)buffer_taxid[2]) {
//			case 2157:
//				archaea++;
//				break;
//			case 2:
//				bacteria++;
//				break;
//			case 2759:
//				eukaryota++;
//				break;
//			case 12884:
//				viroids++;
//				break;
//			case 10239:
//				viruses++;
//				break;
//			default:
//				somethingelse++;
//				break;
//		}


		//# skip "Other" and "Unclassified"
		// doens't work: if (buffer_taxid[2]== 28384 || buffer_taxid[2] == 12908) {
		if (not (buffer_taxid[2]== 2157 || buffer_taxid[2] == 2 || buffer_taxid[2] == 2759 || buffer_taxid[2] == 12884 || buffer_taxid[2] == 10239)) { // superkingdoms
			//delete match_hit;
			match_hit->masked=true;
		}





	}


//cout << "#archaea: " << archaea << endl;
//cout << "#bacteria: " << bacteria << endl;
//cout << "#eukaryota: " << eukaryota << endl;
//cout << "#viroids: " << viroids << endl;
//cout << "#viruses: " << viruses << endl;
//cout << "#somethingelse: " << somethingelse << endl;

}

ClassificationResult * CARMA_BASE::computeClassification(MatchingQuery * matching_query, bool doReduceTaxon) {

//cout << "computeClassification start" << endl;

	vector<Match_hit * > * hits=matching_query->hits;

	int match_type = matching_query->match_type;

// ------------------------------------------------------------------------------------------------------------------
	ClassificationResult * cr;

	ReciprocResults * recRes = 0;


	// count usable hits, those that are not marked
	int unmasked_hits=0;
	int pos=0;
	while (1) {
		if (pos >= matching_query->hits->size()) {
			break;
		}
		if (matching_query->hits->at(pos)->tail) {
			break;
		}
		if (! matching_query->hits->at(pos)->masked) {
			unmasked_hits++;
		}
		if (unmasked_hits > 1) {
			break;
		}
		pos++;
	}


//cerr << "unmasked_hits: " << unmasked_hits << endl;
	if (unmasked_hits > 1) {


		// ------------------------------------------------------------------------------------------------------------------
		// ------------------------------------------------------------------------------------------------------------------

		recRes = doReciprocalStuff(matching_query, match_type);

		// ------------------------------------------------------------------------------------------------------------------
		// ------------------------------------------------------------------------------------------------------------------
		//ReciprocResults * recRes = doReciprocalStuffOnCompleteBlastDB(matching_query);


		Match_hit * match_hit = getFirstUnmaskedHit(matching_query->hits);
		assert(match_hit);

		assert(match_hit->HSPs);
		Blast_HSP * blast_hsp = match_hit->HSPs->front();
		assert(blast_hsp);

		string evalue_str = double2str(match_hit->HSPs->front()->evalue);

		// ------------------------------------------------------------------------------------------------------------------
		// ------------------------------------------------------------------------------------------------------------------
		cr = coreClassificationAlgorithm(matching_query->query, evalue_str, blast_hsp, recRes, doReduceTaxon);
		// ------------------------------------------------------------------------------------------------------------------
		// ------------------------------------------------------------------------------------------------------------------


//		if (cr->error_no_worse_hit == true) {
//
//		}


	} else if (unmasked_hits == 1){
		// single hit
		string property = "single_hit";

		Match_hit * match_hit = getFirstUnmaskedHit(matching_query->hits);
		assert(match_hit);

		string evalue_str = double2str(match_hit->HSPs->front()->evalue);


		tax_id ncbi_tax_id = match_hit->ncbi_tax_id;
		if (doReduceTaxon) {


			getTaxa(match_hit->ncbi_tax_id, &buffer_taxid, this->parent_taxid, NCBI_MAX);
			int last_valid_node = reduceTaxon(match_hit->HSPs->front(), buffer_taxid, this);

			if (last_valid_node == -1) {
				ncbi_tax_id = 0;

			} else {
				ncbi_tax_id = buffer_taxid[last_valid_node];
			}
		}


		cr = new ClassificationResult(matching_query->query, property, ncbi_tax_id, 0, evalue_str);
		cr->error_no_worse_hit = true;

	}

	if (cr == 0 ) {
		// do not use these, very likley to be from wrong superkingdom !!!
		string property = "not_possible";

		string none_evalue = "-";

		cr = new ClassificationResult(matching_query->query, property, 0, 0, none_evalue);



	}

	if (recRes != 0) {
			delete recRes;
	}
//cout << "computeClassification end" << endl;
	return cr;
// ------------------------------------------------------------------------------------------------------------------
//cerr << "f" << endl;


	cerr << "not expected to come here" << endl;
	exit(EXIT_FAILURE);

	//cout << (*mq_it).second->query << endl;
	//exit(0);
	if (match_type == MATCH_TYPE_BLASTN) {
		cerr << "blubb blastn"<< endl;
		//cerr << rm_command << endl;
		exit(1);
	}
//		if (match_type == BLASTN) {
//		cout << rm_command << endl;
//		exit(0);
//		}
		//string grep_command = "grep -c query ";
		//grep_command.append(hits_fasta_file_name);

		//system(grep_command.c_str());
//if ((*mq_it).second->query.compare("134324_218547440")==0) {
						//cerr << previous_line << endl;
						//cerr << ncbi_tax_id << endl;
//							exit(EXIT_FAILURE);
//					}
		//exit(0);

//cerr << "g" << endl;


	if (false) {
		PhyloNode * tree; // fake line

		if (tree->children!=0) {

			if (ALGORITHM_TYPE == LCA) { //LCA

				tax_id * lca = tree->getLCA();
				int node_pos = 0;
				while (lca[node_pos] != -1) {
					//cout << "lca["<< node_pos <<"] =" << lca[node_pos] << endl;
					node_pos++;

					if (node_pos >= MAX_TREE_DEPTH) {
						cerr << "error: MAX_TREE_DEPTH" << endl;
						exit(EXIT_FAILURE);
					}

				}

				tax_id classification_result;
				if (node_pos > 0) {
					//cout << "LCA: " << lca[node_pos-1] << endl;
					classification_result = lca[node_pos-1];
				} else {
					//cout << "no lca." << endl;
					classification_result = 0;
				}

				//if ( (*mq_it)->hits->size() < 6) {
//						cout << (*mq_it)->query << "\tpfam1\tgoterms\t" << classification_result << endl;
//					} else {
//						cout << (*mq_it)->query << "\tpfam2\tgoterms\t" << classification_result << endl;
//					}


				//cout << *((*mq_it).second->hits->front()->sequence) << endl;
				cout << matching_query->query << "\tpfam\tgoterms\t" << classification_result  << "\tLCA" << endl;



			   // if (lca[node_pos-1] > 1 ) exit(0);

				delete [] lca;

			}
			else if (ALGORITHM_TYPE == REPR) {
				int positives_total = 0; // fake line
				list<PhyloNode * > * good_f_measure_nodes = new list<PhyloNode * >;
				tree->computeFMeasure(positives_total, matching_query->best_length , good_f_measure_nodes);
				//cout << "size: " << good_f_measure_nodes->size() << endl;
				//cout << "positives_total: " << positives_total << endl;
				list<PhyloNode * >::iterator it;
//					if ( (*mq_it).second->query.compare("196493_218547440") == 0) {
//						for (it = good_f_measure_nodes->begin(); it != good_f_measure_nodes->end(); it++) {
//							cout << "--------------" << endl;
//								cout << "ncbi : "<< (*it)->ncbi_tax_id << endl;
//								cout << "f-m  : " << (*it)->f_measure << endl;
//								cout << "pos  : " <<(*it)->positives << endl;
//								cout << "non  : " <<(*it)->non_positives << endl;
//								cout << "b_len:"<<(*mq_it).second->best_length << endl;
//						}
//						exit(0);
//					}

				//cout << "classification: " << good_f_measure_nodes->front()->f_measure << endl;

				cout << matching_query->query << "\tpfam\tgoterms\t" << good_f_measure_nodes->front()->ncbi_tax_id << "\tREPR"<< endl;

				delete good_f_measure_nodes;
//exit(0);

			} else {
				cerr << "unknown algorithm mode" << endl;
				exit(EXIT_FAILURE);
			}

		} else {
			cout << matching_query->query << "\tpfam\tgoterms\t" << 0 << endl;
		}
		tree->setHitsZero();
		delete tree;
	}



}






// frame_offset being -3,-2,-1,1,2,3 for each of 6 frames...
string CARMA_BASE::translateDNAFragment(int frame_offset, string& dna_sequence){

	int a;
	int b;
	int c;
	char aminoacid;
	string translation = "";

	int dna_length = dna_sequence.length();

	if (frame_offset > 0) {
		for (int i = (frame_offset-1); i < dna_length-2; i+=3) {

			a = nucleobase_ascii2num[(uchar)dna_sequence[i]];
			b = nucleobase_ascii2num[(uchar)dna_sequence[i+1]];
			c = nucleobase_ascii2num[(uchar)dna_sequence[i+2]];

			if (a >= 0 && b>=0 && c>= 0) {
				int codon = 16*a + 4*b + c;
				aminoacid = codon2aminoacid[codon];
			} else {
				aminoacid = 'X';
			}


			translation.append(1, aminoacid);
		}

	} else {

		int i = dna_length-3;
		while (i % 3 != ((-frame_offset)-1)) {
			i--;
		}
		while (i >= 0) {

			a = nucleobase_ascii2num[(uchar)dna_sequence[i]];
			b = nucleobase_ascii2num[(uchar)dna_sequence[i+1]];
			c = nucleobase_ascii2num[(uchar)dna_sequence[i+2]];

			if (a >= 0 && b>=0 && c>= 0) {
				int codon = 16*a + 4*b + c;
				aminoacid = codon2RCaminoacid[codon];
			} else {
				aminoacid = 'X';
			}


			//translation_R1.append(1, codon2RCaminoacid[codon]);
			translation.append(1, aminoacid);
			i-=3;
		}


	}


	return translation;
}



map<string, string * > * CARMA_BASE::parseStockholm(const char * stockholm_file){
	map<string, string * > * alignment = new map<string, string * >;

	ifstream myfile (stockholm_file);
	if (myfile.is_open())
	{

		string line;

		string description = "";
		string prot_sequence= "";

		while (! myfile.eof() )	{

			std::getline(myfile,line );

			if (line.length() == 0) {
				continue;
			}

			if (line.compare(0,1,"#")==0) {
				continue;
			}

			if (line.compare(0,2,"//")==0) {
				break;
			}

			//cout << line << endl;

			size_t sep = line.find_first_of(' ');
			if (sep == string::npos) {
					cerr << "error parsing stockholm: "<< line << endl;
					exit(EXIT_FAILURE);
			}

			description = line.substr(0, sep);


			sep = line.find_first_not_of(' ', sep+1);
			if (sep == string::npos) {
					cerr << "error parsing stockholm: " << line << endl;
					exit(EXIT_FAILURE);
			}

			prot_sequence = line.substr(sep);

			stringToUpper(prot_sequence);
			stringDotToGap(prot_sequence);
			//cout << "_" << description << "_" << endl;
			//cout << "_" << prot_sequence << "_" << endl;

			map<string, string * >::iterator ali_it;
			ali_it = alignment->find(description);

			if (ali_it == alignment->end()) {
				//cout << "new" << endl;
					alignment->insert( pair<string, string * >(description, new string(prot_sequence)) ) ;
			} else {
				//cout << "old" << endl;
					string * old_seq = ali_it->second;
					old_seq->append(prot_sequence);

					//cout << *old_seq << endl;
					//exit(0);
			}


		}
		myfile.close();
	} else {
		cerr << "Error: Unable to open file " << stockholm_file << endl;
		exit(1);
	}

	return alignment;
}




int CARMA_BASE::computeBlosumScore(string * seq_A, string * seq_B, int start, int end) {

	if (this->blosum_matrix == 0) {
		//string * blosum_file = new string("BLOSUM62");

		this->blosum_matrix = new Class_Substitution_Matrix(this->blosum_file);
	}

	int len_a = seq_A->length();
	int len_b = seq_B->length();

	int pos = start;

//cout << "A: "<< *seq_A << endl;
//cout << "B: "<< *seq_B << endl << endl;
	int sum = 0;
	int overlap_length = 0;
	while (pos < len_a && pos < len_b && pos <= end) {
			char a = (*seq_A)[pos];
			char b = (*seq_B)[pos];

			//cout << a << " " << b << endl;
			//cout << "score: "<< this->blosum_matrix->getScore(a, b) << endl;
			sum += this->blosum_matrix->getScore(a, b);
			pos++;
			overlap_length++;
	}

	if (overlap_length < this->pairwise_blosum_minoverlap) {
			return 0;
	}

	//cout << "Sum: "<< sum<< endl;

	return sum;
}

void CARMA_BASE::setAlignmentParameters(Blast_HSP * fake_blast_hsp, string * seq_query, string * seq_best, int start, int end) {

	if (this->blosum_matrix == 0) {
		//string * blosum_file = new string("BLOSUM62");

		this->blosum_matrix = new Class_Substitution_Matrix(this->blosum_file);
	}

	int len_query = seq_query->length();
	int len_best = seq_best->length();




	int pos = start;

	// alignment parameters
	int positives = 0;
	int identities = 0;
	int gaps=0;
	int length = 0;

//cout << "A: "<< *seq_query << endl;
//cout << "B: "<< *seq_best << endl << endl;
	int sum = 0;
	int overlap_length = 0;
	while (pos < len_query && pos < len_best && pos <= end) {
			char a = (*seq_query)[pos];
			char b = (*seq_best)[pos];

			if ((a=='-') && (b=='-')) {
					pos++;
					continue;
			}

			length++;

			if ( this->blosum_matrix->getScore(a, b) > 0 ) {
				positives++;
			}

			if (a==b ) {
				identities++;
			}

			if ((a=='-') || (b=='-')) {
				gaps++;
			}


			pos++;
	}

	fake_blast_hsp->positives = positives;
	fake_blast_hsp->identities = identities;
	fake_blast_hsp->gaps = gaps;

	fake_blast_hsp->percent_identities = (int)((double)identities/(double)length * 100);

	fake_blast_hsp->query_start = 1; // exact positions are not important here, only length.
	fake_blast_hsp->query_end = length;
}


void CARMA_BASE::putQuestionmarks(string * sequence) {
	size_t pos = 0;

	while ((*sequence)[pos] == '-') {
		(*sequence)[pos] = '?';
		pos++;
	}

	pos = sequence->length()-1;

	//cout << (*sequence)[pos] << endl;

	while ((*sequence)[pos] == '-') {
		//cout << "X" << endl;
		(*sequence)[pos] = '?';
		pos--;
	}
}








void CARMA_BASE::printClassificationResult(ClassificationResult * cl) {
	string taxon = getTaxonStringByTaxId(cl->classification_tax_id, &taxid_to_name, this->parent_taxid, this->names_dmp);

	string goterms ="goterms" ;

	if (cl->property.compare(0,2, "PF") == 0) {

		string pfam_acc = cl->property.substr(0, 7);
		//cout << pfam_acc << endl;
	//	exit(0);
		if (this->pfam_to_GOid == 0) {
				this->pfam_to_GOid = init_PfamToGOid(pfamA_txt_file, gene_ontology_txt_file, this->zcat_bin);
		}
		map<string, string>::iterator it;
		it = pfam_to_GOid->find(pfam_acc);

		if (it == pfam_to_GOid->end()) {
				//cerr << "error: " << pfam_acc << " not found in pfam_to_GOid." << endl;
				//exit(1);
			goterms = "{}";
		} else {

		//cout << "got: " << it->second << endl;
			goterms = "{" + it->second + "}";
		}
		//cout << "hallo " << endl;
		//exit(0);
	}

//if (cl->error_no_worse_hit == true) {
//	cout << "#error_no_worse_hit" << endl;
//}
	cout << cl->query << "\t" << cl->property << "\t" << goterms << "\t" << cl->classification_tax_id << "\t" << taxon << "\t" << cl->evalue << endl;

}










map<string, string * > * CARMA_BASE::readFastaFileIntoMap(string * fasta_file){

	map<string, string * > * sequence_map = new map<string, string * >();
	pair<map<string,string * >::iterator,bool> ret;


	FASTA_Parser * fasta_parser = new FASTA_Parser(*fasta_file, true, zcat_bin);

	string description;
	string * sequence;

	while (fasta_parser->getNextDescriptionLine(description) ) {

		if (description.length() > 0) {
			//processDNAFragment(description, sequence);

			sequence=fasta_parser->getSequence();

			ret=sequence_map->insert( pair<string,string * >(description,sequence) );
			if (ret.second==false) {
				cerr << "Error: fasta description \"" << description << "\" is not unique." << endl;
				exit(1);
			}

		}

	}

	return sequence_map;
}





DP_Matrix::DP_Matrix(int dp_matrix_size) {
	this->dp_matrix_size = dp_matrix_size;


	const char * blosum_file = "BLOSUM62";
	this->score_matrix = readBLOSUM(blosum_file);

	indel_score = -4; // TODO get it from score file directly...
	max_diag = 10;

	dp_matrix = new int * [dp_matrix_size];
	for (int i = 0; i < dp_matrix_size; i++) {
		dp_matrix[i] = new int[dp_matrix_size];
	}

	for (int i = 0; i < dp_matrix_size; i++) {
		dp_matrix[i][0]=-4*i;
		dp_matrix[0][i]=-4*i;
	}

}



int DP_Matrix::compute(string * sequence_A, string * sequence_B) {

	int a_len = sequence_A->length();
	int b_len = sequence_B->length();

	if (a_len > this->dp_matrix_size || b_len > this->dp_matrix_size ) {
		cerr << a_len << endl;
		cerr << b_len << endl;
		cerr << this->dp_matrix_size << endl;
		cerr << "DP-matrix too small" <<endl;
		exit(EXIT_FAILURE);
	}

	int maximum;
	char a;
	char b;
	int indel_a;
	int indel_b;
	int j_start;
	int j_end;

	//for (int i = 1; i <= a_len; i++) {
	//	for (int j = 1; j <= b_len; j++) {

	for (int i = 1; i <= a_len; i++) {

		// j is limited into the diagonal
		j_start = max(1, i-this->max_diag);
		j_end = min(b_len, i+this->max_diag);

		for (int j = j_start; j <= j_end; j++) {
			int diag = abs(i-j);
			if (diag > this->max_diag) {
				continue;
			}

			a = sequence_A->at(i-1);
			b = sequence_B->at(j-1);

			int match = dp_matrix[i-1][j-1] + (int)this->score_matrix[a + 256* b];

//			if (i==j) {
//				cout << a << ", " << b << ":  " << (int)this->score_matrix[a + 256* b] << endl;
//			}


			if (diag < this->max_diag) {
				indel_a = dp_matrix[i-1][j] + this->indel_score;
				indel_b = dp_matrix[i][j-1] + this->indel_score;
				maximum = (indel_a>indel_b)?indel_a:indel_b;
			} else {
				// diag == this->max_diag

				if (i-j < 0) {
					maximum = dp_matrix[i][j-1] + this->indel_score;
				} else {
					maximum = dp_matrix[i-1][j] + this->indel_score;
				}
			}

			maximum=(maximum>match)?maximum:match;

			//if (i==j) {
				//cout << "screibe at i="<< i << " j=" << j  << ": " << maximum << endl;
			//}

			dp_matrix[i][j]= maximum;
		}
	}

//cout << "a_len" <<  a_len << endl;
//cout << "b_len" << b_len << endl;
//cout << "lese: " << dp_matrix[a_len][b_len] << endl;
	return dp_matrix[a_len][b_len];
}


short int * readBLOSUM(const char * blosum_file) {

	string line;

	int matrix_length = 256*256;
	short int * matrix = new short int [matrix_length];

	for (int i = 0 ; i < matrix_length; i++) {
			matrix[i] = SHRT_MAX;
	}

	char * row = 0 ;

	ifstream myfile (blosum_file);
	if (myfile.is_open())
	{
		while (! myfile.eof() )	{
			char * pch;
			getline (myfile,line, '\n');


//cout << line << endl;



			if (line.length() == 0) {
				continue;
			}

			if (line[0] == '#') {
				continue;
			}


			char buf[1024];
			strncpy(buf, line.c_str(), sizeof(buf) - 1);
			buf[sizeof(buf) - 1] = '\0';

			if (line[0] == ' ') {
				row = new char [256];

				pch = strtok (buf," ");
				int i = 0;
				while (pch != NULL)
				  {
					//cout << i  << ": " <<  pch[0] << endl;
					char c = pch[0];

					//if (c<0 || c > 256) {
					if (c<0) {
							cerr << "c invalid" << endl;
							exit(EXIT_FAILURE);
					}

					row[i]=c;
					//printf ("%s\n",pch);
					pch = strtok (NULL, " ");
					i++;
				  }


				continue;
			}

			if (row == 0) {
					cerr << "row=0" << endl;
					exit(EXIT_FAILURE);
			}

			// normal rows:
			pch = strtok (buf," ");

			char c = pch[0];

			//if (c<0 || c > 256) {
			if (c<0) {
				cerr << "c invalid" << endl;
				exit(EXIT_FAILURE);
			}

			pch = strtok (NULL, " ");

			int i = 0;
			while (pch != NULL) {
				string score_str(pch);
				short int score_int = (short int) str2int(score_str);

				//cout << c << "," << row[i] << ": " << score_int << endl;
				matrix[c+256*row[i]]=score_int;
				matrix[c+32+256*row[i]]=score_int;
				matrix[c+256*(row[i]+32)]=score_int;
				matrix[c+32+256*(row[i]+32)]=score_int;
				pch = strtok (NULL, " ");
				i++;
			}



		}

		myfile.close();

		if (row != 0) delete [] row;

	}
	else {
		cerr << "Error: Unable to open file " << blosum_file << endl;
		exit(1);
	}

	return matrix;
}





string CARMA_BASE::getConfigValue(std::map<std::string, std::string> * parameters, const char * label){
	string label_s = label;
	std::map<std::string, std::string>::iterator map_it;
	map_it = parameters->find(label_s);
	if (map_it == parameters->end()) {
		cerr << "error: value for " << label_s << " not found in configuration file." << endl;
		exit(1);
	}

	string value = map_it->second;
	return value;
}

string * CARMA_BASE::getConfigStringValue(std::map<std::string, std::string> * parameters, const char * label){
	string value = getConfigValue(parameters, label);

	return new string(value);
}



double CARMA_BASE::getConfigDoubleValue(std::map<std::string, std::string> * parameters, const char * label){
	string value = getConfigValue(parameters, label);

	double res = strtod(value.c_str(), NULL);

	return res;
}

int CARMA_BASE::getConfigIntValue(std::map<std::string, std::string> * parameters, const char * label){
	string value = getConfigValue(parameters, label);
	int num = str2int(value);

	return num;
}

bool CARMA_BASE::getConfigBoolValue(std::map<std::string, std::string> * parameters, const char * label){

	string value = getConfigValue(parameters, label);

	stringToLower(value);

	if (value.compare("true") == 0 ) {
		return true;
	} else if (value.compare("false") == 0 ) {
		return false;
	}


	cerr << "error: value " << value << ", but expected true/false." << endl;
	exit(1);
	return false;
}


void CARMA_BASE::setBaseVariables() {

	this->query_multiplier=getConfigDoubleValue(parameters, "query_multiplier");


	this->nodes_dmp = getConfigStringValue(parameters, "nodes_dmp");
	//checkFileWarning(this->nodes_dmp);
	this->merged_dmp = getConfigStringValue(parameters, "merged_dmp");
	//checkFileWarning(this->merged_dmp);
	this->names_dmp = getConfigStringValue(parameters, "names_dmp");
	//checkFileWarning(this->names_dmp);

	this->gzip_bin = getConfigStringValue(parameters, "gzip_bin");
	//checkFileWarning(this->gzip_bin);

	this->zcat_bin = getConfigStringValue(parameters, "zcat_bin");
	//checkFileWarning(this->zcat_bin);


	this->blast_nr_database = getConfigStringValue(parameters, "blast_nr_database");
	this->blast_nt_database = getConfigStringValue(parameters, "blast_nt_database");
	this->pfamId2TaxId_file = getConfigStringValue(parameters, "pfamId2TaxId_file");
	//checkFileWarning(this->pfamId2TaxId_file);
	this->blastall_script = getConfigStringValue(parameters, "blastall_script");
	//checkFileWarning(this->blastall_script);
	this->fastacmd_script = getConfigStringValue(parameters, "fastacmd_script");
	//checkFileWarning(this->fastacmd_script);
	this->formatdb_script = getConfigStringValue(parameters, "formatdb_script");
	//checkFileWarning(this->formatdb_script);

	this->hmmfetch_bin = getConfigStringValue(parameters, "hmmfetch_bin");
	this->hmmalign_bin = getConfigStringValue(parameters, "hmmalign_bin");
	this->hmmscan_bin = getConfigStringValue(parameters, "hmmscan_bin");
	this->pfam_A_hmm_file = getConfigStringValue(parameters, "pfam_A_hmm_file");
	//checkFileWarning(this->pfam_A_hmm_file);
	this->pfam_fasta_dir = getConfigStringValue(parameters, "pfam_fasta_dir");
	//checkFileWarning(this->pfam_fasta_dir);
	this->blosum_file = getConfigStringValue(parameters, "blosum_file");
	//checkFileWarning(this->blosum_file);

	this->pfamA_txt_file = getConfigStringValue(parameters, "pfamA_txt_file");
	//checkFileWarning(this->pfamA_txt_file);
	this->gene_ontology_txt_file = getConfigStringValue(parameters, "gene_ontology_txt_file");
	//checkFileWarning(this->gene_ontology_txt_file);


	this->blastx_evalue=getConfigDoubleValue(parameters, "blastx_evalue");
	this->blastn_evalue=getConfigDoubleValue(parameters, "blastn_evalue");
	this->blastp_evalue=getConfigDoubleValue(parameters, "blastp_evalue");

	this->carma_blastx_evalue=getConfigDoubleValue(parameters, "carma_blastx_evalue");
	this->carma_blastn_evalue=getConfigDoubleValue(parameters, "carma_blastn_evalue");
	this->carma_blastp_evalue=getConfigDoubleValue(parameters, "carma_blastp_evalue");

	this->carma_blastx_bitscore=getConfigDoubleValue(parameters, "carma_blastx_bitscore");
	this->carma_blastp_bitscore=getConfigDoubleValue(parameters, "carma_blastp_bitscore");
	this->carma_blastn_bitscore=getConfigDoubleValue(parameters, "carma_blastn_bitscore");

	this->carma_blastx_alignment_length=getConfigIntValue(parameters, "carma_blastx_alignment_length");
	this->carma_blastp_alignment_length=getConfigIntValue(parameters, "carma_blastp_alignment_length");
	this->carma_blastn_alignment_length=getConfigIntValue(parameters, "carma_blastn_alignment_length");

	//this->CARMA_top_percent = getConfigDoubleValue(cfg, "CARMA_top_percent");
	this->LCA_top_percent = getConfigDoubleValue(parameters, "LCA_top_percent");
	this->hmmscan_evalue = getConfigDoubleValue(parameters, "hmmscan_evalue");
	this->pairwise_blosum_minscore = getConfigDoubleValue(parameters, "pairwise_blosum_minscore");
	this->pairwise_blosum_minoverlap = getConfigDoubleValue(parameters, "pairwise_blosum_minoverlap");

	this->use_hard_threshold = getConfigBoolValue(parameters, "use_hard_threshold");

	this->max_blast_description_length = getConfigIntValue(parameters, "max_blast_description_length");
	this->NCBI_MAX = getConfigIntValue(parameters, "NCBI_MAX");

	this->blast_rdp_database = getConfigStringValue(parameters, "blast_rdp_database");

	this->score_match = getConfigIntValue(parameters, "score_match");
	this->score_mismatch = getConfigIntValue(parameters, "score_mismatch");
	this->score_gapopen = getConfigIntValue(parameters, "score_gapopen");
	this->score_gapextension = getConfigIntValue(parameters, "score_gapextension");

}

void CARMA_BASE::readConfiguration(string * config_file){


	namespace pod = boost::program_options::detail;
	std::ifstream config(config_file->c_str());
    if(!config)
    {
        cerr << "I/O error while reading file." << endl;
		cerr << "Configuration file not found. retry in 20 sec.." << endl;
		sleep(20);

		config.open(config_file->c_str());
		if(!config)
		{
			cerr << "I/O error while reading file." << endl;
			cerr << "Configuration file not found. exit now." << endl;

			exit(EXIT_FAILURE);
		}
    }

    //parameters
    std::set<std::string> options;
    //std::map<std::string, std::string> parameters;
    parameters = new std::map<std::string, std::string>;
    options.insert("*");

    try
    {
        for (pod::config_file_iterator i(config, options), e ; i != e; ++i)
        {
            //std::cout << i->string_key <<" "<<i->value[0] << std::endl;
            (*parameters)[i->string_key] = i->value[0];
        }
        //std::cout << (*parameters)["StatLogServer.Path"] << std::endl;
    }
    catch(std::exception& e)
    {
        std::cerr<<"Exception: "<<e.what()<<std::endl;
    }


}





