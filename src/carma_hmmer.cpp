
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


#include "carma_hmmer.hpp"


void CARMA_HMMER::processFastaFile(string * dna_fasta_file, string * blast_egts, bool is_dna){


	if (blast_egts != 0 && is_dna) {
		this->blast_egts=blast_egts;
		blast_egts_map = readFastaFileIntoMap(this->blast_egts);
	}


	FASTA_Parser * fasta_parser = new FASTA_Parser(*dna_fasta_file, true, zcat_bin);

	string descr;
	string * sequence;
	string my_sequence;
	while (fasta_parser->getNextDescriptionLine(descr)) {

//cerr << "descr: " << descr<< endl;
		sequence = fasta_parser->getSequence();
		my_sequence = *sequence;

		delete sequence;

		if (descr.length() > 0 && my_sequence.length() >= 10) {

			processDNAFragment(descr, my_sequence, is_dna);

		}

	}

	delete fasta_parser;

}

void CARMA_HMMER::processDNAFragment(string& descr, string& sequence, bool is_dna){


	string translation_1;
	string translation_2;
	string translation_3;

	string translation_R1;
	string translation_R2;
	string translation_R3;
//cerr << "huhu:" << is_dna << endl;

	if (is_dna) {

		translation_1 = translateDNAFragment(1, sequence);
		translation_2 = translateDNAFragment(2, sequence);
		translation_3 = translateDNAFragment(3, sequence);

		translation_R1 = translateDNAFragment(-1, sequence);
		translation_R2 = translateDNAFragment(-2, sequence);
		translation_R3 = translateDNAFragment(-3, sequence);
	}

	string * translation_blast;
	bool found_blast_translation=false;

	if (this->blast_egts != 0 && is_dna) {


//			size_t start = descr.find_first_not_of(' ', 1);
//			if (start == string::npos) {
//				cerr << "error parsing: " << descr << endl;
//				exit(1);
//			}
//			string pure_descr = descr.substr(start);
			//cout << "pure_descr: " << pure_descr << endl;


			map<string, string * >:: iterator it;

			//string pure_descr = descr.substr(1);

			it = this->blast_egts_map->find(descr);
			if (it != this->blast_egts_map->end()) {
				found_blast_translation=true;
				//cout << "found: " << it->first << endl;
				//cout << "found: " << it->second << endl;
				translation_blast=it->second;
			}
	}

//exit(0);

	char hmmscan_input_file_name [L_tmpnam];


	char * tmpnam_res1 = tmpnam ( hmmscan_input_file_name );


	ofstream hmmscan_input_file_stream (hmmscan_input_file_name);

	if (is_dna) {
		hmmscan_input_file_stream << ">" << descr << "_1" << endl;
		hmmscan_input_file_stream << translation_1 << endl;
		hmmscan_input_file_stream << ">" << descr << "_2" << endl;
		hmmscan_input_file_stream << translation_2 << endl;
		hmmscan_input_file_stream << ">" << descr << "_3" << endl;
		hmmscan_input_file_stream << translation_3 << endl;
		hmmscan_input_file_stream << ">" << descr << "_-1" << endl;
		hmmscan_input_file_stream << translation_R1 << endl;
		hmmscan_input_file_stream << ">" << descr << "_-2" << endl;
		hmmscan_input_file_stream << translation_R2 << endl;
		hmmscan_input_file_stream << ">" << descr << "_-3" << endl;
		hmmscan_input_file_stream << translation_R3 << endl;

		if (found_blast_translation) {
			hmmscan_input_file_stream << ">" << descr << "_blastx" << endl;
			hmmscan_input_file_stream << *translation_blast << endl;
		}

	} else {
		// use protein
		hmmscan_input_file_stream << ">" << descr << "_1" << endl;
		hmmscan_input_file_stream << sequence << endl;
	}



	// example: hmmscan-3 /vol/biodb/pfam24/Pfam-A.hmm /tmp/fileDh36vp

	string hmmscan_output_file_name_string = hmmscan(this->hmmscan_bin, this->hmmscan_evalue, this->pfam_A_hmm_file, hmmscan_input_file_name);

	vector<HmmscanDomainMatch * > * hdms = parse_hmmscan_results(hmmscan_output_file_name_string.c_str());

//cerr << "hdms->size(): "<< hdms->size()<< endl;
	if (hdms->size() > 0) {
		//cout << "XXXX found something... :-)" << endl;

		// 1. search best matching family
		// 2. check for multiple frames
		//     yes -> try frameshift detection
		//     no -> all fine


		vector<HmmscanDomainMatch * >::iterator hdm_it;



		HmmscanDomainMatch * best_hdm = *hdms->begin();

		// search for best alignment
		for (hdm_it = hdms->begin(); hdm_it != hdms->end(); hdm_it++) {
			HmmscanDomainMatch * hdm = *hdm_it;
			//cout << hdm->full_seq_evalue << endl;
			if (hdm->full_seq_evalue < best_hdm->full_seq_evalue) {
				best_hdm = hdm;

			} else if (hdm->full_seq_evalue == best_hdm->full_seq_evalue) {

				if (hdm->i_evalue < best_hdm->i_evalue) {
					best_hdm = hdm;
				} else if (hdm->i_evalue < best_hdm->i_evalue) {
					if (best_hdm->frame == 0) {
						// if frameshift-free translation was as good as blastx-translation, take the frameshift-free one.
						best_hdm = hdm;
					}
				}
			}
		}



		string goterms;
		string pfam_acc = best_hdm->accession.substr(0, 7);
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

		string frame_string = (best_hdm->frame==0)?"X":stringify(best_hdm->frame);

		if (best_hdm->ali_from < best_hdm->ali_to) {

				//cout << ">" << best_hdm->accession << "=+=" << descr.substr(1) << "_" << ((best_hdm->frame >= 0)?1:-1) << "_" << best_hdm->frame << "=+=" << best_hdm->full_seq_evalue << "=+=" << goterms << endl;
				cout << ">" << best_hdm->accession << "=+=" << descr << "_" << ((best_hdm->frame >= 0)?1:-1) << "_" << frame_string << "=+=" << best_hdm->full_seq_evalue << "=+=" << goterms << endl;

				string egt;

				if (is_dna) {
					if (best_hdm->frame == 1) {
						egt = translation_1;
					} else if (best_hdm->frame == 2) {
						egt = translation_2;
					} else if (best_hdm->frame == 3) {
						egt = translation_3;
					} else if (best_hdm->frame == -1) {
						egt = translation_R1;
					} else if (best_hdm->frame == -2) {
						egt = translation_R2;
					} else if (best_hdm->frame == -3) {
						egt = translation_R3;
					} else { // == 0
						egt = *translation_blast;
					}
				} else {
					egt = sequence;
				}
				//cerr << "org: " << egt << endl;
				//cerr << best_hdm->ali_from << endl;
				//cerr << best_hdm->ali_to << endl;

				egt = egt.substr(best_hdm->ali_from-1, (best_hdm->ali_to - best_hdm->ali_from));
				//cout << best_hdm->ali_from << endl;
				//cout << best_hdm->ali_to << endl;
				cout << egt << endl;
		} else {

			cerr << "Error: best_hdm->ali_from >= best_hdm->ali_to , but continue..." << endl;
		}

	}


	DeleteContainerWithPointers(hdms);



	deleteFile(hmmscan_input_file_name);

	deleteFile(hmmscan_output_file_name_string.c_str());

	return;
}


vector<HmmscanDomainMatch * > * CARMA_HMMER::parse_hmmscan_results(const char * hmmscan_output_file_name){

//format of hmmscan ouput: --domtblout
//# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
//#------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
//Transpeptidase       PF00905.15   304 576.1_209693634_3    -             89   3.9e-11   42.1   0.0   1   1   7.3e-15   4.3e-11   42.0   0.0   155   199    14    58     8    73 0.90 Penicillin binding protein transpeptidase domain
//cerr << "1" << hmmscan_output_file_name << endl;
//exit(1);
	vector<HmmscanDomainMatch * > * hmm_matches = new vector<HmmscanDomainMatch * >;

	ifstream myfile (hmmscan_output_file_name);
	if (myfile.is_open())
	{

		string line;
		string query;
		while (! myfile.eof() )	{
			std::getline(myfile,line );


			if (line.length() == 0) {
				continue;
			}

			if (mycompare(line, 0, 1, "#") == 0) {
				continue;
			}

			//cerr << "got: " << line << endl;

			vector<string> * coldata = parse_column_data_line(line, ' ');


			//--- frame ---------------------------------------------

			string queryname = coldata->at(3);
			size_t left = queryname.find_last_of('_');

			if (left == string::npos) {
					cerr << "error: left == string::npos" << endl;
					cerr << queryname << endl;
					exit(EXIT_FAILURE);
			}
			left++;

			string frame_string =  queryname.substr(left);
			int frame;

			if (frame_string.compare(0, 6, "blastx") == 0) {

				frame = 0;

			} else {

				frame = str2int(frame_string);

				if ((frame < -3) || (frame > 3) || (frame == 0)) {
					cerr << "frame number "<< frame_string << " is invalid.." << endl;
					exit(EXIT_FAILURE);
				}
			}

			//cout << "frame: " << frame << endl;

			// --------------------------------------------------------


			HmmscanDomainMatch * hdm = new HmmscanDomainMatch();
			hdm->accession = coldata->at(1);
			hdm->frame = frame;
			hdm->full_seq_evalue = atof(coldata->at(6).c_str());
			hdm->i_evalue = atof(coldata->at(12).c_str());
			hdm->ali_from = str2int(coldata->at(17)); // query coordinates
			hdm->ali_to = str2int(coldata->at(18));

			if (hdm->ali_from < hdm->ali_to) {
				hmm_matches->push_back(hdm);
			} else {
				delete hdm;
			}
			//cout << "i_evalue string:" << coldata->at(12) << endl;
			delete coldata;



		}


	} else {
			cerr << "error reading file: " << hmmscan_output_file_name << endl;
			exit(EXIT_FAILURE);
	}

	return hmm_matches;
}


vector<string> * CARMA_HMMER::getEGTDescription(string& descr) {
	// parsing of EGT-description line
	// example: >PF04262.7=+=25.1_209693634_-1_-2=+=8.7e-06=+={}

	size_t sep1 = descr.find("=+=");

	if (sep1 == string::npos) {
		cerr << "parsing error, maybe not an EGT-file?" << endl;
		cerr << "line: "<< descr << endl;
		exit(EXIT_FAILURE);
	}


	size_t sep2 = descr.find("=+=", sep1+1);

	if (sep2 == string::npos) {
		cerr << "parsing error, maybe not an EGT-file?" << endl;
		exit(EXIT_FAILURE);
	}


	size_t sep3 = descr.find("=+=", sep2+1);

	if (sep3 == string::npos) {
		cerr << "parsing error, maybe not an EGT-file?" << endl;
		exit(EXIT_FAILURE);
	}



	string pfam_descr = descr.substr(1, sep1-1);
//	cout << pfam_descr << endl;

	string read_name = descr.substr(sep1+3, sep2-sep1-3);
//	cout << read_name << endl;

	string evalue_string = descr.substr(sep2+3, sep3-sep2-3);
//	cout << evalue_string << endl;

	string go_ids = descr.substr(sep3+3);
//	cout << go_ids << endl;

	vector<string> * res_vec = new vector<string>;

	res_vec->push_back(pfam_descr);
	res_vec->push_back(read_name);
	res_vec->push_back(evalue_string);
	res_vec->push_back(go_ids);

	if (res_vec->size() != 4) {
			cerr << "error: parsing EGT descr line..." << endl;
			exit(EXIT_FAILURE);
	}

	return res_vec;
}

void CARMA_HMMER::insertSequence(string& description, string& dna_sequence, map<string, vector<Sequence * > * > * sequencesByFamily) {
	Sequence * seq = new Sequence(description, dna_sequence);

	vector<string> * egt_descr_parts = getEGTDescription(description) ;
	string pfam_family = (*egt_descr_parts)[0];
	delete egt_descr_parts;

	//cout << pfam_family << endl;

	map<string, vector<Sequence * > * >::iterator family_it;

	family_it = sequencesByFamily->find(pfam_family);

	vector<Sequence * > * vec = 0;
	if (family_it == sequencesByFamily->end()) {
		vec = new vector<Sequence * >;
		sequencesByFamily->insert( pair<string, vector<Sequence * > * >(pfam_family, vec) );
	} else {
		vec = family_it->second;
	}

	vec->push_back(seq);


}



map<string, vector<Sequence * > * > * CARMA_HMMER::parseEGTFile(string * egt_file) {

	map<string, vector<Sequence * > * > * sequencesByFamily = new map<string, vector<Sequence * > * > ;
//cerr << "hallo welt. "<< endl;
	ifstream myfile (egt_file->c_str());
	if (myfile.is_open())
	{

		string line;

		string description = "";
		string dna_sequence= "";

		while (! myfile.eof() )	{

			std::getline(myfile,line );
//cerr << "read: " << line << endl;

			if (line.length() == 0) {
				continue;
			}

			if (line.compare(0,1,"#")==0) {
				continue;
			}

			if (line.compare(0,1,">")==0) {

				if (description.length() > 0 && dna_sequence.length() >= 10) {
					//processEGTFragment(description, dna_sequence);

					//cerr << "description: " << description << endl;
					//cerr << "dna_sequence: " << dna_sequence << endl;

					insertSequence(description, dna_sequence, sequencesByFamily);


					description = "";
					dna_sequence= "";
					//break;
				}
				description = line;

				continue;
			}

			dna_sequence.append(line);

		}

		if (description.length() > 0 && dna_sequence.length() >= 10) {
			//processEGTFragment(description, dna_sequence);
			insertSequence(description, dna_sequence, sequencesByFamily);
		}


	} else {
		cerr << "Error: Unable to open file " << *egt_file << endl;
		exit(1);
	}


	return sequencesByFamily;


}


// -------------------------------------------
// - fetch family hmm
// - get PFAM alignment
// - add EGTs to alignment
// - hmmalign
// - extract aligned EGTs
// - iterate through EGTs and classify
void CARMA_HMMER::processEGTfamily(string pfam_descr, vector<Sequence * > * vec) {


	// ---------------------------------------------------------------------------
	//fetch family hmm
	string single_hmm_file_name_string = hmmfetch(this->hmmfetch_bin, this->pfam_A_hmm_file, pfam_descr);





	// ---------------------------------------------------------------------------
	// create fasta file for multiple alignment:


	// 1. get PFAM alignment

	char combined_egt_pfam_file_name [L_tmpnam];
	char * tmpnam_res2 = tmpnam ( combined_egt_pfam_file_name );

	string copy_command =("cp ");
	// source
	copy_command.append(*pfam_fasta_dir);
	copy_command.append(pfam_descr);
	copy_command.append(".fas");
	copy_command.append(" ");
	// target
	copy_command.append(combined_egt_pfam_file_name);

	int sys_ret3 = system(copy_command.c_str());
	if (sys_ret3 != 0) {
			cerr << "error calling " << copy_command << endl;
			exit(EXIT_FAILURE);
	}
	if (PRINT_EXTERNAL_CALLS == 1) {
		cerr << "call: " << copy_command << endl;
	}


	// 2. add EGTs to alignment

	ofstream combined_egt_pfam_file_stream (combined_egt_pfam_file_name, ios_base::app);

	int count_egts = vec->size();

	for (int i = 0; i<count_egts; i++) {
		combined_egt_pfam_file_stream << ">EGT_" << i << endl;
		combined_egt_pfam_file_stream << vec->at(i)->data_sequence << endl;
	}

	combined_egt_pfam_file_stream.close();



	// ---------------------------------------------------------------------------
	// hmmalign

	string new_alignment_file_name_string = hmmalign(this->hmmalign_bin, single_hmm_file_name_string, combined_egt_pfam_file_name );





	deleteFile(single_hmm_file_name_string.c_str());
	deleteFile(combined_egt_pfam_file_name);

	// ---------------------------------------------------------------------------
	// ---------------------------------------------------------------------------
	// parse Stockholm for pfam and EGTs:

	map<string, string * >::iterator ali_it_pfam;
	map<string, string * > * pfam_alignment = parseStockholm(new_alignment_file_name_string.c_str());


	deleteFile(new_alignment_file_name_string.c_str());





	map<int, string * > * egt_alignment = new map<int, string * >;




	// ---------------------------------------------------------------------------
	// move EGTs from pfam_alignment to egt_alignment
	for (ali_it_pfam = pfam_alignment->begin(); ali_it_pfam != pfam_alignment->end(); ) {
			putQuestionmarks(ali_it_pfam->second);
			//cout << ali_it_pfam->first << "\t" << *(ali_it_pfam->second) << endl;

			if (ali_it_pfam->first.compare(0,4, "EGT_")  ==  0) {
					//cerr << "found EGT!!!" << endl;
					//cerr << "EGT_: " << ali_it_pfam->first << endl;

					string number_str = ali_it_pfam->first.substr(4);
					int number_int = str2int(number_str);


					egt_alignment->insert(pair<int, string * >(number_int, ali_it_pfam->second));

					pfam_alignment->erase(ali_it_pfam++);

			} else {
				++ali_it_pfam;
			}
	}





	// filter
	bool filter_metagenom = false;


	if (filter_metagenom) {
			cerr << "warning: filter for evaluation purposes is activated" << endl;

			string filter_rank_string = "species";

			//vector<int> filtered_species_vec (


			tax_id filtered_species[]  =	 {
								381754 , //"Pseudomonas aeruginosa PA7",
								158878 ,// "Staphylococcus aureus subsp. aureus Mu50",
								557723 , // "Haemophilus parasuis SH0165",
								413997 , // "Escherichia coli B str. REL606",
								266834 , // "Sinorhizobium meliloti 1021",
								592205 , // "Helicobacter pylori B38",
								316057 , // "Rhodopseudomonas palustris BisB5",
								264462 , // "Bdellovibrio bacteriovorus HD100",
								387093 , // "Sulfurovum sp. NBC37-1",
								167879 , // "Colwellia psychrerythraea 34H",
								579112 , // "Vibrio cholerae M66-2",
								223926 , // "Vibrio parahaemolyticus, RIMD 2210633",
								357809 , // "Clostridium phytofermentans ISDg",
								334380 , // "Orientia tsutsugamushi str. Ikeda",
								176280 , // "Staphylococcus epidermidis ATCC 12228",
								395019 , // "Burkholderia multivorans ATCC 17616",
								43989 , // "Cyanothece sp. ATCC 51142",
								316275 , // "Aliivibrio salmonicida LFI1238",
								272561 , // "Chlamydia trachomatis D/UW-3/CX",
								320388 , // "Burkholderia mallei SAVP1",
								344609 , // "Shigella boydii CDC 3083-94",
								561007 , // "Mycobacterium abscessus ATCC 19977",
								359391 , // "Brucella melitensis biovar Abortus 2308",
								512566 , // "Streptococcus pneumoniae G54",
								110662  // "Synechococcus sp. CC9605"
								};

			vector<tax_id> * filtered_species_vec = new vector<tax_id>(filtered_species,filtered_species+25);

			rank_t filter_rank = (*rank_to_id)[filter_rank_string];

			if (filtered_taxa == 0) {
				filtered_taxa = init_filtered_taxa(filtered_species_vec, filter_rank , rank_to_id, parent_taxid, ranks, NCBI_MAX);
			}


			for (ali_it_pfam = pfam_alignment->begin(); ali_it_pfam != pfam_alignment->end(); ) {



						string description = ali_it_pfam->first;
						size_t slash = description.find_first_of('/');
						if (slash == string::npos) {
							cerr << "(A) parsing error, no slash found:" << endl;
							cerr << description << endl;
							exit(1);
						}

						description=description.substr(0, slash);
						//cout << description << endl;

						map<string, tax_id >::iterator taxid_it;

						// determine ncbi tax_id of family member:
						taxid_it = pfamid2taxid_map->find(description);

						bool filter_taxon = false;

						if (taxid_it != pfamid2taxid_map->end()) {
							filter_taxon = isFilteredTaxon(filter_rank, taxid_it->second, &this->filtered_taxa, this->rank_to_id, this->parent_taxid, this->ranks, this->NCBI_MAX);

						} else {
							// tax_id unknown, filter it away...
							filter_taxon = true;
						}

						if (filter_taxon) {
							pfam_alignment->erase(ali_it_pfam++);
						} else {
							++ali_it_pfam;
						}

			}

	}

	if (egt_alignment->size() != count_egts) {
			cerr << "egt_alignment->size(): " << egt_alignment->size() << endl;
			cerr << "count_egts: " << count_egts << endl;
			cerr << "error, number of egts before hmmalign and after hmmalign is not the same." << endl;
			cerr << "maybe some EGTs got lost or non-unique names were used for the EGTs!?" << endl;


			//string cp_comm = "cp ";
			//cp_comm.append(new_alignment_file_name);
			//cp_comm.append(" /vol/cluster-data/wgerlach/test/classify_output/");
			//cp_comm.append(" .");
			//cerr << cp_comm << endl;
			//system(cp_comm.c_str());

			exit(EXIT_FAILURE);
	}


	// ---------------------------------------------------------------------------
	// process each egt:

	tax_id * best_hit_taxonomy=0;
	tax_id * other_hit_taxonomy=0;

	map<int, string * >::iterator ali_it_egt;
	map<string, string * >::iterator best_alignment;

	for (ali_it_egt = egt_alignment->begin(); ali_it_egt != egt_alignment->end(); ali_it_egt++) {
		//cout << ali_it_egt->first << "\t" << *(ali_it_egt->second) << endl;
		//cout << "EGT -----------------------------------"<< endl;
		//cout << egt_alignment->size() << endl;

		int number_int = ali_it_egt->first;
		string descr = vec->at(number_int)->description;
		vector<string> * egt_descr_parts = getEGTDescription(descr);


		string pfam_descr = (*egt_descr_parts)[0];
		string read_name = (*egt_descr_parts)[1];
		string evalue_string = (*egt_descr_parts)[2];
		//string go_ids = (*egt_descr_parts)[3];

		delete egt_descr_parts;


		int best_score = INT_MIN;
		tax_id best_taxid = 0;



		// search for best matching pfam-sequence
		for (ali_it_pfam = pfam_alignment->begin(); ali_it_pfam != pfam_alignment->end(); ali_it_pfam++) {

			int score = computeBlosumScore(ali_it_egt->second, ali_it_pfam->second, 0, INT_MAX);
			//cout << "score: " << score << endl;
			if (score > best_score) {

					// check for valid taxon:
					//cout << ali_it_pfam->first << endl;
					string description = ali_it_pfam->first;
					size_t slash = description.find_first_of('/');
					if (slash == string::npos) {
						cerr << "(B) parsing error, no slash found:" << endl;
						cerr << description << endl;
						exit(1);
					}

					description=description.substr(0, slash);
					//cout << description << endl;

					map<string, tax_id >::iterator taxid_it;

					// determine ncbi tax_id of family member:
					taxid_it = pfamid2taxid_map->find(description);

					if (taxid_it != pfamid2taxid_map->end()) {

						best_score = score;
						best_alignment = ali_it_pfam;
						best_taxid = taxid_it->second;

					}


			}
			//cout << "score: " << score << endl;
		}
		//cout << "best_score: " << best_score << endl;

		if (best_score < this->pairwise_blosum_minscore) {
		//if (best_score < 0) {

			ClassificationResult * classification_result = new ClassificationResult(read_name, pfam_descr, 0, 0, evalue_string);
			printClassificationResult(classification_result);
			delete classification_result;
			continue;

		}


		// --------------------------------------------
		// get borders of EGT ??xxx??
		size_t start_area = ali_it_egt->second->find_first_not_of('?');

		if (start_area == string::npos) {
			cerr << "no start area found" << endl;
			exit(1);
		}

		size_t stop_area = ali_it_egt->second->find_last_not_of('?');

		if (stop_area == string::npos) {
			cerr << "no stop area found" << endl;
			exit(1);
		}
		// --------------------------------------------




		getTaxa(best_taxid, &best_hit_taxonomy, this->parent_taxid, NCBI_MAX);


		ReciprocResults * recRes = new ReciprocResults();
		recRes->match_type = MATCH_TYPE_PFAM_HMM;
		recRes->query_bitscore = best_score;
		recRes->self_bitscore = computeBlosumScore(best_alignment->second, best_alignment->second, start_area, stop_area);
		recRes->best_hit_taxid = best_taxid;
		recRes->other_hits = new vector<pair<double, rank_t> * >();

		//cout << start_area << " " << stop_area << endl;
		//cout << best_alignment->second->substr(start_area, stop_area-start_area+1) << endl;
		if (recRes->self_bitscore < this->pairwise_blosum_minscore) { // should not happen
			ClassificationResult * classification_result = new ClassificationResult(read_name, pfam_descr, 0, 0, evalue_string);
			printClassificationResult(classification_result);
			delete classification_result;
			delete recRes;
			continue;
		}

		// reciprocal search:
		for (ali_it_pfam = pfam_alignment->begin(); ali_it_pfam != pfam_alignment->end(); ali_it_pfam++) {

			int score = computeBlosumScore(best_alignment->second, ali_it_pfam->second, start_area, stop_area);

			if (score >= this->pairwise_blosum_minscore) {
				string description = ali_it_pfam->first;

				size_t slash = description.find_first_of('/');
				if (slash == string::npos) {
					cerr << "(C) parsing error, no slash found:" << endl;
					cerr << description << endl;
					exit(1);
				}

				description=description.substr(0, slash);
				//cout << description << endl;

				map<string, tax_id >::iterator taxid_it;

				taxid_it = pfamid2taxid_map->find(description);

				if (taxid_it != pfamid2taxid_map->end()) {
					//cout << "best vs all: "<< score << endl;
					//cout << taxid_it->second << endl;

					// add to recRes


					getTaxa(taxid_it->second, &other_hit_taxonomy, this->parent_taxid, NCBI_MAX);

					tax_id lct = getLowestCommonTaxID(best_hit_taxonomy, other_hit_taxonomy, NCBI_MAX);

					//cout << "last common tax_id: " << lct << endl;
					rank_t lct_rank = getLowestRank(lct, parent_taxid, this->ranks, NCBI_MAX);


					pair<double, rank_t> * ppp = new pair<double, rank_t>(score, lct_rank);

					recRes->other_hits->push_back(ppp);


				}
			}

		} // end pfam-for-loop

		//cerr << "start classifi" << endl;

		// prepare fall-back method:
		Blast_HSP * fake_blast_hsp = new Blast_HSP();
		setAlignmentParameters(fake_blast_hsp, ali_it_egt->second, best_alignment->second, start_area, stop_area);


		ClassificationResult * classification_result = coreClassificationAlgorithm(read_name, evalue_string, fake_blast_hsp, recRes, true);

		if (classification_result != 0) {

			classification_result->property = pfam_descr;
			printClassificationResult(classification_result);

		} else {
			classification_result = new ClassificationResult(read_name, pfam_descr, 0, 0, evalue_string);
			printClassificationResult(classification_result);

		}

		if (recRes !=0) {
			delete recRes;
		}

		delete classification_result;
	} // end EGT-for-loop

	DeleteMapWithPointers(pfam_alignment);
	DeleteMapWithPointers(egt_alignment);

	if (best_hit_taxonomy != 0 ){
		delete [] best_hit_taxonomy;
	}
	if (other_hit_taxonomy != 0 ){
		delete [] other_hit_taxonomy;
	}

	return;
}


void CARMA_HMMER::processEGTs(map<string, vector<Sequence * > * > * sequencesByFamily) {

		//cout << "blubb" << endl;



		map<string, vector<Sequence * > * >::iterator family_it;

		for (family_it = sequencesByFamily->begin(); family_it != sequencesByFamily->end(); family_it++) {
				//cout << family_it->first << " " << family_it->second->size() << endl;
				vector<Sequence * > * vec = family_it->second;
				processEGTfamily(family_it->first, vec);
		}


}
