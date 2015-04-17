
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

#include "carma_rna.hpp"

void CARMA_RNA::setVariables(){


	this->rdp_arch_aligned = getConfigStringValue(parameters, "rdp_arch_aligned");
	//checkFileWarning(this->rdp_arch_aligned);

	this->rdp_bact_aligned = getConfigStringValue(parameters, "rdp_bact_aligned");
	//checkFileWarning(this->rdp_bact_aligned);



	this->blast_rdp_evalue = getConfigDoubleValue(parameters, "blast_rdp_evalue");



}



void CARMA_RNA::parse_16S_blast(string file){

	this->rdp_alignment_archea = read_RNA_Alignment(*(this->rdp_arch_aligned));
	this->rdp_alignment_bacteria = read_RNA_Alignment(*(this->rdp_bact_aligned));


	BLAST_M9_Parser * parser = new BLAST_M9_Parser(file, zcat_bin);

	vector<vector<string> * > * blast_hits = 0;
	string query_id;
	while (parser->getNextBLASTHits(&blast_hits, query_id)) {

		//cerr << "q: "<< blast_hits->at(0)->at(0) << endl;
		//cerr << "number of hits: "<< blast_hits->size() << endl;


		if ( blast_hits->size() > 0 ) {
			classify_16S_Fragment(blast_hits);

		} else {
			string prop = "carma_RDP";
			string evalue_none = "-";

			ClassificationResult * classification_result = new ClassificationResult(query_id, prop, 0, 0, evalue_none );
			printClassificationResult(classification_result);
			delete classification_result;
		}


		DeleteContainerWithPointers(blast_hits);
		blast_hits=0;



	}



}



//void CARMA_RNA::classify_16S_Fragment(vector<string> * line_data, vector<string> * other_blast_hits) {
void CARMA_RNA::classify_16S_Fragment(vector<vector<string> * > * blast_hits) {


	map<string, Sequence_with_taxid * > * rdp_alignment = 0;
	map<string, Sequence_with_taxid * >::iterator ali_it;


	vector<string> * line_data = blast_hits->at(0);


	//cerr << "xxx" << line_data->size() <<endl;
	// Fields: (0) Query id, (1)Subject id, (2)% identity, (3)alignment length, (4)mismatches, (5)gap openings, (6)q. start, (7)q. end, (8)s. start, (9)s. end, (10)e-value, (11)bit score
	string query_id = (*line_data)[0];
	string subject = (*line_data)[1];


	vector<string> * data_field = parse_column_data_line(subject, '_');

	if (data_field->size() != 2) {
		cerr << "error: data_field->size() != 2" << endl;
		exit(1);
	}

	subject = data_field->at(0);

	delete data_field;

//cerr << "subj: "<< subject << endl;
//exit(1);

// ********** get alignment borders *********
	int subject_start = str2int((*line_data)[8]);
	int subject_end = str2int((*line_data)[9]);

	string evalue_string = (*line_data)[10];
	//double best_evalue_double = string2double(evalue_string);

	string best_bitscore_str = (*line_data)[11];
	int best_bitscore = str2int(best_bitscore_str);


//cerr << "yyy" << endl;
	if (subject_start > subject_end) {
		cerr << "warning: subject_start > subject_end, subj: "<< subject << endl;
		swap(subject_start, subject_end);
		//cerr << "subj: "<< subject << endl;
		//exit(1);
	}

	if (str2int((*line_data)[6]) > str2int((*line_data)[7])) {
		cerr << "query_start > query_end" << endl;
		exit(1);
	}

//	cerr << "query_id: "<< query_id << endl;

//	cerr << "subj_start_end: "<< subject_start<< " - " << subject_end << endl;


	// ********** find alignment sequence of best hit (subject)*********

	ali_it = this->rdp_alignment_bacteria->find(subject);
	char superkingdom = '-';
	if (ali_it == this->rdp_alignment_bacteria->end()) {

		ali_it = this->rdp_alignment_archea->find(subject);

		if (ali_it == this->rdp_alignment_archea->end()) {
			cerr << "error: Subject \"" << subject << "\" not found in full alignment (bacteria or archea)" << endl;
			exit(1);
		}
		superkingdom = 'a';
		rdp_alignment = this->rdp_alignment_archea;
	} else {
		superkingdom = 'b';
		rdp_alignment = this->rdp_alignment_bacteria;
	}

//	cerr << "Superkingdom: " << superkingdom << endl;

	Sequence_with_taxid * best_hit_alignment_sequence = ali_it->second;

	if (best_hit_alignment_sequence->taxid == 0) {
			cerr << "error: best_hit_alignment_sequence->taxid == 0" << endl;
			exit(1);
	}


// ********** check uniqueness of best hit regarding alignment class (archea and bacteria) *********

	vector<vector<string> * >::iterator other_blast_hits_it;

	bool is_unique = true;
	other_blast_hits_it = blast_hits->begin();
	++other_blast_hits_it; // skip first

	for (; other_blast_hits_it != blast_hits->end(); ++other_blast_hits_it) {
		string other_id = (*other_blast_hits_it)->at(1);

		data_field = parse_column_data_line(other_id, '_');

		if (data_field->size() != 2) {
			cerr << "error: data_field->size() != 2" << endl;
			exit(1);
		}

		other_id = data_field->at(0);

		delete data_field;


		string other_bitscore_str = (*other_blast_hits_it)->at(11);
		int other_bitscore = str2int(other_bitscore_str);
		char other_superkingdom = '-';
		if (other_bitscore == best_bitscore) {

			ali_it = this->rdp_alignment_bacteria->find(other_id);
			if (ali_it == this->rdp_alignment_bacteria->end()) {
				ali_it = this->rdp_alignment_archea->find(other_id);

				if (ali_it == this->rdp_alignment_archea->end()) {
					cerr << "error: Other hit \"" << other_id << "\" not found in full alignment (bacteria or archea)" << endl;
					exit(1);
				} else {
					other_superkingdom = 'a';
				}
			} else {
				other_superkingdom = 'b';
			}

			if (superkingdom != other_superkingdom) {
				is_unique = false;
				//cerr << "NOT UNIQUE !!!" << endl;
				break;
			}

		} else if (other_bitscore < best_bitscore) {
			//is_unique = false;
			break;
		}

	}

	ClassificationResult * classification_result = 0;
	ReciprocResults * recRes = 0;

	if (is_unique) {



		//cerr << "taxid: " << best_hit_alignment_sequence->taxid << endl;
		//cerr << "best hit: " << best_hit_alignment_sequence->substr(best_hit_alignment_sequence->offset, 20) << endl;

		// ******** score of query **********


		int query_score =  getQueryScore(line_data, score_match, score_mismatch, score_gapopen, score_gapextension);


		// ******** get alignment borders *************

		size_t alignment_start=-1;
		for (int i = 0; i < subject_start; ++i) {
			alignment_start = best_hit_alignment_sequence->find_first_not_of("-?", alignment_start+1);
			if (alignment_start == string::npos) {
				cerr << "error: alignment ended too early... (start)" << endl;
				exit(1);
			}
		}
		alignment_start+=best_hit_alignment_sequence->offset; // because the find-function is not aware of the offset
		//cerr << "alignment_start: " << alignment_start << endl;


		int length_of_aligned_subject = subject_end-subject_start+1;
		size_t alignment_end = alignment_start;
		int length=0;
		while (length < length_of_aligned_subject ) {

			if (alignment_end >= (best_hit_alignment_sequence->length()  ) ) {
				cerr << "error: alignment ended too early... (end)" << endl;
				cerr << "alignment_end: " << alignment_end << endl;
				cerr << "best_hit_alignment_sequence->length(): " << best_hit_alignment_sequence->length() << endl;
				cerr << "length: " << length << endl;
				cerr << "length_of_aligned_subject: " << length_of_aligned_subject << endl;

				exit(1);
			}
			char c = (*(best_hit_alignment_sequence))[alignment_end];
			if (c != '-' && c != '?') {
					length++;
					//cerr << c;
			}
			alignment_end++;
		}
		//cerr << endl;
		//cerr << "alignment_end: " << alignment_end << endl;


		// ******************************** reciprocal search **********************************************
		tax_id * best_hit_taxonomy=0;
		tax_id * other_hit_taxonomy=0;

		getTaxa(best_hit_alignment_sequence->taxid, &best_hit_taxonomy, this->parent_taxid, NCBI_MAX);

		recRes = new ReciprocResults();
		recRes->match_type = MATCH_TYPE_RDP;
		recRes->query_bitscore = query_score;
		recRes->self_bitscore = computePairwiseDNAScore(best_hit_alignment_sequence, best_hit_alignment_sequence, alignment_start, alignment_end);
		//cerr << "recRes->self_bitscore: " << recRes->self_bitscore << endl;

		recRes->best_hit_taxid = best_hit_alignment_sequence->taxid;
		//cerr << "recRes->best_hit_taxid: " << recRes->best_hit_taxid << endl;
		recRes->other_hits = new vector<pair<double, rank_t> * >();


		//map<string, Sequence_with_taxid * >::iterator ali_it;
		int good_score=0;
		//cerr << "reciprocal search:" << endl;

		//vector<vector<string> * >::iterator other_blast_hits_it;


		other_blast_hits_it = blast_hits->begin();
		++other_blast_hits_it; // skip first

		for (; other_blast_hits_it != blast_hits->end(); ++other_blast_hits_it) {
			string other_id = (*other_blast_hits_it)->at(1);
			data_field = parse_column_data_line(other_id, '_');

			if (data_field->size() != 2) {
				cerr << "error: data_field->size() != 2" << endl;
				exit(1);
			}

			other_id = data_field->at(0);

			delete data_field;


			//cerr << "other_id: " << other_id << endl;
			ali_it = rdp_alignment->find(other_id);
			if (ali_it == rdp_alignment->end()) {
				cerr << "warning: " << other_id << " not found in rdp_alignment (" << superkingdom <<")..." << endl;
				continue; // might be a hit from other superkingdom
				//cerr << "error: " << other_id << " not found in rdp_alignment..." << endl;
				//exit(1);
			}

		//}

		//for (ali_it = rdp_alignment->begin(); ali_it != rdp_alignment->end(); ++ali_it) {
			Sequence_with_taxid * other_alignment_sequence = ali_it->second;

			int score = computePairwiseDNAScore(best_hit_alignment_sequence, other_alignment_sequence, alignment_start, alignment_end);

			if (score > 0) {
				good_score++;
				//cerr << "score: " << score << endl;
				//cerr << "other_alignment_sequence->taxid: " << other_alignment_sequence->taxid << endl;
				getTaxa(other_alignment_sequence->taxid, &other_hit_taxonomy, this->parent_taxid, NCBI_MAX);

				tax_id lct = getLowestCommonTaxID(best_hit_taxonomy, other_hit_taxonomy, NCBI_MAX);

				//cout << "last common tax_id: " << lct << endl;
				rank_t lct_rank = getLowestRank(lct, parent_taxid, this->ranks, NCBI_MAX);


				pair<double, rank_t> * ppp = new pair<double, rank_t>(score, lct_rank);

				recRes->other_hits->push_back(ppp);

			}

	//		if (best_hit_alignment_sequence == other_alignment_sequence) {
	//			cerr << "self: " << score << endl;
	//
	//		}

			bool showali = false;
			if (showali) {
				cerr << best_hit_alignment_sequence->substr(alignment_start, 50) << endl;
				cerr << other_alignment_sequence->substr(alignment_start, 50) << endl;
			}
				// self score !
			//}


		}
		//cerr << "number of sequences with high score: " << good_score << endl;


		classification_result = coreClassificationAlgorithm(query_id, evalue_string, 0, recRes, false);


	} // end if is_unique



	if (classification_result != 0) {

		//classification_result->property = pfam_descr;
		printClassificationResult(classification_result);

	} else {
		//cerr << "error: no classification result" << endl;
		//exit(1);
		string prop = "carma_RDP";
		classification_result = new ClassificationResult(query_id, prop, 0, 0, evalue_string);
		printClassificationResult(classification_result);

	}

	if (recRes !=0) {
		delete recRes;
	}

	if (classification_result != 0) {
		delete classification_result;
	}
}


int CARMA_RNA::computePairwiseDNAScore(Sequence_with_taxid * seq_A, Sequence_with_taxid * seq_B, int start, int end) {
	int len_a = seq_A->length();
	int len_b = seq_B->length();

	int pos = start;
//cout << "A: "<< *seq_A << endl;
//cout << "B: "<< *seq_B << endl << endl;
	int sum = 0;


	bool gap=false;
	bool gap_on_A = false;


	if ((end < seq_A->offset) || (end < seq_B->offset)) {
		return 0;
	}


	while (pos < len_a && pos < len_b && pos <= end) {
			char a = tolower((*seq_A)[pos]);
			char b = tolower((*seq_B)[pos]);

			// alphabet: actgu-n
			if (a=='u') {
				a = 't';
			}
			if (b=='u') {
				b = 't';
			}

			if ((a=='?') || (b == '?')) {
				// change nothing
			} else if ((a=='-') && (b == '-')) {
				// change nothing
			} else if ((a=='n') && (b == 'n')) {
				gap = false;
				sum+=score_mismatch;
			} else if ((a=='-') || (b == '-')) { // one gap
				if (gap) {
					// check if this is still the same gap
					if ( ((a=='-') && (gap_on_A==true)) || ((b=='-') && (gap_on_A==false))) {
						sum += score_gapextension;
					} else {
						// start new gap on other sequence:
						sum += score_gapopen;
						gap_on_A = !gap_on_A; // toggle gap
					}

				} else {
					// start gap
					sum += score_gapopen;
					gap=true;
					if (a=='-') {
						gap_on_A = true;
					} else {
						gap_on_A = false;
					}
				}
			} else if (a == b) {
				gap = false;
				sum += score_match;
			} else {
				gap = false;
				sum += score_mismatch;
			}


			pos++;

	}


	return sum;
}

map<string, Sequence_with_taxid * > * CARMA_RNA::read_RNA_Alignment(string file){

//cerr << sizeof(unsigned long int)<< endl;
//exit(0);


	map<string, Sequence_with_taxid * > * rdp_alignment = new map<string, Sequence_with_taxid * >;
	//rdp_alignment_bacteria_vec = new vector< string * >;

	//rdp_alignment_bacteria_vec->reserve ( 350000 );


	FASTA_Parser * rdp_parser = new FASTA_Parser(file, false, zcat_bin);
	rdp_parser->expected_sequence_length = 30000;

	string line;
	string * identifier = 0;
	int offset;
	int total_count = 0;
	while ( rdp_parser->getNextDescriptionLine(line) ) {
		total_count++;

		if (total_count >= 50000 && false) {
			cerr << "Warning: file has been read only partial..." << endl;
			break;
		}


		if (total_count % 50000 == 0) {
			cerr << "total: " <<total_count << endl;
		}

		string * rna_seq = rdp_parser->getSequence();

		tax_id taxid = 0;

		//parse_RDP_description_line(line, &identifier, &taxid, &offset, NameToNCBI );
		vector<string> * data_field = parse_column_data_line(line, ' ');
		if (data_field->size() != 2) {
			cerr << "error: data_field->size() != 2" << endl;
			exit(1);
		}
		string part1_2 = data_field->at(0);
		//cerr << part1_2 << endl;
		string part_3 = data_field->at(1);
		//cerr << part_3 << endl;
		delete data_field;

		string offset_string = part_3.substr(6);
		offset = str2int( offset_string );

		//cerr << offset << endl;

		vector<string> * data_field_2 = parse_column_data_line(part1_2, '_');
		if (data_field_2->size() != 2) {
			cerr << "error: data_field_2->size() != 2" << endl;
			exit(1);
		}

		string part_1 = data_field_2->at(0);
		string part_2 = data_field_2->at(1);
		delete data_field_2;


		part_1 = part_1.substr(1);
		identifier = new string(part_1);

		//cerr << *identifier << endl;

		taxid = str2int(part_2);
		//cerr << taxid << endl;
//exit(1);
		if (taxid == 0) {
			cerr << "error: taxid == 0" << endl;
			exit(1);
		}

//cerr <<"identifier: " << *identifier<<  endl;
//cerr <<"offset: " << offset<<  endl;

		Sequence_with_taxid * seq_tax = new Sequence_with_taxid();
		//cerr << "len: " << rna_seq->length() << endl;
		seq_tax->assign(*rna_seq);
		seq_tax->reserve(rna_seq->length());
		seq_tax->offset=offset;
		seq_tax->taxid = taxid;

		//cerr << "cap: " << seq_tax->capacity() << endl;
//		if (identifier->compare("S001043913")==0) {
//				cerr << "S001043913 insert: "<< rna_seq->substr(0,1200) << endl;
//		}

		delete rna_seq;

		map<string, Sequence_with_taxid * >::iterator ali_it = rdp_alignment->find(*identifier);
		if (ali_it != rdp_alignment->end()) {
			cerr << "Sequence already in alignment: " << *identifier << endl;
			exit(1);
		}

		rdp_alignment->insert( pair<string, Sequence_with_taxid * > (*identifier, seq_tax) );

		delete identifier;
		identifier = 0;

		//rdp_alignment_bacteria_vec->push_back(rna_seq);

		//cerr << rna_seq->length() << endl;
		//exit(0);
		//cerr << "\""<<  identifier << "\"" << endl;
		//cout << *identifier << endl;
		//cout << "taxid: " << taxid << endl;
		//cout << rna_seq->substr(0,50) << endl;
		//cout << *rna_seq << endl;
	}
	cerr << "have put alignment into memory: " << total_count << endl;
	delete rdp_parser;

	return rdp_alignment;
}


bool CARMA_RNA::is_specific_taxon(tax_id taxid, map<tax_id, int > * phylum_counts){

	// phylum, genus and species must exist
	// genus should not be "enviornmental samples"


	if (taxid == 0) {
		return false;
	}


	// check for species
	rank_t taxon_rank = getLowestRank(taxid, parent_taxid, ranks, NCBI_MAX);
	if (taxon_rank < species_rank) {
		//unspecific_count++;

		return false;
	}


	// check that pylum exists
	tax_id phylum = getTaxonAtRank(taxid, phylum_rank , parent_taxid, ranks, NCBI_MAX);
	//cout << "phylum: " << (int) phylum << endl;

	rank_t verifiy_rank = getLowestRank(phylum, parent_taxid, ranks, NCBI_MAX);
	if (verifiy_rank == phylum_rank) {

		map<tax_id, int >::iterator phy_it;

		if (phylum_counts != 0) {
			phy_it = phylum_counts->find(phylum);

			if (phy_it != phylum_counts->end()) {
					(*phy_it).second++;
			} else {
				phylum_counts->insert( pair<tax_id, int >(phylum, 1) ) ;
			}
		}
	} else {
		//cout << line;
		//cout << "no phylum: " <<  taxid << endl;
		//no_phylum_count++;
		return false;
	}

	// check that genus exists
	tax_id genus = getTaxonAtRank(taxid, genus_rank , parent_taxid, ranks, NCBI_MAX);


	if (taxid_to_name == 0) {
		taxid_to_name = init_taxid_to_name(names_dmp->c_str());
	}

	map<tax_id, string >::iterator map_it;
	map_it=this->taxid_to_name->find(genus);
	if (map_it==this->taxid_to_name->end()) {
		cerr << "error: map_it==this->taxid_to_name->end()" << endl;
		exit(1);
	}
	string genus_name = map_it->second;
	//cerr << "name: " << genus_name << endl;

	verifiy_rank = getLowestRank(genus, parent_taxid, ranks, NCBI_MAX);
	if (verifiy_rank != genus_rank) {
		//no_genus_count++;
		return false;
	}

	//cerr << "genus: " << (int) genus << endl;
	//cerr << "name: " << genus_name << endl;

	if (genus_name.compare("environmental samples")==0) {
		return false;
		//cerr << "genus: " << (int) genus << endl;
		//cerr << "name: " << genus_name << endl;

	}

	return true;

}



void CARMA_RNA::extractSpecificRDP(string file, int expected_length, bool isAlignment, bool removeGapColumns, vector<int> * taxid_vec, string filter_rank_string){

	rank_t filter_rank;
	if (taxid_vec != 0) {

			cerr << "WARNING: filter activated" << endl;
			cerr << "filter_rank_string: " << filter_rank_string << endl;



			//string filter_rank_string = "species";
			if (filter_rank_string.length() == 0) {
				cerr << "error: filter_rank_string.length() == 0" <<endl;
				exit(1);
			}
			int taxid_num = taxid_vec->size();
			cerr << "taxid_num: " << taxid_num << endl;
			//tax_id filtered_species[taxid_num];
			vector<tax_id> * filtered_species_vec = new vector<tax_id>();

			vector<int>::iterator vec_it;

			for (vec_it=taxid_vec->begin(); vec_it != taxid_vec->end(); vec_it++) {
				filtered_species_vec->push_back( (tax_id) *vec_it );
			}
			cerr << endl;
			//filtered_species.

			map<string, rank_t>::iterator rank_it;
			rank_it = rank_to_id->find(filter_rank_string);

			if (rank_it == rank_to_id->end()) {
				cerr << "error: rank_it == rank_to_id->end()" << endl;
				exit(1);
			}

			filter_rank = rank_it->second;

			if (filtered_taxa == 0) {
				filtered_taxa = init_filtered_taxa(filtered_species_vec, filter_rank , rank_to_id, parent_taxid, ranks, NCBI_MAX);
			}
	}



	int unspecific_count = 0;
	int no_phylum_count = 0;
	int no_genus_count = 0;
	int specifc_count = 0;
	int total_count = 0;

	int alignment_length = 0;

	string * identifier = 0;
	tax_id taxid;
	string descr;

	superkingdom_rank = (*rank_to_id)["superkingdom"];
	phylum_rank = (*rank_to_id)["phylum"];
	class_rank = (*rank_to_id)["class"];
	order_rank = (*rank_to_id)["order"];
	family_rank = (*rank_to_id)["family"];
	genus_rank = (*rank_to_id)["genus"];
	species_rank = (*rank_to_id)["species"];

	//map<string, string * > * rdp_alignment_bacteria = new map<string, string * >;
	//vector< string * > * rdp_alignment_bacteria_vec = new vector< string * >;
	//rdp_alignment_bacteria_vec->reserve ( 300000 );

	map<tax_id, int > * phylum_counts = new map<tax_id, int >;

	set<string> set_of_identifiers;


	FASTA_Parser * rdp_parser = 0;


	vector< bool > * gap_vec = 0;
	if (removeGapColumns) {
		alignment_length = 0;
		// first iteration to find gaps
		rdp_parser = new FASTA_Parser(file, false, zcat_bin);
		if (expected_length > 0) {
			// this is used to tell the strings how much memory will be used (dynamic allocation would slow down)
			rdp_parser->expected_sequence_length = expected_length;
		}


		gap_vec = new vector< bool >;
		gap_vec->reserve(30000);
		for(int i = 0; i < 30000; i++) {
			(*gap_vec)[i] = true; // everything is a gap
		}
		map<string, Sequence_with_taxid * >::iterator ali_it;



		while ( rdp_parser->getNextDescriptionLine(descr) ) {
			total_count++;

			if (total_count >= 100000 && false) {
				cerr << "Warning: file has been read only partial..." << endl;
				break;
			}

			if (total_count % 50000 == 0) {
					cerr << "total (gap-search): " << total_count << endl;
			}
			parse_RDP_description_line(descr, NULL, &taxid, NULL, NameToNCBI );

			bool is_specific = this->is_specific_taxon(taxid, NULL);
			if (is_specific == false) {
				continue;
			}

			string * rna_seq = rdp_parser->getSequence();

			if (rna_seq->length() > alignment_length) {
				alignment_length = rna_seq->length();
			}

			for(int i = 0; i < rna_seq->length(); i++) {
				char c = (*rna_seq)[i];
				if (c != '-' && c != '?' ) {
						(*gap_vec)[i] = false; // this is not a gap
				}
			}

			delete rna_seq;
		}

		int gap_count =0;
		for(int i = 0; i < alignment_length; i++) {
			if ((*gap_vec)[i] == true) {
				gap_count++;
			}
		}
		cerr << "alignment_length: " << alignment_length << endl;
		cerr << "gap_count: " << gap_count << endl;



		delete rdp_parser;
	}




	rdp_parser = new FASTA_Parser(file, false, zcat_bin);
	if (expected_length > 0) {
		// this is used to tell the strings how much memory will be used (dynamic allocation would slow down)
		rdp_parser->expected_sequence_length = expected_length;
	}
	total_count=0;

	bool filter_taxon;

	while ( rdp_parser->getNextDescriptionLine(descr) ) {
		total_count++;
		if (total_count % 50000 == 0) {
				cerr << "tot: " <<total_count << endl;
		}
		parse_RDP_description_line(descr, &identifier, &taxid, NULL, NameToNCBI );

		bool is_specific = this->is_specific_taxon(taxid, phylum_counts);
		if (is_specific == false) {
			continue;
		}

		if (taxid_vec != 0) {
			filter_taxon = isFilteredTaxon(filter_rank, taxid, &this->filtered_taxa, this->rank_to_id, this->parent_taxid, this->ranks, this->NCBI_MAX);
			if (filter_taxon) {
				continue;
			}
		}

		string * rna_seq = rdp_parser->getSequence();
		specifc_count++;
//cerr << "identifier: " <<  *identifier << endl;
		set<string>::iterator set_it;
		set_it = set_of_identifiers.find(*identifier);
		if (set_it != set_of_identifiers.end()) {
			cerr << "error: " << *identifier << " already in set " << endl;
			exit(1);
		}
		set_of_identifiers.insert(*identifier);



//		if (isAlignment) {
//			// replace leading and trailing gaps with quotationmark
//			putQuestionmarks(rna_seq);
//		}



		if (removeGapColumns) {
			int i_write = 0;
			for(int i_read = 0; i_read < rna_seq->length(); i_read++) {

				if ((*gap_vec)[i_read] == false) {
					if (i_write != i_read) {
						(*rna_seq)[i_write]=(*rna_seq)[i_read];
					}
					i_write++;
				}
			}
			rna_seq->resize(i_write);
		}

		if (isAlignment) {
			Sequence_with_taxid seq;
			seq.setAlignmentSequence(*rna_seq);

			//cerr << seq.offset << endl;

			//cout << descr << "; OFFSET" << seq.offset << endl;
			cout << ">" << *identifier << "_" << (int)taxid << " OFFSET" << seq.offset << endl;
			cout << seq << endl;

		} else {

			//cout << descr << endl;
			cout << ">" << *identifier << "_" << (int)taxid << endl;
			cout << *rna_seq << endl;
		}





		delete rna_seq;

	}  // end while





	cerr << "unspecific_count: "<< unspecific_count << endl;
	cerr << "no_phylum_count: "<< no_phylum_count << endl;
	cerr << "no_genus_count: "<< no_genus_count << endl;


	cerr << "specifc_count: "<< specifc_count << endl;
	cerr << "total_count: "<< total_count << endl;

	map<tax_id, int >::iterator phy_it;
	for (phy_it = phylum_counts->begin() ; phy_it != phylum_counts->end(); ++phy_it ) {
		tax_id blubb = (*phy_it).first;
		string phy_name = getTaxonStringByTaxId(blubb, &taxid_to_name, parent_taxid, names_dmp);
		cerr << phy_name << " " <<  blubb << " - " << (*phy_it).second << endl;
	}

}

// expects as input the unaligned RDP sequences, but specfic
void CARMA_RNA::rdp_train_data(string * input_file) {

	//cerr << *input_file << endl;

	string output_taxonomy_file = "taxonomy.txt";
	string output_rdptraining_file = "rdptraining.fas";

	if (FileExists((char *) output_taxonomy_file.c_str())) {
		cerr << "error: file already exists: " << output_taxonomy_file << endl;
		exit(1);
	}

	if (FileExists((char *) output_rdptraining_file.c_str())) {
		cerr << "error: file already exists: " << output_rdptraining_file << endl;
		exit(1);
	}

	ofstream output_rdptraining_stream (output_rdptraining_file.c_str());
	if (! output_rdptraining_stream.is_open() ) {
		cerr << "error openeing rdptraining file for writing" << endl;
		exit(1);
	}


	if (taxid_to_name == 0) {
		taxid_to_name = init_taxid_to_name(names_dmp->c_str());
	}

	superkingdom_rank = (*rank_to_id)["superkingdom"];
	phylum_rank = (*rank_to_id)["phylum"];
	class_rank = (*rank_to_id)["class"];
	order_rank = (*rank_to_id)["order"];
	family_rank = (*rank_to_id)["family"];
	genus_rank = (*rank_to_id)["genus"];
	species_rank = (*rank_to_id)["species"];

	FASTA_Parser * fasta_parser = new FASTA_Parser(*input_file, false, zcat_bin);

	set<tax_id> * used_taxons_set = new set<tax_id>();

	//tax_id * tax_buffer =0;
	string descr;
	string * identifier = 0;
	string * sequence=0;
	tax_id taxid;
	while (fasta_parser->getNextDescriptionLine(descr)) {
		//cerr << descr << endl;

		//parse_RDP_description_line(descr, &identifier, &taxid, 0, NameToNCBI);




		vector<string> * data_field = parse_column_data_line(descr, '_');
		if (data_field->size() != 2) {
			cerr << "error: data_field->size() != 2" << endl;
			exit(1);
		}

		string part_1 = data_field->at(0);
		string part_2 = data_field->at(1);
		delete data_field;


		part_1 = part_1.substr(1);
		identifier = new string(part_1);

		//cerr << *identifier << endl;

		taxid = str2int(part_2);








		//cerr << *identifier << endl;
		//cerr << taxid << endl;
//exit(0);
		//getTaxa(taxid, &tax_buffer, parent_taxid, NCBI_MAX);


		//bool usetaxon = false;
		string taxon="";
		string name;
		//string previous_name="nothing";
		tax_id temp_tax_id = taxid;
		map<tax_id, string>::iterator it;
		int cur_rank_int;
		//int previous_rank_int = -1;
		while (temp_tax_id > 1) {

			it = taxid_to_name->find(temp_tax_id);

			if (it == taxid_to_name->end()) {
					cerr << "error mapping taxid to name: " << (int)temp_tax_id << endl;
					name = "unknown";
					exit(1);
			} else {
					name = int2str((int)temp_tax_id);
					name.append("_");
					name.append(it->second);
			}

			//usetaxon = false;
			rank_t cur_rank = ranks[temp_tax_id];

			cur_rank_int = ranktid_to_rdp_rank_id(cur_rank);
			if (cur_rank_int != -1) {

				//cerr << "previous_rank_int: " << (int) previous_rank_int << endl;
				//cerr << "cur_rank_int: " << cur_rank_int << endl;


				//cerr << "superkingdom_rank: " << (int) superkingdom_rank << endl;
				//cerr << "phylum_rank: " << (int) phylum_rank << endl;
				//cerr << "class_rank: " << (int) class_rank << endl;


				used_taxons_set->insert(temp_tax_id);

				if (taxon.length()==0) {
					taxon = name ;

				} else {
					taxon = name + ";" + taxon;

				}

			}
				//cout << "found: " << it->second << endl;

			temp_tax_id = parent_taxid[temp_tax_id];
			//previous_rank_int = cur_rank_int;
			//previous_name = name;
		}
		taxon = "1_root;" + taxon;
		//cerr <<"t: " << taxon << endl;

		sequence = fasta_parser->getSequence();

//cerr << ">" << *identifier << " " << taxon << endl;
//exit(1);
		output_rdptraining_stream << ">" << *identifier << " " << taxon << endl;
		output_rdptraining_stream << *sequence << endl;


		delete sequence;
		sequence=0;


	} // end while fasta_parser

	delete fasta_parser;

	output_rdptraining_stream.close();


	ofstream output_taxonomy_stream (output_taxonomy_file.c_str());
	if (! output_taxonomy_stream.is_open() ) {
		cerr << "error openeing taxonomy file for writing" << endl;
		exit(1);
	}

	//root
	output_taxonomy_stream << "1" << "*" << "1_root" << "*" << "0" << "*" << "0" << "*" << "norank" << endl;

	set<tax_id>::iterator set_it;

	for(set_it=used_taxons_set->begin(); set_it!=used_taxons_set->end(); ++set_it) {

		tax_id taxid = *set_it;
//cerr << "taxid: "<< taxid << endl;
		string taxname = (*taxid_to_name)[taxid];



		//int int_rank = ranktid_to_rdp_rank_id(taxrank);

		tax_id parent = parent_taxid[taxid];

		int actual_rank = get_actual_RDP_rank(taxid);

		if (parent== 131567) { // cellular organism
			parent = 0;
		}

		if (taxid == 2 || taxid == 2157 || taxid==2759) { // superkingdom has root=1 as parent
			parent = 1;
		}


		while (used_taxons_set->count(parent)==0 && taxid != 2 && taxid != 2157 && taxid!=2759) { // archea=2157 and bacteria=2 and eukaryota=2759
			if (parent == 1) {
				cerr << "parent==1 " << endl;
				exit(1);
			}
			parent = parent_taxid[parent];
			//cerr << "new: " << (int) parent << endl;
//exit(1);
		}

		rank_t taxrank = ranks[taxid];
		string rankname = ranktid_to_rankname(taxrank);

		output_taxonomy_stream << taxid << "*" << taxid<< "_" << taxname << "*" << parent << "*" << actual_rank << "*" << rankname << endl;
	}



	output_taxonomy_stream.close();


	delete used_taxons_set;

	cerr << "created file: " << output_taxonomy_file << endl;
	cerr << "created file: " << output_rdptraining_file << endl;
}



string CARMA_RNA::ranktid_to_rankname(rank_t taxrank){


	if (taxrank == superkingdom_rank) {
		return "superkingdom";
	} else if (taxrank == phylum_rank) {
		return "phylum";
	} else if (taxrank == class_rank) {
		return "class";
	} else if (taxrank == order_rank) {
		return "order";
	} else if (taxrank == family_rank) {
		return "family";
	} else if (taxrank == genus_rank) {
		return "genus";
	} else {
		cerr << "taxrank error" << endl;
		exit(1);
	}
	return "error";
}


int CARMA_RNA::ranktid_to_rdp_rank_id(rank_t taxrank){
	int int_rank;

	if (taxrank == superkingdom_rank) {
		int_rank = 1;
	} else if (taxrank == phylum_rank) {
		int_rank = 2;
	} else if (taxrank == class_rank) {
		int_rank = 3;
	} else if (taxrank == order_rank) {
		int_rank = 4;
	} else if (taxrank == family_rank) {
		int_rank = 5;
	} else if (taxrank == genus_rank) {
		int_rank = 6;
	} else {
		int_rank = -1;
	}

	return int_rank;
}


int CARMA_RNA::get_actual_RDP_rank(tax_id taxid){

	rank_t taxrank;

	tax_id temp_tax_id = taxid;
	int count_rank = 0;

	while (temp_tax_id > 1) {
		taxrank = ranks[temp_tax_id];

		if ( ranktid_to_rdp_rank_id(taxrank) != -1) {
			count_rank++;
		}

		temp_tax_id = parent_taxid[temp_tax_id];
	}

	if (count_rank > 7) {
		cerr << "error: count_rank > 7" << endl;
		exit(1);
	}

	return count_rank;
}



