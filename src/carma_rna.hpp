
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

#ifndef GUARD_CARMA_RNA
#define GUARD_CARMA_RNA



#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>

#include "carma_base.hpp"
#include "biotools.hpp"
#include "biotools_rdp.hpp"


using namespace std;

class CARMA_BASE;


class CARMA_RNA : public CARMA_BASE {

	private:

	rank_t superkingdom_rank;
	rank_t phylum_rank;
	rank_t class_rank;
	rank_t order_rank;
	rank_t family_rank;
	rank_t genus_rank;
	rank_t species_rank;



	public:

	//vector< string * > * rna_alignment_vec;
	map<string, Sequence_with_taxid * > * rdp_alignment_bacteria;
	map<string, Sequence_with_taxid * > * rdp_alignment_archea;

	string * rdp_arch_aligned;
	string * rdp_bact_aligned;

	double blast_rdp_evalue;

	CARMA_RNA(string * config_file) : CARMA_BASE(config_file, 0){
		setVariables();


	}

	void setVariables();


	void extractSpecificRDP(string file,  int expected_length, bool isAlignment, bool removeGapColumns, vector<int> * taxid_vec, string filter_rank_string);

	//void parse_RDP_description_line(string& line, string ** identifier, tax_id* taxid, int * offset, map<string, tax_id> * NameToNCBI );
	void parse_16S_blast(string file);
	void classify_16S_Fragment(vector<vector<string> * > * blast_hits);
	int computePairwiseDNAScore(Sequence_with_taxid * seq_A, Sequence_with_taxid * seq_B, int start, int end);



	map<string, Sequence_with_taxid * > * read_RNA_Alignment(string file);

	bool is_specific_taxon(tax_id taxid, map<tax_id, int > * phylum_counts);

	void rdp_train_data(string * input_file);
	int ranktid_to_rdp_rank_id(rank_t taxrank);
	string ranktid_to_rankname(rank_t taxrank);
	int get_actual_RDP_rank(tax_id taxid);

};


#endif
