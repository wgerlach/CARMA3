
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

#ifndef GUARD_BIOTOOLS
#define GUARD_BIOTOOLS

#include <set>

#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>


#include "PhyloNode.h"




using namespace std;






// ATGC							//   1
const char codon2aminoacid[64] = {	'K', 'N', 'K', 'N',
									'I', 'I', 'M', 'I',
									'R', 'S', 'R', 'S',
									'T', 'T', 'T', 'T',

									'*', 'Y', '*', 'Y',
									'L', 'F', 'L', 'F',
									'*', 'C', 'W', 'C',
									'S', 'S', 'S', 'S',

									'E', 'D', 'E', 'D',
									'V', 'V', 'V', 'V',
									'G', 'G', 'G', 'G',
									'A', 'A', 'A', 'A',

									'Q', 'H', 'Q', 'H',
									'L', 'L', 'L', 'L',
									'R', 'R', 'R', 'R',
									'P', 'P', 'P', 'P'	};

// reverse complement
const char codon2RCaminoacid[64] = {	'F', 'I', 'L', 'V',
										'Y', 'N', 'H', 'D',
										'S', 'T', 'P', 'A',
										'C', 'S', 'R', 'G',

										'L', 'I', 'L', 'V',
										'*', 'K', 'Q', 'E',
										'S', 'T', 'P', 'A',
										'*', 'R', 'R', 'G',

										'F', 'I', 'L', 'V',
										'Y', 'N', 'H', 'D',
										'S', 'T', 'P', 'A',
										'C', 'S', 'R', 'G',

										'L', 'M', 'L', 'V',
										'*', 'K', 'Q', 'E',
										'S', 'T', 'P', 'A',
										'W', 'R', 'R', 'G'	};





map<string, tax_id> * init_NCBI_NAMES(const char * names);
map<tax_id, string> * init_taxid_to_name(const char * names);

void init_NCBI_NODES(const char * nodes, const char * merged, tax_id *& parent_taxid,  rank_t *& ranks, map<string, rank_t> *& rank_to_id, int NCBI_MAX);


void reverseComplementDNA(string & sequence);


int getTaxa(tax_id ncbi_tax_id, tax_id ** buffer_location, tax_id * parent_taxid, int NCBI_MAX);

int get_LCP_Length(tax_id * a, tax_id * b);

rank_t getLowestRank(tax_id taxid, tax_id * parent_taxid, rank_t * ranks, int NCBI_MAX);

tax_id getTaxonAtRank(tax_id taxid, rank_t rank, tax_id * parent_taxid, rank_t * ranks, int NCBI_MAX);

tax_id getLowestCommonTaxID(tax_id * taxonomy_A, tax_id * taxonomy_B, int NCBI_MAX);

int reduceTaxonByRank(tax_id * buffer_taxid, rank_t * ranks, rank_t assignment_rank_id, int NCBI_MAX);

string getTaxonStringByTaxId(tax_id ncbi_tax_id, map<tax_id, string> ** taxid_to_name, tax_id * parent_taxid, string * names_dmp);

set<tax_id> * init_filtered_taxa(vector<tax_id> * filtered_species_vec, rank_t filter_rank , map<string, rank_t> * rank_to_id, tax_id * parent_taxid, rank_t * ranks, int NCBI_MAX);

bool isFilteredTaxon(rank_t filter_rank, tax_id taxid, set<tax_id> ** filtered_taxa, map<string, rank_t> * rank_to_id, tax_id * parent_taxid, rank_t * ranks, int NCBI_MAX);

tax_id taxonname2taxid(string taxon_string, map<string, tax_id> * NameToNCBI);

int mapRankToSimpleRank(rank_t rank);

string getFASTADescription (string & line);



class FASTA_Parser {

	public:

	FILE *fp;
	std::istream * rna_stream;
	boost::iostreams::stream_buffer<boost::iostreams::file_descriptor_source> * fpstream;

	string description_line;
	bool have_read_descr;
	int expected_sequence_length;
	bool getOnlyIdentfier;
	//bool file_closed;

	FASTA_Parser(string command, bool getOnlyIdentfier);
	FASTA_Parser(string file, bool getOnlyIdentfier, string * zcat);
	void start(string command, bool getOnlyIdentfier);

	~FASTA_Parser();
	string * getSequence();
	bool getNextDescriptionLine(string& descr);


};

class Sequence {
        public:
        string description;
        string data_sequence;

        Sequence(string description, string data_sequence);
};


class Sequence_with_taxid: public string {

	public:
	int offset;
	//string * sequence;
	tax_id taxid;

	string substr(size_t pos, size_t len);
	char operator[](size_t i);
	size_t length() const;
	void setAlignmentSequence(string& seq);

};

#endif
