
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

#ifndef GUARD_CARMA_BLASTNP
#define GUARD_CARMA_BLASTNP


#include "carma_base.hpp"
#include "biotools.hpp"
#include "biotools_rdp.hpp"


using namespace std;

class CARMA_BASE;


class CARMA_BLASTNP : public CARMA_BASE {

	private:



	string * rdp_unaligned;


//	map<string, tax_id> * rdp2taxid_mapping;

	public:
	map<string, string * > * query_sequences;
	map<int,  pair<tax_id, string * > * > * blast_dump;

	CARMA_BLASTNP(string * config_file);
//	CARMA_BLASTNP(string * config_file) : CARMA_BASE(config_file){
//		this->query_sequences=0;
//		rdp_unaligned=0;
//		rdp2taxid_mapping = 0;
//		setVariables();
//	}

	~CARMA_BLASTNP();
	//void reads_query_sequences(string * dna_input_file);
	void setVariables();
	void parse_blast_M9(string * blast_file, string * fasta_input_file, char database, int match_type);
	void classify_blastnp_results(vector<vector<string> * > * blast_hits, string * blast_database, char database, int match_type);
	MatchingQuery * blastn_m8_to_CARMA3(vector<vector<string> * > * blast_hits, string * blast_database, char database, int match_type);
	map<string, tax_id> * get_rdp2taxid_mapping(string * file);
	map<int,  pair<tax_id, string * > * > *  initBlastDump();

};


#endif
