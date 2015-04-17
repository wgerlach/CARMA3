
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

#ifndef GUARD_CARMA_HMMER
#define GUARD_CARMA_HMMER


#include "carma_base.hpp"
#include "biotools.hpp"
//#include "biotools_rdp.hpp"


using namespace std;

class CARMA_BASE;
class HmmscanDomainMatch;


class CARMA_HMMER : public CARMA_BASE {

	private:
	string * blast_egts;

	public:

	CARMA_HMMER(string * config_file) : CARMA_BASE(config_file, 0){
		this->blast_egts=0;
	}


	void processFastaFile(string * dna_fasta_file, string * blast_egts, bool is_dna);
	void processDNAFragment(string& descr, string& sequence, bool is_dna);
	vector<HmmscanDomainMatch * > * parse_hmmscan_results(const char * hmmscan_output_file_name);
	vector<string> * getEGTDescription(string& descr);
	void insertSequence(string& description, string& dna_sequence, map<string, vector<Sequence * > * > * sequencesByFamily);
	map<string, vector<Sequence * > * > * parseEGTFile(string * egt_file);
	void processEGTs(map<string, vector<Sequence * > * > * sequencesByFamily);
	void processEGTfamily(string pfam_descr, vector<Sequence * > * vec);

};



class HmmscanDomainMatch {
		public:
		string accession;
		int frame;
		double full_seq_evalue;
		double i_evalue;
		int ali_from;
		int ali_to;
};




#endif
