
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

#ifndef GUARD_BIOTOOLS_BLAST
#define GUARD_BIOTOOLS_BLAST


#include "global.hpp"
#include "PhyloNode.h"
#include "carma_base.hpp"

using namespace std;


//enum class Match_Type:short { BLASTX, BLASTN, BLASTP, TBLASTX, TBLASTN, PFAM_HMM};
//enum class Match_Type: short;
//enum Match_Type: short {};
class Blast_HSP;
class Match_hit;

Match_hit * getFirstUnmaskedHit(vector<Match_hit * > * hits);

string * extractSequenceID(string& complete_line);
int extractGiFromBlastDB(string * blast_id, map<string, tax_id> * NameToNCBI, int match_type, string * blast_database, string * fastacmd_script);

tax_id extractTaxonomyInformationById(string * id, int match_type, string * blast_database, string * fastacmd_script);
tax_id extractTaxonomyInformation(string& complete_line, map<string, tax_id> * NameToNCBI, bool use_fastacmd, int match_type, string * blast_database, string * fastacmd_script);

double extractBLASTBitScore(string& line);
double extractBLASTEvalue(string& line);
size_t extractIdentities(string& line);
size_t extractAlignmentLength(string& line);
size_t extractPercentIdentities(string& line);
size_t extractPositives(string& line);
size_t extractGaps(string& line);

map<string, ifstream::pos_type > * indexBlastFile(const char * blastx_file);


int extractFrame(string& line);
int extractSequenceStartNumber(string& line);
int extractSequenceEndNumber(string& line);
string extractSequence(string& line);
int getFrameShiftNum(string& line);

string extractSequenceFromBlastDB(Match_hit * match_hit, Blast_HSP * blast_hsp, int match_type, int best_query_start, int best_query_end, string * blast_nr_database, string * blast_nt_database, char * path, int max_blast_description_length, string * fastacmd_script);
string extractExactSubstringFromBlastDB(string * sequence_id, int start, int end,  int match_type, string * blast_database, string * fastacmd_script);

int getQueryScore(vector<string> * line_data, int score_match, int score_mismatch, int score_gapopen, int score_gapextension);

int runBLAST(string * blastall_script, string * blast_database, double evalue, string * input, string temp_output_file, char blast_type, bool use_gzip, string * gzip_bin);



class BLAST_M9_Parser {

	public:

	FILE *fp;
	std::istream * blast_stream;
	boost::iostreams::stream_buffer<boost::iostreams::file_descriptor_source> * fpstream;

	string current_query;
	vector<string> * first_line_data;
	bool end_of_file;

	BLAST_M9_Parser(string file, string * zcat_bin);

	bool getNextBLASTHits(vector<vector<string> * > ** blast_hits, string& query_id);
	~BLAST_M9_Parser();

};


class Fastacmd_Class {
	private:
	map<string, pair<tax_id, string * > > * data;
	int count_new;
	int count_old;
	public:
	Fastacmd_Class();
	~Fastacmd_Class();

	string fastacmd_sequence(string id, int start, int end, string * blast_database, string * fastacmd_script);

};

#endif
