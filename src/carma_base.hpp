
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

#ifndef GUARD_CARMA_BASE
#define GUARD_CARMA_BASE

#include "global.hpp"

#include "PhyloNode.h"
#include "biotools.hpp"
#include "biotools_blast.hpp"
#include "tools.hpp"
#include "pfam.hpp"
#include "matrix.hpp"

using namespace std;

// foward declarations
class Sequence;
class MatchingQuery;
class ReciprocResults;
class ClassificationResult;

class Class_Substitution_Matrix;
class Blast_HSP;
class Match_hit;
class Fastacmd_Class;

//if you change this numbering, change it also in blast.cpp ! There haven beenn problems with forward declaration of enums.
//enum Match_Type { BLASTX=1, BLASTN=2, BLASTP=3, TBLASTX=4, TBLASTN=5, PFAM_HMM=6, BLAST=7, UNKNOWN=8};
enum Algo_Type { LCA=1, REPR=2, CARMA3=4, BESTBLAST=5};


class CARMA_BASE {
	public:

	std::map<std::string, std::string> * parameters;

	tax_id * parent_taxid;
	rank_t * ranks;
	map<string, rank_t> * rank_to_id;
	map<string, tax_id> * NameToNCBI;
	map<tax_id, string> * taxid_to_name;
	map <string, string> * pfam_to_GOid;

	int score_match;
	int score_mismatch;
	int score_gapopen;
	int score_gapextension;

	map<string, string * > * blast_egts_map;

	int * nucleobase_ascii2num;

	char * path;
	int max_blast_description_length;

	tax_id * buffer_taxid;


	string * pfamA_txt_file;
	string * gene_ontology_txt_file;



	string * blosum_file;

	int NCBI_MAX;

	map<string, ifstream::pos_type > * blast_index;

	double query_multiplier;

	double blastx_evalue;
	double blastn_evalue;
	double blastp_evalue;

	double carma_blastx_evalue;
	double carma_blastn_evalue;
	double carma_blastp_evalue;

	double carma_blastx_bitscore;
	double carma_blastp_bitscore;
	double carma_blastn_bitscore;

	int carma_blastx_alignment_length;
	int carma_blastp_alignment_length;
	int carma_blastn_alignment_length;

//	double blastx_min_bitscore_for_best_hit;
	double blastx_min_bitscore_for_hits;

//	double blastn_min_bitscore_for_best_hit;
	double blastn_min_bitscore_for_hits;

	double CARMA_top_percent;
	double LCA_top_percent;

	double pairwise_blosum_minscore;
	double pairwise_blosum_minoverlap;

	double hmmscan_evalue;

	bool use_hard_threshold;




	string * blast_nr_database;
	string * blast_nt_database;
	string * blast_rdp_database;

	string * blastall_script;
	string * fastacmd_script;
	string * formatdb_script;

	string * nodes_dmp;
	string * merged_dmp;
	string * names_dmp;

	string * pfamId2TaxId_file;

	string * hmmalign_bin;
	string * hmmfetch_bin;
	string * hmmscan_bin;
	string * pfam_A_hmm_file;
	string * pfam_fasta_dir;

	string * config_file;

	string * gzip_bin;
	string * zcat_bin;

	Fastacmd_Class * fastacmd_class;

	Class_Substitution_Matrix * blosum_matrix;

	map<string, tax_id > * pfamid2taxid_map;

	set<tax_id> * filtered_taxa;

	// -------------------------------------------------------

	CARMA_BASE(string * config_file, string * config_overlay_string);
	virtual ~CARMA_BASE();
	//string getTaxonStringByTaxId(tax_id taxid);




	//void classifyBlastXandBlastN(bool use_hard_threshold);
	void classifyBlastXOnly(string * blast_file, string * blast_egts, bool use_hard_threshold, enum Algo_Type algo_type, bool local_temp);
	MatchingQuery * parseBlastEntry(ifstream &myfile, ifstream::pos_type streamposition,  map<string, tax_id> * NameToNCBI, bool use_hard_threshold, bool report_only_first_hsp);





	ClassificationResult * filterMatchingQuery(MatchingQuery * matching_query, enum Algo_Type algo_type);
	void filterTopPercent(MatchingQuery * matching_query, double _top_percent);
	void filterOtherAndUnclassified(MatchingQuery * matching_query);
	void filterByEvalue(MatchingQuery * matching_query, double evalue);
	void filterByBitscore(MatchingQuery * matching_query, double bitscore);
	void filterByAlignmentLength(MatchingQuery * matching_query, int alignment_length);
	void filterReduceTaxa(MatchingQuery * matching_query);

	map<string, string * > * parseStockholm(const char * stockholm_file);




	map<string, string * > * readFastaFileIntoMap(string * fasta_file);




	int computeBlosumScore(string * seq_A, string * seq_B, int s, int e);

	ClassificationResult * computeClassification(MatchingQuery * mq, bool doReduceTaxon);

	void printClassificationResult(ClassificationResult * cl);

	ReciprocResults * doReciprocalStuff(MatchingQuery * matching_query, int match_type);
	ReciprocResults * doReciprocalStuffOnCompleteBlastDB(MatchingQuery * matching_query);

	ClassificationResult * coreClassificationAlgorithm(string query_name, string evalue_str, Blast_HSP * blast_hsp, ReciprocResults * recRes, bool doReduceTaxon);

	//string extractSequenceFromBlastDB(Match_hit * match_hit, Blast_HSP * blast_hsp, enum Match_Type match_type, int best_query_start, int best_query_end);



	void putQuestionmarks(string * sequence);


	ReciprocResults * parseReciprocBLASTP(tax_id best_hit_tax_id, rank_t * ranks, char * blastp_results_file_name, int match_type);
	ReciprocResults * parseReciprocBLASTNP(tax_id best_hit_tax_id, rank_t * ranks, char * reciprok_blast_results_file_name, int match_type);

	tax_id computeLCA(vector<Match_hit * > * hits);


	void setAlignmentParameters(Blast_HSP * fake_blast_hsp, string * seq_query, string * seq_best, int start, int end);

	string translateDNAFragment(int frame_offset, string& dna_sequence);

	string MatchTypeToString(int match_type);

	void readConfiguration(string * config_file);
	string getConfigValue(std::map<std::string, std::string> * parameters, const char * label);
	string * getConfigStringValue(std::map<std::string, std::string> * parameters, const char * label);
	double getConfigDoubleValue(std::map<std::string, std::string> * parameters, const char * label);
	int getConfigIntValue(std::map<std::string, std::string> * parameters, const char * label);
	bool getConfigBoolValue(std::map<std::string, std::string> * parameters, const char * label);
	void setBaseVariables();

};


class ReciprocResults {
	public:
	double query_bitscore;
	double self_bitscore;
	int match_type;

	tax_id best_hit_taxid;

	vector<pair<double, rank_t> * > * other_hits;

	ReciprocResults()
		: query_bitscore(-1), self_bitscore(-1), other_hits(0) {};

	~ReciprocResults(){
		if (this->other_hits!=0) {
			DeleteContainerWithPointers(this->other_hits);
		}
	}

};

class MatchingQuery {
    public:

	int match_type;

	double best_hit_bitscore;
	int best_length;
	string query;
	string * query_sequence;
	long position_in_blast_file;
	bool top_percent_exhausted;

    vector<Match_hit * > * hits;




	MatchingQuery()
		: best_hit_bitscore(-10000), hits(0), best_length(-1), query_sequence(0), position_in_blast_file(-1), top_percent_exhausted(false) {};
    ~MatchingQuery(){
		DeleteContainerWithPointers(hits);
	};

};


class Match_hit {
    public:
        tax_id ncbi_tax_id;
		string * blast_id;
		//double best_hsp_bitscore;
		vector<Blast_HSP * > * HSPs;
		bool masked;
		bool tail;

		Match_hit()
			: ncbi_tax_id(0), blast_id(0), HSPs(0), masked(false), tail(false) {};

		~Match_hit() {
			if (blast_id != 0) {
					delete blast_id;
			}
			DeleteContainerWithPointers(HSPs);
			HSPs=0;
		};
};





class Blast_HSP {
	public:
		double bitscore;
		double evalue;

		int positives;
		int identities;
		int percent_identities;
		int alignment_length;
		int gaps;
		//int gapopen;
		//int mismatches;

		string * query_sequence;
		string * subject_sequence;

		int query_start;
		int query_end;
		int subject_start;
		int subject_end;

		int frame;

		int query_frameshifts;

		Blast_HSP()
			: bitscore(0), evalue(0), positives(-1), identities(-1), percent_identities(-1), alignment_length(-1), gaps(-1), query_sequence(0), subject_sequence(0), query_start(-1), query_end(-1), subject_start(-1), subject_end(-1), frame(0), query_frameshifts(0) {};

		~Blast_HSP(){
			if (query_sequence != 0) {
					delete query_sequence;
			}
			if (subject_sequence != 0) {
					delete subject_sequence;
			}
		};
};


class ClassificationResult {
	public:
	string query;
	string property;
	string evalue;

	tax_id classification_tax_id;
	rank_t classification_rank;

	bool error_no_worse_hit;

	ClassificationResult(){
		classification_tax_id = 0;
		classification_rank = 0;
		error_no_worse_hit = false;
	}
	ClassificationResult(string query, string property, tax_id classification_tax_id, rank_t classification_rank, string evalue)
		: query(query), property(property), classification_tax_id(classification_tax_id), classification_rank(classification_rank), evalue(evalue), error_no_worse_hit(false) {};

};

class DP_Matrix {
	public:
	int dp_matrix_size;
	int ** dp_matrix;
	short int * score_matrix;
	int indel_score;
	int max_diag;

	DP_Matrix(int dp_matrix_size);
	//void init(string& sequence_A, string& sequence_B);
	int compute(string * sequence_A, string * sequence_B);
};

int reduceTaxon(Blast_HSP * blast_hsp, tax_id * buffer_taxid, CARMA_BASE * basic);

short int * readBLOSUM(const char * blosum_file);



#endif
