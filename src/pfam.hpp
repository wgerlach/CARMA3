#ifndef GUARD_PFAM
#define GUARD_PFAM




#include "tools.hpp"
#include "PhyloNode.h"

using namespace std;


map<string, tax_id > * initPfamId2TaxId(string * pfamId2TaxId_file, string * zcat_bin);
map<string, string> * init_PfamToGOid (string * pfamA_txt_file, string * gene_ontology_txt_file, string * zcat_bin);


string hmmfetch(string * hmmfetch_bin, string * pfam_A_hmm_file, string pfam_descr);
string hmmalign(string * hmmalign_bin, string single_hmm_file_name_string, const char * combined_egt_pfam_file_name );
string hmmscan(string * hmmscan_bin, double hmmscan_evalue, string * pfam_A_hmm_file, char * hmmscan_input_file_name);


#endif
