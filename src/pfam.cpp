
#include "pfam.hpp"


map<string, tax_id > * initPfamId2TaxId(string * pfamId2TaxId_file, string * zcat_bin){


	map<string, tax_id > * pfamid2taxid_map = new map<string, tax_id >();

	FILE *fp;
	int status;
	int max_line_length = 1000;
	char path[max_line_length+1];

	string zcat_command = *zcat_bin;
	zcat_command.append(" -f ");
	zcat_command.append(*pfamId2TaxId_file);

	fp = popen(zcat_command.c_str(), "r");
	if (fp == NULL) {
		cerr << "error" << endl;
		exit(EXIT_FAILURE);
	}

	string line;
	string name;
	string tax_id_string;
	int lines = 0;
	while (fgets(path, max_line_length, fp) != NULL) {
			line = path;
			lines++;
			int line_len = line.length();
			if (line_len > 0 ) {
					//cout << "c:" << (int)line[line_len-1] << endl;

					if (line[line_len-1] == '\n') {
						//cout << "match" << endl;
						//line[line_len-1] == 'X';
						line.erase(line_len-1, 1);
					}
			}

			size_t comma = line.find_first_of(',');
			if (comma == string::npos) {
					cerr << "parsing error, no comma found."<< endl;
					exit(1);
			}

			//cout << line << endl;

			name = line.substr(0, comma);
			tax_id_string = line.substr(comma+1);

			tax_id taxid = (tax_id) str2int(tax_id_string);

			//cout << line.substr(0, comma) << "."<< endl;
			//cout << line.substr(comma+1) << "."<< endl;

			pfamid2taxid_map->insert( pair<string, tax_id >(name, taxid) );

			//exit(0);
	}
	pclose(fp);

	if (lines < 10) {
			cerr << "something wrong! Maybe file does not exist or is empty?" << endl;
			exit(EXIT_FAILURE);
	}

	return pfamid2taxid_map;

}


map<string, string> * init_PfamToGOid (string * pfamA_txt_file, string * gene_ontology_txt_file, string * zcat_bin) {

	map<string, string> * pfam_to_GOid = new map<string, string>;

	map<int, string> * auto_pfamA2pfam_acc = new map<int, string>;

	//string * pfamA_txt_file = new string( "/vol/biodb/pfam24/pfamA.txt.gz");
	//string * gene_ontology_txt_file = new string("/vol/biodb/pfam24/gene_ontology.txt");


	FILE *fp;
	int status;
	int max_line_length = 100000;
	char path[max_line_length+1];

	string zcat_command = *zcat_bin;
	zcat_command.append(" -f ");
	zcat_command.append(*pfamA_txt_file);

	fp = popen(zcat_command.c_str(), "r");
	if (fp == NULL) {
		cerr << "error" << endl;
		exit(EXIT_FAILURE);
	}

	string line;
	string name;
	string tax_id_string;
	int lines = 0;
	while (fgets(path, max_line_length, fp) != NULL) {
			line = path;
			lines++;
			int line_len = line.length();
			if (line_len > 0 ) {
					//cout << "c:" << (int)line[line_len-1] << endl;

					if (line[line_len-1] == '\n') {
						//cout << "match" << endl;
						//line[line_len-1] == 'X';
						line.erase(line_len-1, 1);
					}
			}

			size_t apo1 = line.find_first_of('\'');
			if (apo1 == string::npos) {
					cerr << "parsing error, no apostrophe found."<< endl;
					exit(1);
			}


			size_t apo2 = line.find_first_of('\'', apo1+1);
			if (apo2 == string::npos) {
					cerr << "parsing error, no apostrophe found."<< endl;
					exit(1);
			}


			size_t apo3 = line.find_first_of('\'', apo2+1);
			if (apo3 == string::npos) {
					cerr << "parsing error, no apostrophe found."<< endl;
					exit(1);
			}

			size_t apo4 = line.find_first_of('\'', apo3+1);
			if (apo4 == string::npos) {
					cerr << "parsing error, no apostrophe found."<< endl;
					exit(1);
			}


//			cout << apo1 << endl;
//			cout << apo2 << endl;
//			cout << apo3 << endl;
//			cout << apo4 << endl;

			string auto_pfamA_string = line.substr(apo1+1, apo2-apo1 -1);
			string pfamA_acc = line.substr(apo3+1, apo4-apo3 -1);

			int auto_pfamA_int = str2int(auto_pfamA_string);

			//(*pfam_acc2auto_pfamA)[pfamA_acc] = auto_pfamA_int;
			(*auto_pfamA2pfam_acc)[auto_pfamA_int] = pfamA_acc;
			//cout << auto_pfamA_string << endl;
			//cout << pfamA_acc << endl;

	}
	pclose(fp);

	if (lines < 10) {
			cerr << "something wrong! Maybe file does not exist or is empty?" << endl;
			exit(EXIT_FAILURE);
	}

	// ----------------------------------------------


	zcat_command = *zcat_bin;
	zcat_command.append(" -f ");
	zcat_command.append(*gene_ontology_txt_file);

	fp = popen(zcat_command.c_str(), "r");
	if (fp == NULL) {
		cerr << "error" << endl;
		exit(EXIT_FAILURE);
	}

	//string line;
	//string name;
	//string tax_id_string;
	lines = 0;
	while (fgets(path, max_line_length, fp) != NULL) {
			line = path;
			lines++;
			int line_len = line.length();
			if (line_len > 0 ) {
					//cout << "c:" << (int)line[line_len-1] << endl;

					if (line[line_len-1] == '\n') {
						//cout << "match" << endl;
						//line[line_len-1] == 'X';
						line.erase(line_len-1, 1);
					}
			}

			size_t apo1 = line.find_first_of('\'');
			if (apo1 == string::npos) {
					cerr << "parsing error, no apostrophe found."<< endl;
					exit(1);
			}


			size_t apo2 = line.find_first_of('\'', apo1+1);
			if (apo2 == string::npos) {
					cerr << "parsing error, no apostrophe found."<< endl;
					exit(1);
			}


			size_t apo3 = line.find_first_of('\'', apo2+1);
			if (apo3 == string::npos) {
					cerr << "parsing error, no apostrophe found."<< endl;
					exit(1);
			}

			size_t apo4 = line.find_first_of('\'', apo3+1);
			if (apo4 == string::npos) {
					cerr << "parsing error, no apostrophe found."<< endl;
					exit(1);
			}


//			cout << apo1 << endl;
//			cout << apo2 << endl;
//			cout << apo3 << endl;
//			cout << apo4 << endl;

			string auto_pfamA_string = line.substr(apo1+1, apo2-apo1 -1);
			string go_id = line.substr(apo3+1, apo4-apo3 -1);

			int auto_pfamA_int = str2int(auto_pfamA_string);

			//(*pfam_acc2auto_pfamA)[pfamA_acc] = auto_pfamA_int;

			//cout << auto_pfamA_string << endl;
			//cout <<  go_id << endl;

			map<int, string>::iterator it;

			//(*auto_pfamA2pfam_acc)[auto_pfamA_int]
			it = auto_pfamA2pfam_acc->find(auto_pfamA_int);

			if (it == auto_pfamA2pfam_acc->end()) {
					cerr << "error: auto_pfamA2pfam_acc, element not found" << endl;
					exit(1);
			}
			string pfam_acc = it->second;
			//cout << "pfam: " <<  pfam_acc << endl;

			map<string, string>::iterator pf2go_it;
			pf2go_it = pfam_to_GOid->find(pfam_acc);
			if (pf2go_it == pfam_to_GOid->end()) {
				pfam_to_GOid->insert( pair<string, string>(pfam_acc, go_id));
			} else {
				//cout << "map: " <<  pf2go_it->second << endl;
				pf2go_it->second.append(",");
				pf2go_it->second.append(go_id);
				//cout << "map2: " <<  pf2go_it->second << endl;
			}



	}
	pclose(fp);

	if (lines < 10) {
			cerr << "something wrong! Maybe file does not exist or is empty?" << endl;
			exit(EXIT_FAILURE);
	}



	delete auto_pfamA2pfam_acc;

	return pfam_to_GOid;
}


string hmmfetch(string * hmmfetch_bin, string * pfam_A_hmm_file, string pfam_descr) {

	// hmmfetch example:
	// /home/wgerlach/hmmer/hmmer-3.0rc2-linux-intel-x86_64/binaries/hmmfetch  /vol/biodb/pfam24/Pfam-A.hmm PF02698.10

	char single_hmm_file_name [L_tmpnam];
	char * tmpnam_res1 = tmpnam ( single_hmm_file_name );


	string hmmfetch_command = *(hmmfetch_bin);
	hmmfetch_command.append(" -o ");
	hmmfetch_command.append(single_hmm_file_name);
	hmmfetch_command.append(" ");
	hmmfetch_command.append(*pfam_A_hmm_file);
	hmmfetch_command.append(" ");
	hmmfetch_command.append(pfam_descr);
	hmmfetch_command.append(" > /dev/null"); // otherwise I get useless information..

	if (PRINT_EXTERNAL_CALLS == 1) {
		cerr << "call: " << hmmfetch_command << endl;
	}

	int sys_res = system(hmmfetch_command.c_str());
	if (sys_res != 0) {
			cerr << "error calling hmmfetch" << endl;
			cerr << "call was: " << hmmfetch_command << endl;
			exit(EXIT_FAILURE);
	}

	string single_hmm_file_name_string = single_hmm_file_name;

	return single_hmm_file_name_string;
}


string hmmalign(string * hmmalign_bin, string single_hmm_file_name_string, const char * combined_egt_pfam_file_name ) {

	// hmmalign example:
	//string hmmalign_bin = "/home/wgerlach/hmmer/hmmer-3.0rc2-linux-intel-x86_64/binaries/hmmalign";
	//string hmmalign_bin = "hmmalign-3";

	char new_alignment_file_name [L_tmpnam];
	char * tmpnam_res3 = tmpnam ( new_alignment_file_name );


	string hmmalign_command = *(hmmalign_bin);
	hmmalign_command.append(" --trim --amino --informat FASTA --outformat Stockholm ");
	hmmalign_command.append(" -o ");
	hmmalign_command.append(new_alignment_file_name);
	hmmalign_command.append(" ");
	hmmalign_command.append(single_hmm_file_name_string);
	hmmalign_command.append(" ");
	hmmalign_command.append(combined_egt_pfam_file_name);


	if (PRINT_EXTERNAL_CALLS == 1) {
		cerr << "call: " << hmmalign_command << endl;
	}

	int sys_ret2 = system(hmmalign_command.c_str());
	if (sys_ret2 != 0) {
		cerr << "error calling " << hmmalign_command << endl;
		exit(EXIT_FAILURE);
	}

	string new_alignment_file_name_string = new_alignment_file_name;
	return new_alignment_file_name_string;
}


string hmmscan(string * hmmscan_bin, double hmmscan_evalue, string * pfam_A_hmm_file, char * hmmscan_input_file_name) {

// example: hmmscan-3 /vol/biodb/pfam24/Pfam-A.hmm /tmp/fileDh36vp


	char hmmscan_output_file_name [L_tmpnam];
	char * tmpnam_res2 = tmpnam ( hmmscan_output_file_name );


	string hmmscan_command = *hmmscan_bin;

	hmmscan_command.append(" -E ");
	hmmscan_command.append(stringify(hmmscan_evalue));
	hmmscan_command.append(" --noali --notextw ");
	hmmscan_command.append(" --domtblout ");
	hmmscan_command.append(hmmscan_output_file_name);
	hmmscan_command.append(" -o /dev/null ");
	hmmscan_command.append(*pfam_A_hmm_file);
	hmmscan_command.append(" ");
	hmmscan_command.append(hmmscan_input_file_name);


	if (PRINT_EXTERNAL_CALLS == 1) {
		cerr << "call: " << hmmscan_command << endl;
	}

	int ret_hmmscan = system(hmmscan_command.c_str());

	if (ret_hmmscan != 0) {
		cerr << "error calling: " << hmmscan_command << endl;
		exit(1);
	}


	string hmmscan_output_file_name_string = hmmscan_output_file_name;
	return hmmscan_output_file_name_string;
}
