
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


#include "main.h"




using namespace std;

string * local_stderr_file;
string * final_stderr_file;

void preExit ()
{

	fflush (stderr);
	fclose (stderr);

	if (local_stderr_file != 0) {
		copyFile(local_stderr_file->c_str(), final_stderr_file->c_str());
	}


}



void usage() {
		cout << endl;
		cout << "**** CARMA3 ****" << endl;
		cout << endl;
		cout << "Usage examples: " << endl;
		cout << "    carma --blast --type n --database p --input dna.fas --output dna.blastx" << endl;
		cout << "    carma --classify-blast --type n --database p --input dna.blastx --output dna.tax" << endl;
		cout << endl;
		cout << endl;
}

 void sighandler ( int signum )
 {
 printf ( "SIGINT received. Inside sighandler\n" ) ;
 }



int set_cloexec_flag (int desc, int value)
     {
       int oldflags = fcntl (desc, F_GETFD, 0);
       /* If reading the flags failed, return error indication now. */
       if (oldflags < 0)
         return oldflags;
       /* Set just the flag we want to set. */
       if (value != 0)
         oldflags |= FD_CLOEXEC;
       else
         oldflags &= ~FD_CLOEXEC;
       /* Store modified flag word in the descriptor. */
       return fcntl (desc, F_SETFD, oldflags);
     }


int main(int argc, char* argv[])
{
	local_stderr_file = 0;
	final_stderr_file = 0;
	//signal ( SIGINT,  sighandler ) ;
	//signal ( SIGKILL,  sighandler ) ;
	//signal ( SIGHUP,  sighandler ) ;
	//signal ( SIGTERM,  sighandler ) ;


	// ---------------------------------------- parameter check --------------------------



	FILE * blabla;
	string temp_output_file="";

	string * blast_egts_final = 0;

	string * input_file = 0;
	string * fasta_input_file = 0;
	string * error_file = 0;
	string * output_file = 0;
	string * config_file = new string("carma.cfg");

	string * config_overlay_string = 0;

	vector<int> * taxid_vec = 0;

	char database='-';
	char input_type = 'p';
	bool use_gzip = false;
	bool local = false;
	int delay = 0;
	namespace po = boost::program_options;


	enum EnumMode 				{ Ehmmer,  Eblast,   Eclassifyhmmer,  Eclassifyblast,  Eclassifyrdp,   EcreatespecificRDPseq,     EcreatespecificRDPali    ,  Erdptraindata };
	const char* const modes[] = {"hmmer", "blast",  "classify-egt",  "classify-blast", "classify-rdp", "create_specific_RDP_seq", "create_specific_RDP_ali", "rdp_data"};
	int opt_mode = -1;
	//std::vector<std::string> v = {"Hello", "World"};





	// Declare the supported options.
	po::options_description homology("Homology search modes");
	homology.add_options()
		(modes[0], 	"HMMER search against Pfam\n"
					"  use with --type (input: fasta, output: EGTs)")
		(modes[1],	"blast DNA/Protein/RNA against NT/NR/RDP\n"
					"  use with --type and --database (input: fasta)\n"
					"  optional for blastx: --blast-egts\n"
					"  optional for blastp: --gzip")
	//    ("help", "produce help message")
	//    ("compression", po::value<int>(), "set compression level")
	;


	po::options_description classification("Classification modes");
	classification.add_options()
		(modes[2], 	"classify HMMER results (input: EGTs)\n"
					"  optional: --blast-egts")
		(modes[3],	"classify BLAST results(input: BLAST)\n"
					"  use with --type\n"
		            "  BLASTx results must be in default output format\n"
		            "  BLASTn/p results must be in m8/9 output format\n"
					"  BLASTn/p: use with --fasta-input\n"
					"  BLASTn: use with --database")
		(modes[4], 	"classify RDP-BLAST results(input: BLAST m8/9)")
	;


	po::options_description optional("Optional");
	optional.add_options()
		("type", po::value<char>(), 			"input type\n"
												"  n - nucleotide (DNA/RNA)\n"
												"  p - protein")
		("blast-egts", po::value<string>(), 	"input/output file for BLAST based EGTs\n"
												"  optional for --classify-blastx and --search-hmmer")
		("fasta-input", po::value<string>(), 	"additional input file for BLASTn/p based classification\n"
												"  required for --classify-blast in case of blastn and\n  blastp")
		("database", po::value<char>(), 		"choose database to blast against, possible values:\n"
												"  n - NCBI nucleotide sequence database NT\n"
												"  p - NCBI protein sequence database NR\n"
												"  r - Ribosomal Database Project\n"
												"  required for --blast and --classify-blast")
		("local",								"write to local machine, use with --output")
		("gzip",                        		"use gzip")
	;

	po::options_description preprocessing("Preprocessing modes");
	preprocessing.add_options()
		(modes[5], 	"preprocessing") //--create_specific_RDP_seq
		(modes[6], 	"preprocessing") //--create_specific_RDP_ali
		(modes[7], "experimental") // rdp train
	//	("taxids", po::value< vector<int> >(), "experimental_t")
		("taxids", po::value<vector<int> >()->multitoken(), "experimental")
		("filter_rank", po::value<string>(), "experimental")
//value<vector<float> >()->multitoken()

	;

	po::options_description general("General");
	general.add_options()
		("input",	po::value<string>(),								"input file (required)")
		("output",	po::value<string>(),								"output file (default STDOUT)")
		("error",	po::value<string>(),								"error file (default STDERR)")
		("config",	po::value<string>()->default_value("carma.cfg"),	"configuration file")
		("config_overlay",	po::value<string>(),						"overwrite config (e.g. carma_blastx_evalue=0.1)\n"
																		"  comma-separated, experimental")
		("delay",	po::value<int>(),									"delay of n*3 sec")
		("help",														"produce help message")
	//    ("compression", po::value<int>(), "set compression level")
	;




	po::options_description all("Allowed options");
	all.add(homology).add(classification).add(optional).add(preprocessing).add(general);

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, all), vm);
	po::notify(vm);


	// ---------------------------------------- parameter syntax check -------------------------

	if (vm.count("output")) {
		output_file  = new string(vm["output"].as<string>() );
	}

	if (vm.count("delay")){

		delay = vm["delay"].as<int>();
		int delay_work = delay%500;
		delay_work++;
		delay_work=delay_work*2;


		int delay_check = delay%50;
		delay_check++;
		sleep(delay_check);
		if (output_file != 0) {
			if (FileExists((char *) output_file->c_str())) {
				if (not local) {
					cerr << "File \"" << *output_file << "\" already exists. Skip this job." << endl;
				}
				exit(0);
			}
		}

		sleep(delay_work);

	}

	if (vm.count("local")) {
		local = true;
	}

	if (vm.count("error")) {
		error_file = new string(vm["error"].as<string>() );

		if (local) {
			final_stderr_file = new string(*error_file);

			FILE * temporary_error_file = tmpfile();
			int temporary_error_file_no=fileno(temporary_error_file);

			string temp_error_file_string = "/proc/self/fd/";
			temp_error_file_string.append(int2str(temporary_error_file_no));

			local_stderr_file = new string (temp_error_file_string);
			blabla =freopen(local_stderr_file->c_str(),"w",stderr);
		} else {

			blabla =freopen(error_file->c_str(),"w",stderr);
		}


		atexit (preExit);

	}



	if (vm.count("help")) {
		usage();
		cout << all << endl;
		exit(0);
	}



	if (vm.count("input")) {
		input_file  = new string(vm["input"].as<string>() );
	} else {
		usage();
		cout << all << endl;
		cerr << "ERROR: --input is missing" << endl;
		exit(1);
	}



	if (vm.count("type")) {
		input_type  = vm["type"].as<char>();
	} else {
		usage();
		cout << all << endl;
		cerr << "ERROR: --type is missing" << endl;
		exit(1);
	}

	if (vm.count("config")) {
		if (config_file != 0) {
			delete config_file;
		}
		config_file = new string(vm["config"].as<string>() );
	}



	if (vm.count("config_overlay")) {
		config_overlay_string = new string(vm["config_overlay"].as<string>() );
	}


	if (vm.count("fasta-input")) {
		fasta_input_file = new string(vm["fasta-input"].as<string>() );
	}

	for (int i=0; i<=7; i++) {
		if (vm.count(modes[i])) {
			if (opt_mode != -1) {
				cerr << "error: you can use only one mode at once (--" << modes[i] << " , --" << modes[opt_mode] << ")" <<endl;
			}
			opt_mode = i;
		}
	}

	if (opt_mode == -1) {
		usage();
		cout << all << endl;
		exit(1);
	}

	if (vm.count("database")) {
		database = vm["database"].as<char>();
		if ((database != 'n') && (database != 'p') && (database != 'r')) {
			cerr << "error: --database argument " << database << " unknown" << endl;
			exit(1);
		}

	} else {
		if (opt_mode == Eblast || opt_mode == Eclassifyblast) {
			cerr << "error: argument --database missing" << endl;
			exit(1);
		}

	}

	if (vm.count("blast-egts")) {
		blast_egts_final = new string(vm["blast-egts"].as<string>() );
	}

	if (vm.count("gzip")) {
		use_gzip = true;
	}



	if (vm.count("taxids"))
	{

		vector<int> files = vm["taxids"].as< vector<int> >();
		taxid_vec = new vector<int>(files.begin(), files.end());

	//            for(int i=0; i<files.size(); ++i) {
	//                cout << files[i]<<",";
	//            }
	//            cout << "\n";
	//		exit(1);
	}
	//exit(1);

	string filter_rank = "";
	if (vm.count("filter_rank")) {
		filter_rank = vm["filter_rank"].as<string>();
	}
// ---------------------------------------- semantic check --------------------------

	cerr << "argv: " ;
	for (int i = 0; i<argc; i++) {
		cerr << argv[i] << " ";
	}
	cerr << endl;


	if (local) {
		// copy input_file

		FILE * temporary_input_file = tmpfile();
		int temporary_input_file_no=fileno(temporary_input_file);

		string temp_input_file_string = "/proc/self/fd/";
		temp_input_file_string.append(int2str(temporary_input_file_no));

		
		if (not FileExists((char *) temp_input_file_string.c_str())) {
			cerr << "error: could not open temp_input_file " << temp_input_file_string << endl;
			if (vm.count("error")) {
				preExit();
			}
			exit(1);
		}
		//cerr << "got1: " << *input_file << endl;
		//cerr << "got1: " << temp_input_file_string << endl;
		bool copy_ret = copyFile(input_file->c_str(), temp_input_file_string.c_str());
		if (copy_ret == false) {
			cerr << "error copying file from " << *input_file << " to "  << temp_input_file_string << endl;
			exit(1);
		}

		input_file->assign(temp_input_file_string);

	}



//	if (local && (blast_egts_final != 0)) {
//
//
//		FILE * temporary_blastegts_file = tmpfile();
//		int temporary_blastegts_file_no=fileno(temporary_blastegts_file);
//		string temp_blastegts_file_string = "/proc/self/fd/";
//		temp_blastegts_file_string.append(int2str(temporary_blastegts_file_no));
//		blast_egts_temp = new string (temp_blastegts_file_string);
//
//		// for (opt_mode == Eclassifyblast)  copy in other direction
//		if (opt_mode == Ehmmer) {
//
//			bool copy_ret = copyFile(blast_egts_final->c_str(), blast_egts_temp->c_str());
//			if (copy_ret == false) {
//				cerr << "error copying file" << endl;
//				exit(1);
//			}
//
//		}
//	}

	if (output_file != 0) {

//		char * outfile_dir = new char[255];
//		strcpy(outfile_dir,  output_file->c_str());
//
//		outfile_dir = dirname(outfile_dir);
//		size_t disk_space = getFreeDiskSpace(outfile_dir); // bytes
//		disk_space /= (1024*1024); // MB
//
//		if (disk_space < 100) {
//			cerr << "Warning: Only " << disk_space << " MB free at mountpoint \"" << outfile_dir << "\" !" << endl;
//			flush(cerr);
//		}

		if (FileExists((char *) output_file->c_str())) {
				cerr << "File \"" << *output_file << "\" already exists. Skip this job." << endl;
				exit(0);
		}

		temp_output_file = *output_file;
		temp_output_file.append(".part");
		if (FileExists((char *) temp_output_file.c_str())) {
			remove(temp_output_file.c_str());
		}

		if (local) {



			FILE * temporary_output_file = tmpfile();
			int filedesc=fileno(temporary_output_file);
			//temp_output_file = "/dev/fd/";
			temp_output_file = "/proc/self/fd/";
			temp_output_file.append(int2str(filedesc));
			//cerr << "got2: " << temp_output_file << endl;
			//cerr << "file: " << temp_output_file << endl;
		}

		if (opt_mode != Eblast) {
			blabla =freopen(temp_output_file.c_str(),"w",stdout);
		}
	}



	if (opt_mode==Eclassifyblast && ((input_type=='n' && database != 'p') || input_type=='p') )  {
		if (fasta_input_file == 0) {
			cerr << "error: --fasta-input is missing" << endl;
			exit(1);
		}

		if (database == '-' && (opt_mode==Eblast || opt_mode==Eclassifyblast)) {
			cerr << "error: --database is missing or wrong" << endl;
			exit(1);
		}
	}



// ----------------------------------------  --------------------------


	CARMA_BASE   * basic = 0;
	CARMA_RNA    * basic_rna = 0;
	CARMA_BLASTNP * basic_blastnp = 0;
	CARMA_HMMER  * basic_hmmer = 0;

	if ( ((opt_mode==Eblast) && (database== 'r') ) || opt_mode == Eclassifyrdp || opt_mode == EcreatespecificRDPseq || opt_mode == EcreatespecificRDPali || opt_mode == Erdptraindata) {
		basic = new CARMA_RNA(config_file);
		basic_rna = (CARMA_RNA *) basic;
	} else if (opt_mode==Eclassifyblast && (not (input_type=='n' && database=='p'))) {
		basic = new CARMA_BLASTNP(config_file);
		basic_blastnp = (CARMA_BLASTNP *) basic;
	} else if (opt_mode == Ehmmer || opt_mode == Eclassifyhmmer) {
		basic = new CARMA_HMMER(config_file);
		basic_hmmer = (CARMA_HMMER *) basic;
	} else {
		basic = new CARMA_BASE(config_file, config_overlay_string);
	}


// ---------------------------------------- pre-processing --------------------------




if ( opt_mode != Eblast && opt_mode != Ehmmer ) {
	cerr << "read NCBI data (nodes.dmp, merged.dmp)...";
	flush(cerr);
	init_NCBI_NODES(basic->nodes_dmp->c_str(), basic->merged_dmp->c_str(), basic->parent_taxid, basic->ranks, basic->rank_to_id, basic->NCBI_MAX);
	cerr << "   done." << endl << endl;
	flush(cerr);


	cerr << "read NCBI data (nodes.dmp)...";
	flush(cerr);
	basic->NameToNCBI = init_NCBI_NAMES(basic->names_dmp->c_str());
	cerr << "   done." << endl << endl;
	flush(cerr);
}



// -----------------------------------------------------------------------------------


	switch (opt_mode) {
		case Eclassifyrdp:
			{
				//CARMA_RNA * basic_rna = (CARMA_RNA *) basic;

				basic_rna->parse_16S_blast(*input_file);

				break;
			}
		case EcreatespecificRDPseq:
			{
				cerr << "parse RDP Sequences file..." << endl;
				basic_rna->extractSpecificRDP(*input_file, 2000, false, false, taxid_vec, filter_rank);
				delete taxid_vec;
				break;
			}
		case EcreatespecificRDPali:
			{
				cerr << "parse RDP Alignment file..." << endl;
				basic_rna->extractSpecificRDP(*input_file, 26000, true, true, taxid_vec, filter_rank);
				delete taxid_vec;
				break;
			}
		case Ehmmer:
			{
				bool is_dna;
				if (input_type=='n') {
					is_dna= true;
				} else {
					is_dna= false;
				}
				basic_hmmer->processFastaFile(input_file, blast_egts_final, is_dna);
				break;
			}
		case Eclassifyhmmer:
			{
				cerr << "read "<< *(basic_hmmer->pfamId2TaxId_file) << "...";
				flush(cerr);
				basic_hmmer->pfamid2taxid_map = initPfamId2TaxId(basic_hmmer->pfamId2TaxId_file, basic_hmmer->zcat_bin);
				cerr << "   done." << endl << endl;
				flush(cerr);


				map<string, vector<Sequence * > * > * sequencesByFamily = basic_hmmer->parseEGTFile(input_file);

				basic_hmmer->processEGTs(sequencesByFamily);
				break;
			}
		case 101:
			{
				// best BLAST hit
				cerr << "index BLAST file..." << endl;
				flush(cerr);
				// index blastx
				basic->blast_index = indexBlastFile(input_file->c_str());
				cerr << "   done." << endl << endl;
				flush(cerr);

				cerr << "parse BLAST and classify reads..." << endl;
				flush(cerr);


				cerr << "   ...." << endl;
				basic->classifyBlastXOnly(input_file, blast_egts_final, basic->use_hard_threshold, BESTBLAST, local);
				cerr << "   done." << endl << endl;
				flush(cerr);
				break;
			}
		case 102:
			{
				// LCA
				cerr << "index BLAST file..." << endl;
				flush(cerr);
				basic->blast_index = indexBlastFile(input_file->c_str());
				cerr << "   done." << endl << endl;
				flush(cerr);

				cerr << "parse BLAST and classify reads..." << endl;
				flush(cerr);


				cerr << "   ...." << endl;
				basic->classifyBlastXOnly(input_file, blast_egts_final, basic->use_hard_threshold, LCA, local);
				cerr << "   done." << endl << endl;
				flush(cerr);
				break;
			}
		case 103:
			{
				// F-measure...
				break;
			}
		case 104:
			{


				cerr << "index BLASTX file..." << endl;
				flush(cerr);
				// index blastx
				basic->blast_index = indexBlastFile(input_file->c_str());
				cerr << "   done." << endl << endl;
				flush(cerr);


				cerr << "index BLASTN file..." << endl;
				flush(cerr);
				// index blastx
				basic->blast_index = indexBlastFile(input_file->c_str());
				cerr << "   done." << endl << endl;
				flush(cerr);




				cerr << "parse BLAST and classify reads..." << endl;
				flush(cerr);

				//map<string, MatchingQuery * > * matching_queries_blastx = basic->parseBlast(blastx_file->c_str(), NameToNCBI, blastx_min_bitscore_for_hits, top_percent);
				cerr << "   ...." << endl;
				//basic->classifyBlastXandBlastN(basic->use_hard_threshold);
				//map<string, MatchingQuery * > * matching_queries_blastn = basic->parseBlast();
				cerr << "   done." << endl << endl;
				flush(cerr);

				//cerr << "merge..." << endl;
				//flush(cerr);
				//map<string, MatchingQuery * > * matching_queries = mergeHits(matching_queries_blastx, matching_queries_blastn);
				//cerr << "   done." << endl << endl;
				//flush(cerr);


				//DeleteMapWithPointers(matching_queries_blastx);
				//DeleteMapWithPointers(matching_queries_blastn);

				//cerr << "compute classification..." << endl;
				//flush(cerr);
				//computeClassification(matching_queries_blastn, parent_taxid, 50, 0.1);
				//basic->computeClassification(matching_queries_blastx, blastx_min_bitscore_for_best_hit, top_percent);
				//cerr << "   done." << endl << endl;
				//flush(cerr);
				break;
			}
		case Eclassifyblast:
			{
				if (input_type=='n' && database == 'p' ) {
					// classification from blastx output
					cerr << "index BLASTX file..." << endl;
					flush(cerr);
					// index blastx
					basic->blast_index = indexBlastFile(input_file->c_str());
					cerr << "   done." << endl << endl;
					flush(cerr);


					cerr << "parse BLAST and classify reads..." << endl;
					flush(cerr);


					cerr << "   ...." << endl;
					basic->classifyBlastXOnly(input_file, blast_egts_final, basic->use_hard_threshold, CARMA3, local);
					cerr << "   done." << endl << endl;
					flush(cerr);
				} else if (input_type=='n' && (database == 'n' || database == 'r')) {
					// classification from blastn (NT or RDP) output
					basic_blastnp->parse_blast_M9(input_file, fasta_input_file, database, MATCH_TYPE_BLASTN);
				} else if (input_type=='p' && database == 'p' ) {
					// classification from blastp output
					basic_blastnp->parse_blast_M9(input_file, fasta_input_file, database, MATCH_TYPE_BLASTP);
				} else {
					cerr << "error: this combination of database and input_type is not supported" << endl;
					exit(1);
				}

				break;
			}
		case Eblast: // blastx
			{
				int blastret = 1;
				if (input_type == 'n' && database == 'p') {
					//blastx NR
					blastret = runBLAST(basic->blastall_script, basic->blast_nr_database, basic->blastx_evalue, input_file, temp_output_file, 'x', false, 0);
				} else if (input_type == 'n' && database == 'n') {
					//blastn NT
					blastret = runBLAST(basic->blastall_script, basic->blast_nt_database, basic->blastn_evalue, input_file, temp_output_file,'n', false, 0);
				} else if (input_type == 'n' && database == 'r') {
					//blastn RDP
					blastret = runBLAST(basic_rna->blastall_script, basic_rna->blast_rdp_database, basic_rna->blast_rdp_evalue, input_file, temp_output_file, 'r', false, 0);
				} else if (input_type == 'p' && database == 'p') {
					//blastp NR
					blastret = runBLAST(basic->blastall_script, basic->blast_nr_database, basic->blastp_evalue, input_file, temp_output_file,'p', use_gzip, basic->gzip_bin);
				} else {
					cerr << "error: this combination of database and input_type is not supported" << endl;
					exit(1);
				}
				
				
				if (blastret != 0) {
					cerr << "error: runBLAST function returned with error" << endl;
					if (vm.count("error")) {
						preExit();
					}
					exit(1);
				}
				
				break;
			}

		case Erdptraindata:
			{
				//CARMA_RNA * basic_rna = (CARMA_RNA *) basic;
				basic_rna->rdp_train_data(input_file);
				break;
			}
		default:
			{
				cerr << "error: should not reach here..." << endl;
				break;
			}
	} // end switch


	if (output_file != 0) {
		if (opt_mode != Eblast) {
			fclose (stdout);
		}
//system("ls -la /dev/fd/");
//		system("cat /dev/fd/3");
//		system("ls -la /dev/fd/");
//		int oldflags = fcntl (filedesc, F_GETFD, 0);

//		cerr << "oldflags: " << oldflags << endl;


		if (local) { // first copy, then rename
			string output_file_part = *output_file;
			output_file_part.append(".part");
			bool copy_ret = copyFile(temp_output_file.c_str(), output_file_part.c_str());
			if (copy_ret == false) {
				cerr << "error copying file from " << temp_output_file << " to " << output_file_part << endl;
				exit(1);
			}
			sleep(2);
			int rename_ret = rename ( output_file_part.c_str(), output_file->c_str());
			if (rename_ret != 0) {
					cerr << "error on renameing file" << endl;
			}
		} else {

			int rename_ret = rename ( temp_output_file.c_str(), output_file->c_str());
			if (rename_ret != 0) {
					cerr << "error on renameing file" << endl;
			}
		}

	}




	if (error_file != 0) {
		fclose (stderr);
		error_file = 0;
	}


	delete basic;


	if (config_file != 0) {
		delete config_file;
	}
	if (input_file != 0) {
		delete input_file;
	}
	if (output_file != 0) {
		delete output_file;
	}

	if (error_file != 0) {
		delete error_file;
	}
	if (blast_egts_final != 0) {
		delete blast_egts_final;
	}

	if (fasta_input_file != 0) {
		delete fasta_input_file;
	}

	return 0;
}
