
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

#ifndef GUARD_GLOBAL
#define GUARD_GLOBAL


//#define NDEBUG


#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <signal.h>

#include <map>
#include <vector>
#include <list>
#include <algorithm>
#include <cctype>
#include "limits.h"
#include <limits>
#include "float.h"
#include <sys/stat.h>
#include <cmath>
#include <stdexcept>
#include <libgen.h>
#include <set>
#include <assert.h>
#include <getopt.h>
#include <sys/statvfs.h>
#include <fcntl.h>


//#include <libconfig.h++>

#include <boost/filesystem.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>

#include <boost/config.hpp>
#include <boost/program_options/detail/config_file.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>




#define HIT_MIN_SUPPORT 1    // minimum hits a read must have before it is classified.


#define PRINT_UNCLASSIFIED_READS 1     // default 1

#define PRINT_UNPARSABLE_BLAST_DESCRIPTIONS 0  // default 0

#define PRINT_UNKNOWN_EMBL_CDS_SPECIES 0  // default 0
#define PRINT_UNKNOWN_SWISSPROT 0  // default 0

#define PRINT_EXTERNAL_CALLS 0

#define ALGORITHM_TYPE 4


 // this was done because forward declarations are not possible with enum
#define MATCH_TYPE_BLASTX 1
#define MATCH_TYPE_BLASTN 2
#define MATCH_TYPE_BLASTP 3
#define MATCH_TYPE_TBLASTX 4
#define MATCH_TYPE_TBLASTN 5
#define MATCH_TYPE_PFAM_HMM 6
#define MATCH_TYPE_BLAST 7
#define MATCH_TYPE_RDP 8
#define MATCH_TYPE_UNKNOWN 9



typedef unsigned char uchar;


#endif
