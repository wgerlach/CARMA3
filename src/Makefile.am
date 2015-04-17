bin_PROGRAMS = carma

carma_SOURCES = main.cpp main.h global.hpp PhyloNode.cpp PhyloNode.h matrix.cpp matrix.hpp biotools.cpp biotools.hpp tools.cpp tools.hpp pfam.hpp pfam.cpp biotools_blast.cpp biotools_blast.hpp biotools_rdp.cpp biotools_rdp.hpp carma_rna.cpp carma_rna.hpp carma_base.cpp  carma_base.hpp carma_blastnp.cpp carma_blastnp.hpp carma_hmmer.cpp carma_hmmer.hpp

carma_CXXFLAGS = -m64 $(BOOST_CPPFLAGS)
carma_LDFLAGS = -m64 $(BOOST_LDFLAGS)

carma_LDADD = $(BOOST_IOSTREAMS_LIB) $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB)
