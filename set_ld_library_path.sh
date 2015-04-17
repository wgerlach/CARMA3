#!/bin/sh





# SGE setup
# modify these variable to select the required SGE version
SGE_ROOT=/vol/codine-6.2
SGE_CELL=default

. $SGE_ROOT/$SGE_CELL/common/settings.sh


# unfortunatly, no LD_LIBRARY_PATH is set in settings.sh
ARCH=`$SGE_ROOT/util/arch`
if [ $ARCH = "sol-sparc" ]; then
    LD_LIBRARY_PATH=$SGE_ROOT/lib/sol-sparc:$LD_LIBRARY_PATH
else
    if [ $ARCH = "sol-sparc64" ]; then
        LD_LIBRARY_PATH=$SGE_ROOT/lib/sol-sparc:$LD_LIBRARY_PATH
    else
        # all other cluster hosts are x86,
        # and we need the 32 bit library for perl
        LD_LIBRARY_PATH=$SGE_ROOT/lib/sol-x86:$LD_LIBRARY_PATH
    fi
fi
export LD_LIBRARY_PATH

