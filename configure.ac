AC_INIT(src/main.cpp)

PACKAGE=carma
VERSION=3.0b

####CFLAGS = -ggdb

AM_INIT_AUTOMAKE($PACKAGE, $VERSION)

AC_PROG_CXX
#AC_PROG_LIBTOOL

AX_BOOST_BASE([1.44.0])
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_IOSTREAMS
AX_BOOST_FILESYSTEM
AX_BOOST_SYSTEM


##########################################################################
# debug compilation support
##########################################################################

AC_MSG_CHECKING([whether to build with debug information])
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug],
        [enable debug data generation (def=no)])],
    [debugit="$enableval"],
    [debugit=no])
AC_MSG_RESULT([$debugit])

if test x"$debugit" = x"yes"; then
    AC_DEFINE([DEBUG],[],[Debug Mode])
    AM_CXXFLAGS="$AM_CXXFLAGS -g -Wall  -Wno-uninitialized -O0" # -Werror
else
    AC_DEFINE([NDEBUG],[],[No-debug Mode])
    AM_CXXFLAGS="$AM_CXXFLAGS -O3"
fi

AC_SUBST([AM_CXXFLAGS])
##########################################################################

#AC_SUBST([AM_LDFLAGS])
#AC_SUBST([AM_CPPFLAGS])

#PKG_CHECK_MODULES([LIBCONFIGXX], [libconfig++ >= 1.4],,
#	AC_MSG_ERROR([libconfig++ 1.4 or newer not found.])
#)

###CXXFLAGS="$CXXFLAGS -std=gnu++0x"



AC_OUTPUT(Makefile src/Makefile)
