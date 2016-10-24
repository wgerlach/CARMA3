#! /bin/sh

set -x
set -e

aclocal -I .
automake --add-missing
autoconf

