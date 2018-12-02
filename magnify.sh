#!/bin/bash

magnify_source="$(dirname $(readlink -f $BASH_SOURCE))"

rootfile="$1" ; shift
if [[ "$rootfile" =~ :// ]] ; then
    echo "Loading URL $rootfile"
else
    rootfile=$(readlink -f "$rootfile")
fi
sign="${1:-0}"
# startdir=$(pwd)
# rebin="${2:-4}"

# echo "Loading frame \"$frame\" rebin \"$rebin\""

cd $magnify_source/scripts

# echo $rootfile
# echo $frame
# echo $rebin

root -l loadClasses.C Magnify.C'("'"$rootfile"'", '$sign')'
