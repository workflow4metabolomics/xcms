#!/bin/bash

TARGET=$1
#TARGET="toolshed"
#TARGET="testtoolshed"
#TARGET="local"


realpath() {
    [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
}

ROOTPATH=$(dirname $(realpath "$0"))

cd $ROOTPATH/msnbase_readmsdata
planemo shed_update -t $TARGET

cd $ROOTPATH/xcms_xcmsset
planemo shed_update -t $TARGET

cd $ROOTPATH/xcms_plot_chromatogram
planemo shed_update -t $TARGET

cd $ROOTPATH/xcms_merge
planemo shed_update -t $TARGET

cd $ROOTPATH/xcms_group
planemo shed_update -t $TARGET

cd $ROOTPATH/xcms_retcor
planemo shed_update -t $TARGET

cd $ROOTPATH/xcms_fillpeaks
planemo shed_update -t $TARGET

cd $ROOTPATH/xcms_summary
planemo shed_update -t $TARGET

cd $ROOTPATH/xcms_export_samplemetadata
planemo shed_update -t $TARGET

cd $ROOTPATH/xcms_repository_suite
planemo shed_update -t $TARGET
