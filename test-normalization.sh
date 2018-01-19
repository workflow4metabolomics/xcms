#!/bin/bash

# Variables
script_dir="${0%/*}"
normalization_script=$script_dir/galaxy/normalization/NmrNormalization_wrapper.R
data_dir=$script_dir/galaxy/normalization/test-data

# Run script
Rscript $normalization_script dataMatrix $data_dir/MTBLS1_bucketedData.tabular scalingMethod Total graphType Overlay logOut norm.log dataMatrixOut $script_dir/matrixout.tabular || exit 1

# Test output file
diff $data_dir/MTBLS1_bucketedData_normalized.tabular $script_dir/matrixout.tabular || exit 2
