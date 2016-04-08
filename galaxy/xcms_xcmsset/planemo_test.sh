planemo conda_init
planemo conda_install .
planemo test --install_galaxy --conda_dependency_resolution --galaxy_branch "dev"

#All 2 test(s) executed passed.
#abims_xcms_xcmsSet[0]: passed
#abims_xcms_xcmsSet[1]: passed




source /w/galaxy/dev/shed_tools_tool_dependency_dir/R/3.1.2/iuc/package_r_3_1_2/1ca39eb16186/env.sh
source /w/galaxy/dev/shed_tools_tool_dependency_dir/bioconductor-xcms/1.44.0/lecorguille/package_bioconductor_xcms_1_44_0/0c38f7d43e08/env.sh
planemo test --install_galaxy --galaxy_branch "dev"

#All 2 test(s) executed passed.
#abims_xcms_xcmsSet[0]: passed
#abims_xcms_xcmsSet[1]: passed