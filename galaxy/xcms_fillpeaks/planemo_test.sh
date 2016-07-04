# Example of planemo command to launch test


# -- Use of installed package environments
# after having installing package on a local galaxy instance
source /w/galaxy/dev/shed_tools_tool_dependency_dir/R/3.1.2/iuc/package_r_3_1_2/1ca39eb16186/env.sh
source /w/galaxy/dev/shed_tools_tool_dependency_dir/bioconductor-xcms/1.44.0/lecorguille/package_bioconductor_xcms_1_44_0/0c38f7d43e08/env.sh
planemo test --install_galaxy

#All 1 test(s) executed passed.
#abims_xcms_fillPeaks[0]: passed


# -- Use of conda dependencies
planemo conda_init --conda_prefix /tmp/mc
planemo conda_install --conda_prefix /tmp/mc . 
planemo test --install_galaxy --conda_prefix /tmp/mc --conda_dependency_resolution

#All 1 test(s) executed passed.
#abims_xcms_fillPeaks[0]: passed


# -- Use of shed_test
planemo shed_test --install_galaxy --galaxy_branch "dev" -t testtoolshed
#All 1 test(s) executed passed.
#testtoolshed.g2.bx.psu.edu/repos/lecorguille/xcms_fillpeaks/abims_xcms_fillPeaks/2.0.6[0]: passed
