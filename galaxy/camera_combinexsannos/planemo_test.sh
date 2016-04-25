# Example of planemo command to launch test

# Note: --galaxy_branch "dev" is set to deal with zip file


# -- Use of conda dependencies
planemo conda_init --conda_prefix /tmp/mc
planemo conda_install --conda_prefix /tmp/mc . 
planemo test --install_galaxy --conda_prefix /tmp/mc --conda_dependency_resolution --galaxy_branch "dev"

#All 1 test(s) executed passed.
#abims_CAMERA_combinexsAnnos[0]: passed


# -- Use of shed_test
planemo shed_test --install_galaxy -t testtoolshed

#All 1 test(s) executed passed.
#testtoolshed.g2.bx.psu.edu/repos/mmonsoor/camera_combinexsannos/abims_CAMERA_combinexsAnnos/2.0.4[0]: passed
