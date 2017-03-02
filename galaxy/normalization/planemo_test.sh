planemo conda_init
planemo conda_install .
planemo test --install_galaxy --conda_dependency_resolution

#All 1 test(s) executed passed.
#NmrNormalization[0]: passed

planemo shed_test --install_galaxy -t testtoolshed

#All 1 test(s) executed passed.
#testtoolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_normalization/NmrNormalization/1.0.1[0]: passed

