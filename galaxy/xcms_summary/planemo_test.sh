planemo conda_init
planemo conda_install .
#Linking packages ...
#Error: ERROR: placeholder '/root/miniconda3/envs/_build_placehold_placehold_placehold_placehold_placehold_p' too short in: ncurses-5.9-4



source /w/galaxy/dev/shed_tools_tool_dependency_dir/R/3.1.2/iuc/package_r_3_1_2/1ca39eb16186/env.sh
source /w/galaxy/dev/shed_tools_tool_dependency_dir/bioconductor-camera/1.22.0/lecorguille/package_bioconductor_camera_1_22_0/22cec61d66c2/env.sh
planemo test --install_galaxy

#All 1 test(s) executed passed.
#abims_xcms_summary[0]: passed


planemo shed_test --install_galaxy --galaxy_branch "dev" 

#All 2 test(s) executed passed.
#testtoolshed.g2.bx.psu.edu/repos/lecorguille/xcms_xcmsset/abims_xcms_xcmsSet/2.0.8[0]: passed
#testtoolshed.g2.bx.psu.edu/repos/lecorguille/xcms_xcmsset/abims_xcms_xcmsSet/2.0.8[1]: passed

