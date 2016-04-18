source /w/galaxy/dev/shed_tools_tool_dependency_dir/R/3.1.2/iuc/package_r_3_1_2/1ca39eb16186/env.sh
source /w/galaxy/dev/shed_tools_tool_dependency_dir/bioconductor-camera/1.22.0/lecorguille/package_bioconductor_camera_1_22_0/22cec61d66c2/env.sh
planemo test --install_galaxy --galaxy_branch "dev"

#All 1 test(s) executed passed.
#abims_CAMERA_annotateDiffreport[0]: passed


planemo shed_test --install_galaxy -t testtoolshed --galaxy_branch "dev"

#All 1 test(s) executed passed.
#testtoolshed.g2.bx.psu.edu/repos/lecorguille/camera_annotate/abims_CAMERA_annotateDiffreport/2.1.3[0]: passed
