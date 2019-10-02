conda create -n r-ipo --use-local r-ipo r-batch
. ~/miniconda2/bin/activate r-ipo
planemo test --install_galaxy --galaxy_branch "dev" --update_test_data


# -- Use of conda dependencies
planemo conda_init --conda_prefix /tmp/mc
planemo conda_install --conda_prefix /tmp/mc .
planemo test --install_galaxy --conda_prefix /tmp/mc --conda_dependency_resolution

planemo test --conda_auto_init --conda_auto_install --conda_dependency_resolution --conda_prefix /tmp/mc --galaxy_branch dev  > planemo.out 2>&1
