#!/usr/bin/env bash
tv_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "trimViz repo is located in: $tv_dir"
echo "creating conda environment..."
conda env create -f "${tv_dir}/trimviz_conda.yml"
conda activate trimViz2025

echo "setting up R packages via renv..."
mkdir "${tv_dir}/Rlibs"
Rscript "${tv_dir}/setup_r_env.R" "build" "${tv_dir}" "${tv_dir}/Rlibs"
sleep 1
echo "testing to see if can use R libs..."
Rscript "${tv_dir}/setup_r_env.R" "test" "${tv_dir}" "${tv_dir}/Rlibs"
