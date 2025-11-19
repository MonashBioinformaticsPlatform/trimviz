#!/usr/bin/env bash
tv_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
tv_dir=$(realpath $tv_dir)
echo "trimViz repo is located in: $tv_dir"
echo "creating conda environment..."
if command -v mamba &> /dev/null; then
  mamba env create -f "${tv_dir}/trimviz_conda_minimal.yml"
else
  conda env create -f "${tv_dir}/trimviz_conda_minimal.yml"
fi

# Initialize conda in this shell session (need renv in R - /setup_r_env.R locates this via Sys.getenv("CONDA_PREFIX") )
eval "$(conda shell.bash hook)"
conda activate trimViz2025_renv
echo "setting up R packages via renv..."
#Rscript "${tv_dir}/setup_r_env.R"  "build_from_scratch" "${tv_dir}"
Rscript "${tv_dir}/setup_r_env.R"  "build" "${tv_dir}"
sleep 1
echo "testing to see if we can use R libs..."
Rscript "${tv_dir}/setup_r_env.R" "test" "${tv_dir}"
echo "Done. Remember to use: \"conda activate trimViz2025_renv\" before calling trimviz."
