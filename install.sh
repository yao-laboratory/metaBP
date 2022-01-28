conda create -n metabp_env

source activate metabp_env

mamba install python=3.9 

mamba install pandas

mamba install -c bioconda/label/cf201901 eggnog-mapper

conda deactivate