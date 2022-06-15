module load mamba
conda create -n metabp_env
source activate metabp_env
mamba install python=3.9 
mamba install pandas
mamba install -c bioconda eggnog-mapper
conda deactivate
