module load mamba
conda create -n metabp_env
source activate metabp_env
chmod +x bbmap/bbmap/*.sh
chmod +x plass/bin/plass
chmod +x mmseqs/bin/mmseqs
chmod +x clustalo/clustalo-1.2.4-Ubuntu-x86_64
mamba install python=3.9 
mamba install pandas
mamba install -c bioconda eggnog-mapper
conda deactivate
