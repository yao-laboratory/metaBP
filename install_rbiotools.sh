source activate metabp_env

module load mamba

mamba install -c r r

mamba install -c r r-rcpp r-rcurl r-gplots r-gridextra

mamba install -c bioconda bioconductor-msa bioconductor-biostrings r-fmsb r-rentrez r-pheatmap r-seqinr 

mamba install -c conda-forge r-stringr r-grimport r-data.table r-ape r-ggplot2 
 
conda deactivate