# metaBP Mutation Pipeline

## Set-Up

In order to use the tool, users must first clone the code and run the commands in the `install.sh` file provided in the repository. This will allow users access to the mutation pipeline. It is up to the user to ensure anaconda, mamba and java are available. The user must also make sure that the permissions for the build files for the 3rd party tools provided (RBiotools, bbmap, clustalo, mmseqs and plass) are set to "executable."

Once the environment is set up with the necessary packages and the file permissions have been set, the tool is ready for use. Users must ensure the metabp_env environment is active before any commands are run, as well as having java loaded and ready for use.

### RBiotools Set-Up
 
If users want to use the RBiotools version of Linclust (the linear clustering tool used in the pipeline), the commands in the `install_rbiotools.sh` file must also be run. 

Then, the following steps must be followed: 
- Activate the anaconda environment: `conda activate metabp_env`
- Activate the R environment: `R`
- Then use the following command:
```
  install.packages("/path/to/RBiotools_0.5.1.tar.gz", lib="path/to/conda/.conda/envs/new_env/lib/R/library", repos=NULL, type="source")
 ```
 - If installed correctly, this command should work without any errors:
 ```
  library(RBiotools)
 ```
 
## Usage
The tool can be used with two kinds of input, either sequencing read files (.fastq.gz format) or with pre-assembled sequences (.fasta format). 

A sample command using sequencing read files would look something like this:
```
 python annotation_pipeline.py call_mutations -i1 readfile_1.fastq -i2 readfile_2.fastq -o metabp_results -clust 0 -mw 10
```

The parameters for the tool are as follows:
- `-i1` the file path of the first sequence read file (.fastq.gz) to use in paired-end assembly
- `-i2` the file path of the second sequence read file (.fastq.gz) to use in paired-end assembly
- `-s` the file path of the .fasta file of sequences if PLASS assembly is to be bypassed
- `-o` the file path for the output directory
- `-clust` the version of linclust used for clustering. 0 for the original version of Linclust, 1 for the RBiotools version 
- `-mw` the window on either side of a detected mutation within which to check that the similarity threshold (90%) is met. Default = 10

If a file path to a .fasta sequence file is provided using the `-s` flag, PLASS assembly is automatically bypassed.

## Output Files
Output files will be located in the directory provided using the `-o` flag. The output files will include the following:
- assembly.fas: assembled sequences outputed by PLASS
- clusters_all_seqs.fasta, clusters_cluster.tsv, clusters_rep_seq.fasta: Linclust output files. Clusters contain sequences of all lengths.  
- short_prots.fasta: sequences filtered to only include short sequences (length <= 100)
- fragments.fasta: sequences considered to be fragments during the filtering process
- long_proteins.fasta: sequences longer than 100 amino acids in length which were filtered out
- clusters_short_prots.fasta: Clusters filtered to only include short sequences (length <= 100)
- mutations.txt: text file with short protein sequence clusters and mutation information

### Mutation Formatting:
Within each cluster, mutation information is found using the representative sequence as a reference. 

For instance, if a sequence does not match up at location 45, the representative would be annotated with `45A` (the amino acid in the position) and the mutated sequence will be annotated with the differing amino acid eg: `45I`. 

Unmutated sequences (sequences that match the representative sequence) are annotated with the same mutation information as the representative sequence.
