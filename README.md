# metaBP

## Set-Up

In order to use the tool, users must first clone the code and run the commands in the `install.sh` file provided in the repository. This will allow users access to the mutation pipeline. It is up to the user to ensure anaconda, mamba and java are available. The user must also make sure that the permissions for the build files for the 3rd party tools provided (RBiotools, bbmap, clustalo, mmseqs and plass) are set to "executable."

Once the environment is set up with the necessary packages and the file permissions have been set, the tool is ready for use. Users must ensure the metabp_env environment is active before any commands are run, as well as having java loaded and ready for use.

### RBiotools Set-Up:
 
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
 


