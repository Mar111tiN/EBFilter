# EBFilter3.7
empirical bayesian mutation filtering adopted from Genomon-Project EBFilter.

This fork features:
* adaptation to Python 3.7
* fully vectorized computation relying on pandas dataframes
* in-memory computation almost entirely forgoing intermediate physical files

## installation
+ dependencies are found in env/eb-env.yml
+ optimally used in conda environment:
  * `$ conda env create -f env/eb-env.yml -m <your_env_name>`
+ alternatively, install python=3.7 samtools pandas numpy scipy vcf
+ for jupyter notebook, install jupyter 
  
## usage
+ EBfilter is run from within python as follows:
  * place code folder in your working folder
  * for top-level placement import via `from code import run
  * required arguments:
    * targetMutationFile: the .vcf or .anno containing the mutations
    * targetBamPath: path to the tumor bam file (+.bai)
    * controlBamPathList: text list of path to PoN bam files (+ .bai)
    * outputPath
    
## testing
+ run a jupyter notebook (requires jupyter installation with ipython3 kernel, comes with `$ pip install jupyter`)
+ start the testData.ipynb and run the steps
