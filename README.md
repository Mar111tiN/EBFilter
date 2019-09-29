# EBFilter3.7
empirical bayesian mutation filtering adopted from Genomon-Project EBFilter.

This fork features:
* adaptation to Python 3.7
* fully vectorized computation relying on pandas dataframes
* in-memory computation almost entirely forgoing intermediate phyysical files

## installation
+ dependencies are found in env/eb-env.yml
+ optimally used in conda environment:
  * `$ conda env create -f env/eb-env.yml -m <your_env_name>`
  * `$ conda activate <your_env_name>`
+ install in environment
  * `$ git clone https://github.com/friend1ws/EBFilter.git`
  * `$ cd EBFilter`
  * `$ python setup.py build`
  * `$ python setup.py install`
  
## usage
+ within python code:
  * `from ebfilter import run`
  * required arguments:
    * targetMutationFile: the .vcf or .anno containing the mutations
    * targetBamPath: path to the tumor bam file (+.bai)
    * controlBamPathList: text list of path to PoN bam files (+ .bai)
    * outputPath


## testing
+ from CLI:
  * run `$ EBscore -t 3 testdata/input.anno testdata/tumor.bam testdata/PoN_list.txt output/testdata_EB.csv`

### cache mode
+ from CLI:
  * generate the cache file: `$ makeEBcache -t 3 -bed_file testdata/input.bed testdata/list_normal_sample.txt output/testdata_cache`
  * alternatively (but not recommended) generate the cache file for the entire region covered by the bam files in the PoN list:
  * `$ makeEBcache -t 3 -force_caching testdata/PoN_list.txt output/testdata_cache`

  * run EBscore using th
