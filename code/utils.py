import os
import sys
import re
import pandas as pd

def validate_region(region):
    '''
    returns True if region 
    '''
    region_simple = re.compile(r"^[^ \t\n\r\f\v,]+:\d+\-\d+")
    # region format check
    if region:
        region_match = region_exp.match(region)
        if region_match:
            return True

def validate(mut_file, tumor_bam, pon_list):
    # file existence check
    if not os.path.exists(mut_file):
        sys.stderr.write(f"No target mutation file: {mut_file}")
        sys.exit(1)

    if not os.path.exists(tumor_bam):
        sys.stderr.write(f"No target bam file: {tumor_bam}")
        sys.exit(1)

    if not os.path.exists(tumor_bam + ".bai") and not os.path.exists(os.path.splitext(tumor_bam)[0] + '.bai'):
        sys.stderr.write(f"No index for target bam file: {tumor_bam}")
        sys.exit(1)

    if not os.path.exists(pon_list):
        sys.stderr.write(f"No control list file: {pon_list}")
        sys.exit(1)

    with open(pon_list) as hIN:
        for file in hIN:
            file = file.rstrip()
            if not os.path.exists(file):
                sys.stderr.write(f"No control bam file: {file}")
                sys.exit(1)
            if not os.path.exists(file + ".bai") and not os.path.exists(os.path.splitext(file)[0] + '.bai'):
                sys.stderr.write(f"No index for control bam file: {file}")
                sys.exit(1)


def make_region_list(mut_df, out_path, threads):
    # make bed file for mpileup
    # with a pandas dataframe
    # better to open the original file as pandas in this function
    region_pd = mut_df.iloc[:,:5].copy()
    region_pd.iloc[:,1] = mut_df.iloc[:,1] - 1 - (mut_df.iloc[:,4] == '-')
    region_pd.iloc[:,2] = mut_df.iloc[:,1] - (mut_df.iloc[:,4] == '-')
    # outpath: AML033-D.csv --> AML033-D_0.region_list.bed
    outpath = os.path.splitext(out_path)[0]
    if threads > 1: # if thread not -1
        outpath += f"_{os.getpid()}"
    outpath += ".region_list.bed"
    region_pd.iloc[:,:3].to_csv(outpath, sep='\t', header=None, index=False)
    return outpath    
