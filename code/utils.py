import os
import sys
import re
import pandas as pd
from functools import partial

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

sign_re = re.compile(r'\^.|\$')
indel_simple = re.compile(r'[\+\-]([0-9]+)')

def clean_up_df(mut_df, pon_count):
    '''
    removes indels and start/end signs from pileup data
    '''
    # globally remove start/end signs
    for i in range(pon_count+1):
        read = f"read{i}"
        mut_df[read] = mut_df[read].str.replace(sign_re, '')

    is_indel = (mut_df['Ref'] == '-') | (mut_df['Alt'] == '-')


    def remove_indels(read, length):
        '''
        removes indel traces in read using indel length argument
        '''

        # construct the regexp string using the indel length
        indel_re = re.compile(r"([ACGTNacgtn])([\+\-])([0-9]+)([ACGTNacgtn]{" + str(length) + "})")
        # do the indel_re substitution on target and control reads
        return indel_re.sub(r'\1', read)

    # function to remove indels
    def clean_indels(i, row):
        # search for the indel length in target read
        indel_length = indel_simple.search(row['read0']).group(1)
        # in every column, remove the indel traces
        for i in range(i+1):
            row[f"read{i}"] = remove_indels(row[f"read{i}"], indel_length)
        return row

    # create partial function for use in df.apply
    # clean_indels_i = partial(clean_indels, pon_count)
    # apply partial clean_indels to remove indel traces in pileup
    mut_df[is_indel] = mut_df[is_indel].apply(partial(clean_indels, pon_count), axis=1)
    # in case there are indels only in the control files
    def clean_this_row(i, row):
        read = row[f"read{i}"]
        print('index: ', i, 'read: ', read)
        m = indel_simple.search(read)
        # check has to be done because apply is applied one extra time (internally checking for optimization?)
        if m:
            length = m.group(1)
            row[f"read{i}"] = remove_indels(row[f"read{i}"], length)
        return row

    for i in range(pon_count):        
        read = f"read{i+1}"
        Q = f"Q{i+1}"
        # boolean mask for non-fitting read-Q-pairs
        not_paired = mut_df[read].str.len() != mut_df[Q].str.len()
        wrong_rows = mut_df[not_paired]
        if len(wrong_rows.index):
            mut_df[not_paired] = wrong_rows.apply(partial(clean_this_row, i+1), axis=1)
    return mut_df


def cleanup_badQ(mut_df, pon_count, filters):
    filter_string = r"([" + filters + "])"
    filter_re = re.compile(filter_string)

    def remove_badQ(i,row):
        Q = row[f"Q{i}"]
        read = row[f"read{i}"]
        while filter_re.search(Q):
            m = filter_re.search(Q)
            Q = filter_re.sub('', Q, count=1)
            pos = max(m.start(),m.end()-1)
            read = read[:pos] + read[pos+1:]
        return row   

    is_snp = (mut_df['Ref'] != '-') & (mut_df['Alt'] != '-')
    for i in range(pon_count):
        Q = f"Q{i}"
        has_badQ = mut_df[Q].str.contains(filter_re) & mut_df[Q].str.len() != 1
        bad_df = mut_df[is_snp & has_badQ]
        mut_df[is_snp & has_badQ] = bad_df.apply(partial(remove_badQ, i), axis=1)

    return mut_df

