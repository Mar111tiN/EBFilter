import os
import sys
import csv
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

def validate(file, message):
    '''
    file existence checks
    '''
    if not os.path.exists(file):
        sys.stderr.write(f"{message}: {file}")
        sys.exit(1)
    else:
        return file


def validate_bam(bam_file):
    '''
    file existence for bam files and accompanying bai files
    '''
    validate(bam_file, "No control bam file")
    if not os.path.exists(bam_file + ".bai") and not os.path.exists(os.path.splitext(bam_file)[0] + '.bai'):
        sys.stderr.write(f"No index for control bam file: {bam_file}")
        sys.exit(1)
    return bam_file


def validate_pon(pon_list):
    '''
    file existence check for pon_list and the containing bam (and bai) files
    '''
    validate(pon_list, "No control list file")
    with open(pon_list) as file_list:
        for file in file_list:
            bam_file = file.rstrip()
            validate_bam(bam_file)
    return pon_list


def read_anno_csv(mut_file, state):
    '''
    reads in the mutation file and resets the relevant header columns required for dataframe operations
    --> returns the dataframe and the original header names
    if no header is detected, the relevant files are names as needed and additional columns are named other1, other2,...
    '''
    def to_int(Chr_name):
        '''
        converts all number chromosomes to int
        '''
        try:
            return int(Chr_name)
        except ValueError:
            return Chr_name       


    with open(mut_file, 'r') as input_file:
        has_header = csv.Sniffer().has_header(input_file.read(1024))

    sep = state['sep']
    with open(state['log'],'w+') as log:
        print(f'Loading annotation file {mut_file} into dataframe', file=log)


        if has_header:
            print(f'Header detected', file=log)
            anno_df = pd.read_csv(mut_file, sep=sep, converters={0:to_int, 1:to_int, 2:to_int})
            org_columns = anno_df.columns
            anno_df.columns = ['Chr','Start','End','Ref', 'Alt'] + list(anno_df.columns[5:])
        else:
            anno_df = pd.read_csv(mut_file, sep=sep, header=None, converters={0:to_int, 1:to_int, 2:to_int})
            org_columns = None
            rest_columns = [f'other{i+1}' for i in range(len(anno_df.columns) - 5)]
            anno_df.columns = ['Chr','Start','End','Ref', 'Alt'] + rest_columns
            print(f'Adding header: Chr\tStart\tEnd\tRef\tAlt\tother1...', file=log)
    return (anno_df.sort_values(['Chr', 'Start']), org_columns)


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


def clean_up_df(mut_df, pon_count):
    '''
    removes indels and start/end signs from pileup data
    '''
    # regexps for start/end sign
    sign_re = re.compile(r'\^.|\$')
    # regexps for indels with a group for getting indel length
    indel_simple = re.compile(r'[\+\-]([0-9]+)')
    # globally remove start/end signs
    for i in range(pon_count+1):
        read = f"read{i}"
        mut_df[read] = mut_df[read].str.replace(sign_re, '')

    is_indel = (mut_df['Ref'] == '-') | (mut_df['Alt'] == '-')


    def remove_indels(read, length):
        '''
        removes indel traces in read using indel length argument
        indels on pos strand are turned into '-'
        indels on neg strand are turned into '_'
        '''
        # construct the regexp string using the indel length
        pos_indel_re = re.compile(r"([ACGTN])([\+\-])([0-9]+)([ACGTNacgtn]{" + str(length) + "})")
        neg_indel_re = re.compile(r"([acgtn])([\+\-])([0-9]+)([ACGTNacgtn]{" + str(length) + "})")
        # do the indel_re substitution on target and control reads
        read = pos_indel_re.sub('-', read)
        return neg_indel_re.sub('_', read)

    # function to remove indels
    def clean_indels(i, row):
        # search for the indel length in target read
        m = indel_simple.search(row['read0'])
        if m:
            indel_length = m.group(1)
        else:
            sys.stderr.write(f"No indel detected in read pileup at position {row['Start']}.\t Please check validity of annotation file!")
            sys.exit(1)
        # in every column, remove the indel traces
        for i in range(i+1):
            row[f"read{i}"] = remove_indels(row[f"read{i}"], indel_length)
        return row

    # apply partial clean_indels to remove indel traces in pileup
    mut_df[is_indel] = mut_df[is_indel].apply(partial(clean_indels, pon_count), axis=1)
    # in case there are indels only in the control files
    def clean_this_row(i, row):
        read = row[f"read{i}"]
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
    '''
    removes base qualities below given threshold ( state['Q'] )..
        and corresponding read bases
    this function should be redundant as the -Q option..
        already applies during mpileup
    '''
    # regexps for the badQ filtering
    filter_string = r"([" + filters + "])"
    filter_re = re.compile(filter_string)

    def remove_badQ(i,row):
        '''
        used as partial callable within the pd.apply
        '''
        Q = row[f"Q{i}"]
        read = row[f"read{i}"]
        while filter_re.search(Q):
            m = filter_re.search(Q)
            Q = filter_re.sub('', Q, count=1)
            pos = max(m.start(),m.end()-1)
            read = read[:pos] + read[pos+1:]
        return row   
    # set boolean mask for snp (we ignore indels)
    is_snp = (mut_df['Ref'] != '-') & (mut_df['Alt'] != '-')
    for i in range(pon_count):
        # go through all Q columns and look for bad qualities
        Q = f"Q{i}"
        has_badQ = mut_df[Q].str.contains(filter_re) & mut_df[Q].str.len() != 1
        bad_df = mut_df[is_snp & has_badQ]
        mut_df[is_snp & has_badQ] = bad_df.apply(partial(remove_badQ, i), axis=1)

    return mut_df

