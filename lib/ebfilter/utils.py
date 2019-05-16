import os
import sys
import csv
import re
import pandas as pd
from functools import partial
import subprocess
from subprocess import Popen, PIPE, DEVNULL
from io import StringIO

# ################# GLOBALS #############################################
# regexps for start/end sign
sign_re = re.compile(r'\^.|\$')
# regexps for indels with a group for getting indel length
indel_simple = re.compile(r'[\+\-]([0-9]+)')
region_simple = re.compile(r"^[^ \t\n\r\f\v,]+:\d+\-\d+")


# ############ MISC ############################################
################################################################

def split_bam(chrom, pon_folder, pon_row):
    '''
    creates a sub bam file (+bai) of the specified chromosomal region per pon_list row
    returns the name of the bam file for creating a sub pon_list
    '''

    bam_file = pon_row[0]
    # create the name for the sub-bam file as org_bam_name_chr?.bam (and .bai)
    bam_out = os.path.join(pon_folder, f"{os.path.splitext(os.path.basename(bam_file))[0]}_{str(chrom)}.bam")
    split_bam_cmd = ["samtools", "view", "-b", "-o", bam_out, bam_file, str(chrom)]
    bam_index_cmd = ["samtools", "index", bam_out]
    subprocess.check_call(split_bam_cmd)
    subprocess.check_call(bam_index_cmd)
    return bam_out


def delete_pom_bams(bam_file):
    '''
    delete a bam and its accompanying bai file
    '''
    subprocess.check_call(['rm', bam_file, f"{bam_file}.bai"])


def sort_chr(dict):
    '''
    sorts all types of chrom lists
    '''

    chr = dict['chr'].replace('Chr', '').replace('chr', '')
    assigner = {'X':50, 'Y':60, 'M':70, '*':80}
    try:
        chr = int(chr)
    except ValueError:
        if chr in ['X', 'Y', 'M', '*']:
            chr = assigner[chr]
        else:
            chr = 100
    return chr


def make_region_list(mut_df, config):
    '''
    make bed file for mpileup from mut_df
    # better to open the original file as pandas in this function
    '''

    region_pd = mut_df.iloc[:, :5].copy()
    region_pd.iloc[:, 1] = mut_df.iloc[:, 1] - 1 - (mut_df.iloc[:, 4] == '-')
    region_pd.iloc[:, 2] = mut_df.iloc[:, 1] - (mut_df.iloc[:, 4] == '-')
    # outpath: AML033-D.csv --> AML033-D_0.region_list.bed
    outpath = os.path.splitext(config['output_path'])[0]
    if config['threads'] > 1:     # if thread not -1
        outpath += f"_{os.getpid()}"
    outpath += ".region_list.bed"
    region_pd.iloc[:, :3].to_csv(outpath, sep='\t', header=None, index=False)
    return outpath

# ###################### VALIDATIONS ###############################
####################################################################


def validate_region(region):
    '''
    returns True if region 
    '''

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


def validate_pon(pon_list, config):
    '''
    file existence check for pon_list and the containing bam (and bai) files
    returns a tuple of a dict containing the pon_list and the pon_df and the chr_list of containing chroms
    '''

    print(f"Validating PoN list {pon_list}..")
    pon_df = pd.read_csv(validate(pon_list, "No PanelOfNormals list file"), header=None)
    pon_df[0].apply(validate_bam)
    config['pon_chr'] = pon2chr_list(pon_df)
    return {'list': pon_list, 'df': pon_df}


def validate_bed(bed_file, config):
    '''
    check for existence of bed_file and if chroms are compatible with PONs
    '''

    bed_file = validate(bed_file, 'No bed file found. Exiting..')
    bed_chr = bed2chr_list(bed_file)
    if not set(bed_chr).issubset(set(config['pon_chr'])):
        sys.stderr.write(f"Bed file {bed_file} contains chroms not found in PanelOfNormals. Exiting..")
        sys.exit(1)
    return (bed_file, bed_chr)


# ####################### I/O ####################################
##################################################################

def read_anno_csv(mut_file, config):
    '''
    reads in the mutation file and resets the relevant header columns required for dataframe operations
    --> returns the dataframe and the original header names
    if no header is detected, the relevant files are names as needed and additional columns are named other1, other2,...
    '''

    def to_int(Chr_name):
        '''
        converts all number chroms to int
        '''
        try:
            return int(Chr_name)
        except ValueError:
            return Chr_name

    def check_columns(mut_file, anno_df, config):
        '''
        raises error, if only one column was detected
        '''

        row, col = anno_df.shape
        if col == 1:
            sys.stderr.write(f"Only one column detected in {mut_file} - I am guessing wrong separator ( {config['sep']} )?")
            sys.exit(1)

    with open(mut_file, 'r') as f:
        has_header = csv.Sniffer().has_header(f.readline())

    sep = config['sep']

    if has_header:
        print(f'Header detected')
        anno_df = pd.read_csv(mut_file, sep=sep, dtype={'Chr':str}, converters={1: to_int, 2: to_int})
        org_columns = anno_df.columns
        check_columns(mut_file, anno_df, config)
        anno_df.columns = ['Chr', 'Start', 'End', 'Ref', 'Alt'] + list(anno_df.columns[5:])

    else:
        anno_df = pd.read_csv(mut_file, sep=sep, header=None, dtype={'Chr':str}, converters={1: to_int, 2: to_int})
        check_columns(mut_file, anno_df, config)
        org_columns = None
        rest_columns = [f'other{i+1}' for i in range(len(anno_df.columns) - 5)]
        anno_df.columns = ['Chr', 'Start', 'End', 'Ref', 'Alt'] + rest_columns
    # retrieve chrom list occurring in anno_file
    anno_chr = [str(chrom) for chrom in anno_df.iloc[:, 0].unique()]
    if set(anno_chr).issubset(set(config['pon_chr'])):
        # active chroms are only the ones found in the anno file and in config['chr']
        config['chr'] = list(set(anno_chr) & set(config['chr']))
        # if chrom is given, reduce the anno file to the provided chroms
        anno_df = anno_df.query(f"Chr in {config['chr']}")
    else:   # raise error if unknown chroms occur in pon file
        sys.stderr.write(f"Found chromosomes in anno file {mut_file} that are not found in PanelOfNormals. Exiting..")
        sys.exit(1)
    # return tuples of anno_df, original columns and chroms occurring in anno_chr
    return (anno_df.sort_values(['Chr', 'Start']), org_columns, anno_chr)


# ################### CHROMOSOME EXTRACTIONS ######################
###################################################################

def bam2chr_list(bam_file):
    '''
    creates a list of chrom names for the input bam
    '''

    bam_stats_cmd = ['samtools', 'idxstats', bam_file]
    bam_stats = Popen(bam_stats_cmd, stdout=PIPE, stderr=DEVNULL)
    bam_stats_string = StringIO(bam_stats.communicate()[0].decode('utf-8'))
    stats_df = pd.read_csv(bam_stats_string, sep='\t', header=None)
    non_empty = stats_df[stats_df[2] > 100]
    return [str(chrom) for chrom in list(non_empty[0].T)]


def pon2chr_list(pon_df):
    '''
    generate a chrom list from the pon_list
    '''

    chr_set = set()
    pon_df[0].apply(lambda bam_file: chr_set.update(bam2chr_list(bam_file)))
    return [str(chrom) for chrom in list(chr_set)]


def bed2chr_list(bed_file):
    '''
    generates a chrom list from a annotated mutation file
    '''

    bed_df = pd.read_csv(bed_file, sep='\t', dtype={0: str}, header=None, skiprows=10)
    # return the list of unique values from the first row (Chr row)
    return [str(chrom) for chrom in bed_df.iloc[:, 0].unique()]


def anno2chr_list(anno_file):
    '''
    generates a chrom list from a annotated mutation file
    '''

    anno_df = pd.read_csv(anno_file, sep='\t', dtype={0: str}, header=None)
    # return the list of unique values from the first row (Chr row)
    return [str(chrom) for chrom in anno_df.iloc[:, 0].unique()]


# ################# CACHING ####################################
################################################################


def validate_cache(cache_folder, config):
    '''
    file existence check for cache folder and the containing cache files
    returns the validated cache_folder
    '''

    # check existence of cache folder
    if not os.path.isdir(cache_folder):
        sys.stderr.write(f"Cache folder {cache_folder} cannot be found! Exiting..")
        sys.exit(1)

    # check existence of cache files for the chroms in the pon folder and store in cache_chr
    cache_file_tuples = [(os.path.join(cache_folder, f"{chrom}.cache"), chrom) for chrom in config['pon_chr']]
    cache_chr = []
    # store existing cache_files in cache_chr
    for cache_file_tuple in cache_file_tuples:
        if os.path.isfile(cache_file_tuple[0]):
            cache_chr.append(cache_file_tuple[1])
    if len(cache_chr) == 0:
        sys.stderr.write(f"Cache folder {cache_folder} contains no usable cache files! Exiting..")
        sys.exit(1)
    return (cache_folder, cache_chr)


def check_cache_files(config):
    '''
    checks all active chromosoms for a preexisting cache file and only returns the list of missing cache files
    '''

    not_cached = []
    for chrom in config['chr']:
        cache_file = os.path.join(config['cache_folder'], f"{chrom}.cache")
        if os.path.isfile(cache_file):
            print(f"Cache file {cache_file} found. Does not need to be generated.")
        else:
            not_cached.append(chrom)
    return not_cached


def check_pileup_files(config):
    '''
    checks, whether pileup files already exist and returns the chrom list for non_existing pileups
    returns list of missing pileup files and list of pileup_dicts for existing pileups
    '''

    not_piled_up = []
    already_piledup_dicts = []
    for chrom in config['chr']:
        pileup_file = os.path.join(config['pileup_folder'], f"cache_{chrom}.pileup")
        if os.path.isfile(pileup_file):
            pile_len = len(pd.read_csv(pileup_file, sep='\t').index)
            print(f"Pileup file {pileup_file} found. Does not need to be built again.")
            already_piledup_dicts.append({'file': pileup_file, 'chr': chrom, 'pileup_len': pile_len})
        else:
            not_piled_up.append(chrom)

    return not_piled_up, already_piledup_dicts


def clean_read_column(read_series):
    '''
    removes per column all reads with start/end signs
    '''

    # should include indel removal as well?
    return read_series.str.replace(sign_re, '')


def cleanup_df(mut_df, pon_count, config):
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

    def clean_this_row(i, row):
        read = row[f"read{i}"]
        m = indel_simple.search(read)
        # check has to be done because apply is applied one extra time (internally checking for optimization?)
        if m:
            length = m.group(1)
            row[f"read{i}"] = remove_indels(row[f"read{i}"], length)
        return row

    # only clean up other columns if pon pileup is used (non cache_mode)
    if not config['cache_mode']:
        # apply partial clean_indels to remove indel traces in pileup
        mut_df[is_indel] = mut_df[is_indel].apply(partial(clean_indels, pon_count), axis=1)
        # in case there are indels only in the control files

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
    removes base qualities below given threshold ( config['Q'] )..
        and corresponding read bases
    this function should be redundant as the -Q option..
        already applies during mpileup
    '''
    # regexps for the badQ filtering
    filter_string = r"([" + filters + "])"
    filter_re = re.compile(filter_string)

    def remove_badQ(i, row):
        '''
        used as partial callable within the pd.apply
        '''
        Q = row[f"Q{i}"]
        read = row[f"read{i}"]
        while filter_re.search(Q):
            m = filter_re.search(Q)
            Q = filter_re.sub('', Q, count=1)
            pos = max(m.start(), m.end()-1)
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


def get_AB_df(chrom, config):
    '''
    load and reformat the AB-cache file for a chromosome
    '''

    cache_file = os.path.join(config['cache_folder'], f"{chrom}.cache")
    # the merger columns
    cols = ['Chr', 'Start', 'Alt']
    # load the file and set Chr and Start to index for proper setting of multi-index
    # print(f"Loading cache {cache_file}..")
    AB_df = pd.read_csv(cache_file, sep=',', compression='gzip', dtype={'Chr': str, 'Start': int}).set_index(cols[:2])
    AB_columns = pd.MultiIndex.from_product([['A', 'C', 'T', 'G'], ['+', '-'], ['a', 'b']], names=['var', 'strand', 'param'])
    # set multi-indexed columns
    AB_df.columns = AB_columns
    # stack the var column level for merge with the snp_df
    AB_df = AB_df.stack('var')
    # reduce the column index level to 1
    AB_df.columns = AB_df.columns.droplevel(0)
    # unset the indices and transfer to columns
    AB_df = AB_df.reset_index()
    # rename columns for merge
    AB_df.columns = cols + ['a+', 'b+', 'a-', 'b-']
    return AB_df


# ################# PROCESSING ##################################
################################################################

def show_command(command_list, config):
    '''
    prints the command line if debugging is active
    '''

    if config['debug_mode']:
        print('\033[1m', '$ ' + ' '.join(command_list), '\033[0m')
    return
