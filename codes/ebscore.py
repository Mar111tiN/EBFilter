from .eb import get_count_df_snp
from .beta_binomial import fit_bb, bb_pvalues, fisher_combination
from .anno import anno2pileup
from . import utils
import pandas as pd
import numpy as np
import math
from functools import partial


def get_EBscore_from_AB(pen, row):
    '''
    get the EBscore from cached AB parameters
    no fitting is needed as parameters are precomputed and stored in row[5:9]
    '''

    # if row['Chr'] in [str(chrom) for chrom in [3,7,10,15,18,'X']]:
    #     print(row.name, row['Chr'], row['Start'])
    # we only get the snp count_df, using the mut_df 'Alt' as var and adjust for AB_df with column 9
    count_df = get_count_df_snp(row, row['Alt'].upper(), 10)
    bb_params = {}
    if row['depth0'] == np.nan:
        print(row.name, 'Exitinb')
        return None

    # feed-in the AB params coming with the row
    bb_params['p'] = [row['a+'], row['b+']]
    bb_params['n'] = [row['a-'], row['b-']]

    # get the dataSeries for read0
    target_df = count_df.loc['read0']
    p_values = bb_pvalues(bb_params, target_df)

    # ########### FISHER COMBINATION #########################
    # perform Fisher's combination methods for integrating two p-values of positive and negative strands
    EB_pvalue = fisher_combination(p_values)
    EB_score = 0
    if EB_pvalue < 1e-60:
        EB_score = 60
    elif EB_pvalue > 1.0 - 1e-10:
        EB_score = 0
    else:
        EB_score = -round(math.log10(EB_pvalue), 3)
    return EB_score


def get_EB_from_cache(snp_df, tumor_bam, config, chrom):
    '''
    takes a shortened annotation file (snps only) and loads the valid AB parameters and gets the EB score for these lines
    '''

    # ############################# GET THE AB parameters
    AB_df = utils.get_AB_df(chrom, config)
    # merge the BB params into the snp_df
    snpAB_df = pd.merge(snp_df, AB_df, on=['Chr', 'Start', 'Alt'])
    snpAB_df = anno2pileup(snpAB_df, tumor_bam, None, None, config)
    # remove start/end signs
    utils.cleanup_df(snpAB_df, 0, config)
    snpAB_df.dropna(inplace=True)
    snpAB_df['EB_score'] = snpAB_df.apply(partial(get_EBscore_from_AB, config['fitting_penalty']), axis=1)
    # return raw df with minimal columns
    return snpAB_df.loc[:, ['Chr', 'Start', 'EB_score']]


def generate_EBscore(mut_file, config):
    print(f'Loading annotation file {mut_file}')
    # create anno_df, store original other info and get anno_chr list
    anno_df, original_columns, config['anno_chr'] = utils.read_anno_csv(mut_file, config)
    # create small working copy
    mut_df = anno_df[anno_df.columns[:5]].copy()
    print(f'{len(mut_df.index)} variants detected.')
    # ############## CACHE MODE ###############
    EB_dfs = []
    insert = ''
    indel_count = 0

    if config['cache_mode']:
        # ################ CACHE RETRIEVAL ####################################
        # do multi-threaded merging of the different AB_files per chromosome

        # get the SNPs for the annotation
        snp_df = mut_df.query('not (Ref == "-" or Alt == "-")')
        # get the valid chromosomes
        print(f'{len(snp_df.index)} SNVs detected - pulling BB parameters from cache and piling up the target bam')
        cache_chromes = list(set(config['chr']) & set(config['cache_chr']))
        # init multicore
        getEB_pool = Pool(threads)
        # threads are mapped to the pool of chromosomes
        # outputs df_tuples of (has_EB, no_EB)
        EB_dfs += getEB_pool.map(partial(get_EB_from_cache, snp_df, tumor_bam, config), cache_chromes)
        # EBcached_outputs will be concatenated with the annoworker dfs
        getEB_pool.close()
        getEB_pool.join()

    # ################ NOT CACHE MODE // INDELS IN CACHE MODE ###############
        # reduce mut_df to the remaining indel occurrences and uncached chroms for anno worker
        mut_df_indel = mut_df.query('Ref == "-" or Alt == "-"')
        mut_df_chr_no_cache = mut_df.query(f'Chr not in {config["cache_chr"]}')

        # get lengths for summary
        indel_count = len(mut_df_indel.index)
        chrom_not_in_cache_count = len(mut_df_chr_no_cache.index)
        # make list of remaining mut_dfs for concat below
        remaining_mut_dfs = [mut_df_chr_no_cache, mut_df_indel]

        # get the index of SNVs missing in cache (for whatever reason beside missing chromosome)
        cache_EB_df = pd.concat(EB_dfs)
        missing_in_cache_index = cache_EB_df.query('EB_score != EB_score').index
        missing_in_cache_count = len(missing_in_cache_index)
        if missing_in_cache_count:
            # generate mut_df of missing SNVs from cache (using index of cache_EB_df with EBscore == NaN) and add to df list 
            remaining_mut_dfs += mut_df.iloc(missing_in_cache_index)
            EB_dfs = [cache_EB_df.query('EB_score == EB_score')]

        # gather up the remaining mutations for non_cached computation
        mut_df = pd.concat(remaining_mut_dfs)
        print(mut_df[:53])

        if missing_in_cache_count or chrom_not_in_cache_count:
            insert = f" and {missing_in_cache_count + chrom_not_in_cache_count} SNVs not found in EBcache"
        print("Cache_mode off")

    print(f"Computing EBscores for {indel_count} indels{insert}.")
    config['cache_mode'] = False

    # if all are snps, mut_df is now empty after caching and would raise error in if clause
    mut_len = len(mut_df.index)
    if mut_len:
        #  setup the processor pool
        anno_pool = Pool(threads)
        # set the minimum number of lines for one thread to 2000
        split_factor = min(math.ceil(mut_len / 2000), threads)

        # mut_split is the argument pool for anno_pool
        mut_split = np.array_split(mut_df, split_factor)
        print(f"{dt.now().strftime('%H:%M:%S')} Piling up target and PoN for EBscore computation")
        # run partial function anno_partial with mut_df as remaining argument to iterate over for multiprocessing
        out_dfs = anno_pool.map(partial(anno.worker, tumor_bam, pon_dict, region, config), mut_split)  # mut_split is the iterable df_pool
        anno_pool.close()
        anno_pool.join()
        # add the anno.worker output to EB_dfs
        EB_dfs += out_dfs
    # concat EB_dfs to one
    EB_df = pd.concat(EB_dfs)
    EB_df = EB_df.sort_values([EB_df.columns[0], EB_df.columns[1]])

    # add the original columns back to the annotation file
    final_df = pd.merge(left=anno_df, right=EB_df, how='outer', on=['Chr', 'Start'])
    if original_columns is not None:
        final_df.columns = list(original_columns) + ['EB_score']
    output_path = config['output_path']
    print(f'Writing annotation file {output_path} with EBscores to disc..')
    final_df.to_csv(output_path, sep='\t', index=False)

    # ############### DEBUG #####################################
    if config['debug_mode']:
        out_file = output_path.replace('eb', "eb_only")
        EB_df.to_csv(out_file, sep=config['sep'], index=False)
    # ############# ##############################################
    return f"{dt.now().strftime('%H:%M:%S')} EBscore is finished!"
