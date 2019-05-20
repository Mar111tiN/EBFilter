#! /usr/bin/env python
import sys
import re
import os
import pandas as pd
import numpy as np
import math
import subprocess
import multiprocessing
from multiprocessing import Pool
from functools import partial

from . import utils
from .utils import show_output
from . import anno
from .cache import generate_cache
from .ebscore import generate_EBscore


def main(args, make_cache=False):
    '''
    validates files and refers to respective functions
    '''

    # ############## ARGUMENTS and CONFIG ##################################
    config = utils.get_config(args)
    config['generate_cache'] = make_cache

    # ######################################################################
    # #################### --> GENERATE CACHE ##############################
    if config['generate_cache']:

        # set the makeEB specific config from args and set pon_dict
        config, pon_dict = utils.get_makeEB_config(args, config)

        show_output(generate_cache(pon_dict, config), time=True, color='success')
        return

    # #########################################################################
    # #################### EBscore ############################################
    else:
        # get arguments for EBscore
        config, pon_dict, mut_file, is_anno, tumor_bam, region = utils.get_EBscore_args(args, config)
        show_output(f"Starting EBscore for file {os.path.basename(mut_file)}!", time=True)

        # #################### EBscore ANNO ###################
        if is_anno:
            success = generate_EBscore(mut_file, tumor_bam, pon_dict, region, config)
            show_output(success, time=True, color='success')
            return

        else:   # is_anno is False
            show_output('EBFilter is not YET compatible with .vcf files. ', color='warning')
            return
