#! /usr/bin/env python

import collections
import pysam, re, sys

ReIndel = re.compile('([\+\-])([0-9]+)([ACGTNacgtn]+)')
ReStart = re.compile('\^.')
ReEnd = re.compile('\$')

# autovivification
Tree = lambda: collections.defaultdict(Tree)


def varCountCheck(var, depth, baseBar, qualBar, base_qual_thres, verbose):

    # this should be moved to global (just one evaluation should be enough...
    filter_quals = ''
    for qual in range( 33, 33 + base_qual_thres ):
        filter_quals += str( unichr( qual ) )


    if depth == 0: return [0,0,0,0]

    # indel = Tree()
    insertion_p, insertion_n, deletion_p, deletion_n = 0, 0, 0, 0
    insertion_vb_p, insertion_vb_n, deletion_vb_p, deletion_vb_n = 0, 0, 0, 0

    deleted = 0
    iter = ReIndel.finditer(baseBar)
    for m in iter:
        site = m.start()
        type = m.group(1)
        indelSize = m.group(2)
        varChar = m.group(3)[0:int(indelSize)]

        """
        # just leave these codes for the case where I want to evaluate indels in more detail....
        if not (type in indel and varChar.upper() in indel[type]):
            indel[type][varChar.upper()]['+'] = 0
            indel[type][varChar.upper()]['-'] = 0
       
        strand = '+' if varChar.islower() else '-' 
        indel[type][varChar.upper()][strand] += 1
        """

        var_match = False
        if var[0] == "-":
            # for deletion, just checking the size is OK
            if len(var[1:]) == len(varChar):
                var_match = True
        elif var[1:].upper() == varChar.upper():
            var_match = True


        strand = '+' if varChar.isupper() else '-'
        if type == "+":
            if strand == "+":
                insertion_vb_p += 1
                if var_match: insertion_p += 1
            else:
                insertion_vb_n += 1
                if var_match: insertion_n += 1
        else:
            if strand == "+":
                deletion_vb_p += 1
                if var_match: deletion_p += 1
            else:
                deletion_vb_n += 1
                if var_match: deletion_n += 1

        # print type + '\t' + var + '\t' + varChar.upper() + '\t' + strand


        baseBar = baseBar[0:(site - deleted)] + baseBar[(site + int(indelSize) + len(indelSize) + 1 - deleted):]
        deleted += 1 + len(indelSize) + int(indelSize)


    baseBar = ReStart.sub('', baseBar)
    baseBar = ReEnd.sub('', baseBar)

    # error check
    if len(baseBar) != len(qualBar):
        print >> sys.stderr, baseBar + '\n' + qualBar
        print >> sys.stderr, "lengths of bases and qualities are different!"
        sys.exit(1)


    base_num = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "a": 0, "c": 0, "g": 0, "t": 0, "n": 0}

    # for indel check we ignore base qualities
    if var in "ACGTacgt":
        for base, qual in zip(baseBar, qualBar):
            if not ( qual in filter_quals ):
                if base in "ATGCNatgcn":
                    base_num[base] += 1
    else:
        for base in baseBar:
            if base in "ATGCNatgcn":
                base_num[base] += 1


    depth_p = base_num["A"] + base_num["C"] + base_num["G"] + base_num["T"] + base_num["N"]
    depth_n = base_num["a"] + base_num["c"] + base_num["g"] + base_num["t"] + base_num["n"]


    misMatch_p = 0
    misMatch_n = 0
   
    if var in "ACGTacgt":
        misMatch_p = base_num[var.upper()]
        misMatch_n = base_num[var.lower()]    
    else:
        if var[0] == "+":
            if verbose == True:
                misMatch_p, misMatch_n = insertion_vb_p, insertion_vb_n
            else:
                misMatch_p, misMatch_n = insertion_p, insertion_n
        elif var[0] == "-":
            if verbose == True:
                misMatch_p, misMatch_n = deletion_vb_p, deletion_vb_n
            else:
                misMatch_p, misMatch_n = deletion_p, deletion_n
        else:
            print >> sys.stderr, var + ": input var has wrong format!"
            sys.exit(1)


    return [misMatch_p, depth_p, misMatch_n, depth_n]


