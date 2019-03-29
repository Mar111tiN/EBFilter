#! /usr/bin/env python

import process_vcf
import process_anno
import get_eb_score
import sys, os, subprocess, math, re, multiprocessing 
import vcf, pysam, numpy

region_exp = re.compile('^([^ \t\n\r\f\v,]+):(\d+)\-(\d+)')

def EBFilter_worker_vcf(targetMutationFile, targetBamPath, controlBamPathList, outputPath, mapping_qual_thres, base_qual_thres, filter_flags, is_loption, region, debug_mode):

    controlFileNum = sum(1 for line in open(controlBamPathList, 'r'))

    ##########
    # generate pileup files
    process_vcf.vcf2pileup(targetMutationFile, outputPath + '.target.pileup', targetBamPath, mapping_qual_thres, base_qual_thres, filter_flags, False, is_loption, region)
    process_vcf.vcf2pileup(targetMutationFile, outputPath + '.control.pileup', controlBamPathList, mapping_qual_thres, base_qual_thres, filter_flags, True, is_loption, region)
    ##########

    ##########
    # load pileup files
    pos2pileup_target = {}
    pos2pileup_control = {}

    hIN = open(outputPath + '.target.pileup')
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        pos2pileup_target[F[0] + '\t' + F[1]] = '\t'.join(F[3:])
    hIN.close()

    hIN = open(outputPath + '.control.pileup')
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        pos2pileup_control[F[0] + '\t' + F[1]] = '\t'.join(F[3:])
    hIN.close()
    ##########

    ##########
    # get restricted region if not None
    if is_loption == True and region != "":
        region_match = region_exp.match(region)
        reg_chr = region_match.group(1)
        reg_start = int(region_match.group(2))
        reg_end = int(region_match.group(3))
    ##########

    vcf_reader = vcf.Reader(open(targetMutationFile, 'r'))
    vcf_reader.infos['EB'] = vcf.parser._Info('EB', 1, 'Float', "EBCall Score", "EBCall", "ver0.2.0")
    vcf_writer =vcf.Writer(open(outputPath, 'w'), vcf_reader)


    for vcf_record in vcf_reader:
        current_pos = str(vcf_record.CHROM) + '\t' + str(vcf_record.POS) 

        if is_loption == True and region != "":
            if reg_chr != vcf_record.CHROM: continue
            if int(vcf_record.POS) < reg_start or int(vcf_record.POS) > reg_end: continue

        F_target = pos2pileup_target[current_pos].split('\t') if current_pos in pos2pileup_target else []
        F_control = pos2pileup_control[current_pos].split('\t') if current_pos in pos2pileup_control else []

        current_ref = str(vcf_record.REF)
        current_alt = str(vcf_record.ALT[0])
        var = ""
        if len(current_ref) == 1 and len(current_alt) == 1:
            var = current_alt
        else:
            if len(current_ref) == 1:
                var = "+" + current_alt[1:]
            elif len(current_alt) == 1:
                var = "-" + current_ref[1:]

        EB_score = "." # if the variant is complex, we ignore that
        if not var == "":
            EB_score = get_eb_score.get_eb_score(var, F_target, F_control, base_qual_thres, controlFileNum)

        # add the score and write the vcf record
        vcf_record.INFO['EB'] = EB_score
        vcf_writer.write_record(vcf_record)

    vcf_writer.close()


    # delete intermediate files
    if debug_mode == False:
        subprocess.check_call(["rm", outputPath + '.target.pileup'])
        subprocess.check_call(["rm", outputPath + '.control.pileup'])


def EBFilter_worker_anno(targetMutationFile, targetBamPath, controlBamPathList, outputPath, mapping_qual_thres, base_qual_thres, filter_flags, is_loption, region, debug_mode):

    controlFileNum = sum(1 for line in open(controlBamPathList, 'r'))

    ##########
    # generate pileup files
    process_anno.anno2pileup(targetMutationFile, outputPath + '.target.pileup', targetBamPath, mapping_qual_thres, base_qual_thres, filter_flags, False, is_loption, region)
    process_anno.anno2pileup(targetMutationFile, outputPath + '.control.pileup', controlBamPathList, mapping_qual_thres, base_qual_thres, filter_flags, True, is_loption, region)
    ##########

    ##########
    # load pileup files
    pos2pileup_target = {}
    pos2pileup_control = {}

    hIN = open(outputPath + '.target.pileup')
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        pos2pileup_target[F[0] + '\t' + F[1]] = '\t'.join(F[3:])
    hIN.close()

    hIN = open(outputPath + '.control.pileup')
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        pos2pileup_control[F[0] + '\t' + F[1]] = '\t'.join(F[3:])
    hIN.close()
    ##########

    ##########
    # get restricted region if not None
    if is_loption == True and region != "":
        region_match = region_exp.match(region)
        reg_chr = region_match.group(1)
        reg_start = int(region_match.group(2))
        reg_end = int(region_match.group(3))
    ##########

    hIN = open(targetMutationFile, 'r')
    hOUT = open(outputPath, 'w')

    for line in hIN:

        field = line.rstrip('\n').split('\t')
        chr, pos, pos2, ref, alt = field[0], field[1], field[2], field[3], field[4]
        if alt == "-": pos = str(int(pos) - 1)

        if is_loption == True and region != "":
            if reg_chr != chr: continue
            if int(pos) < reg_start or int(pos) > reg_end: continue

        field_target = pos2pileup_target[chr + '\t' + pos].split('\t') if chr + '\t' + pos in pos2pileup_target else []
        field_control = pos2pileup_control[chr + '\t' + pos].split('\t') if chr + '\t' + pos in pos2pileup_control else [] 

        var = ""
        if ref != "-" and alt != "-":
            var = alt
        else:
            if ref == "-":
                var = "+" + alt
            elif alt == "-":
                var = "-" + ref

        EB_score = "." # if the variant is complex, we ignore that
        if not var == "":
            EB_score = get_eb_score.get_eb_score(var, field_target, field_control, base_qual_thres, controlFileNum)

        # add the score and write the vcf record
        print >> hOUT, '\t'.join(F + [str(EB_score)])

    hIN.close()
    hOUT.close()


    # delete intermediate files
    if debug_mode == False:
        subprocess.check_call(["rm", outputPath + '.target.pileup'])
        subprocess.check_call(["rm", outputPath + '.control.pileup'])



def main(args):

    # should add validity check for arguments
    targetMutationFile = args.targetMutationFile
    targetBamPath = args.targetBamPath
    controlBamPathList = args.controlBamPathList
    outputPath = args.outputPath

    mapping_qual_thres = args.q
    base_qual_thres = args.Q
    filter_flags = args.ff
    thread_num = args.t
    is_anno = True if args.f == 'anno' else False
    is_loption = args.loption
    region = args.region
    debug_mode = args.debug

    # region format check
    if region != "":
        region_match = region_exp.match(region)
        if region_match is None:
            print >> sys.stderr, "Wrong format for --region ({chr}:{start}-{end}): " + region
            sys.exit(1)

    # file existence check
    if not os.path.exists(targetMutationFile):
        print >> sys.stderr, "No target mutation file: " + targetMutationFile
        sys.exit(1)

    if not os.path.exists(targetBamPath):
        print >> sys.stderr, "No target bam file: " + targetBamPath
        sys.exit(1)

    if not os.path.exists(targetBamPath + ".bai") and not os.path.exists(re.sub(r'bam$', "bai", targetBamPath)):
        print >> sys.stderr, "No index for target bam file: " + targetBamPath
        sys.exit(1)


    if not os.path.exists(controlBamPathList):
        print >> sys.stderr, "No control list file: " + controlBamPathList 
        sys.exit(1)

    with open(controlBamPathList) as hIN:
        for file in hIN:
            file = file.rstrip()
            if not os.path.exists(file):
                print >> sys.stderr, "No control bam file: " + file 
                sys.exit(1)

            if not os.path.exists(file + ".bai") and not os.path.exists(re.sub(r'bam$', "bai", file)):
                print >> sys.stderr, "No index control bam file: " + file 
                sys.exit(1)

    if thread_num == 1:
        # non multi-threading mode
        if is_anno == True:
            EBFilter_worker_anno(targetMutationFile, targetBamPath, controlBamPathList, outputPath, mapping_qual_thres, base_qual_thres, filter_flags, is_loption, region, debug_mode)
        else: 
            EBFilter_worker_vcf(targetMutationFile, targetBamPath, controlBamPathList, outputPath, mapping_qual_thres, base_qual_thres, filter_flags, is_loption, region, debug_mode)
    else:
        # multi-threading mode
        ##########

        if is_anno == True:
            # partition anno files
            process_anno.partition_anno(targetMutationFile, outputPath + ".tmp.input.anno.", thread_num)

            jobs = []
            for i in range(thread_num):
                process = multiprocessing.Process(target = EBFilter_worker_anno, args = \
                    (outputPath + ".tmp.input.anno." + str(i), targetBamPath, controlBamPathList, outputPath + "." + str(i), mapping_qual_thres, base_qual_thres, filter_flags, is_loption, region, debug_mode))
                jobs.append(process)
                process.start()
        
            # wait all the jobs to be done
            for i in range(thread_num):
                jobs[i].join()
        
            # merge the individual results
            process_anno.merge_anno(outputPath + ".", outputPath, thread_num)
        
            # delete intermediate files
            if debug_mode == False:
                for i in range(thread_num):
                    subprocess.check_call(["rm", outputPath + ".tmp.input.anno." + str(i)])
                    subprocess.check_call(["rm", outputPath + "." + str(i)])

        else:
            # partition vcf files
            process_vcf.partition_vcf(targetMutationFile, outputPath + ".tmp.input.vcf.", thread_num)

            jobs = []
            for i in range(thread_num):
                process = multiprocessing.Process(target = EBFilter_worker_vcf, args = \
                    (outputPath + ".tmp.input.vcf." + str(i), targetBamPath, controlBamPathList, outputPath + "." + str(i), mapping_qual_thres, base_qual_thres, filter_flags, is_loption, region, debug_mode))
                jobs.append(process)
                process.start()

            # wait all the jobs to be done
            for i in range(thread_num):
                jobs[i].join()

            # merge the individual results
            process_vcf.merge_vcf(outputPath + ".", outputPath, thread_num)

            # delete intermediate files
            if debug_mode == False:
                for i in range(thread_num):
                    subprocess.check_call(["rm", outputPath + ".tmp.input.vcf." + str(i)])
                    subprocess.check_call(["rm", outputPath + "." + str(i)])



