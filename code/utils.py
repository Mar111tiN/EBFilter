import os
import sys
import re

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


def make_region_list(anno_path, sep):
    # make bed file for mpileup
    # with a pandas dataframe
    # better to open the original file as pandas in this function 
    # mut_pd = pd.read_csv(annopath)
    # mut_pd[1] = mut_pd[1] - 1 - (mut_pd[3] == '-')
    # mut_pd[2] = mut_pd[2] - (mut_pd[3] == '-')
    # mut_pd[selected columns].to_csv(f"{anno_path}.region_list.bed", sep='\t')
    out_path = f"{anno_path}.region_list.bed"
    with open(anno_path) as file_in:
        with open(out_path, 'w') as file_out:
            for line in file_in:
                field = line.rstrip('\n').split(sep)
                # fix for files with header
                if field[0] == 'Chr':
                    continue
                loc = int(field[1]) - (field[4] == "-")  # -1 if field 4 == '-' eg. deletion 
                print(field[0], (loc - 1), loc, file=file_out, sep='\t')
