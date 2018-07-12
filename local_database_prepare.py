"""
Author: ArchieYoung <yangqi2@grandomics.com>
Time: Thu Jul  5 09:24:07 CST 2018
"""
import sys
import argparse
import os
from multiprocessing import Pool
from glob import iglob

from sv_vcf import SV

def vcf_to_db_bed(_args):
    
    vcf, min_support_reads, out_dir, sv_id_prefix = _args
    with open(vcf, "r") as io:
        lines = io.readlines()

        # chromosomes to be kept
        main_chr = (list(range(1, 23)) + ["X", "Y", "MT", "M", "chrX", "chrY",
            "chrM", "chrMT"] + ["chr" + str(i) for i in range(1, 23)])

        # output bedlines
        bed_lines = []

        # check if id is unique
        id_dict = {}
        #chrom1,pos1,chrom2,pos2
        previous_sv_breakpoint = ["NA", "NA", "NA", "NA"]
        for line in lines:
            #skip comment lines
            if line.strip()[0] == "#":
                continue
            sv = SV(line)

            # filter
            if sv.chrom1 not in main_chr or sv.chrom2 not in main_chr:
                continue

            if int(sv.re) < min_support_reads:
                continue

            # remove 'chr' in chromosome id
            sv.chrom1 = sv.chrom1.replace("chr", "")
            sv.chrom2 = sv.chrom2.replace("chr", "")

            # rename sv id
            if sv_id_prefix:
                sv.id = "_".join([sv_id_prefix, sv.id])

            if sv.id not in id_dict:
                id_dict[sv.id] = 1
            else:
                raise RuntimeError("Duplicated SV ID in you VCF "
                "file {}".format(sv.id))
            
            sv_breakpoint = [sv.chrom1, sv.pos1, sv.chrom2, sv.pos2]
            # remove duplicate adjacency BND record in picky vcf
            # Exactly the same
            if ((sv.svtype == "TRA" or sv.svtype == "INV") and
                    sv_breakpoint[:4] == previous_sv_breakpoint[:4]):
                continue
            # just swap breakpoint1 and breakpoint2, still the same
            if ((sv.svtype == "TRA" or sv.svtype == "INV") and
                    sv_breakpoint[:2] == previous_sv_breakpoint[2:] and
                    sv_breakpoint[2:] == previous_sv_breakpoint[:2]):
                previous_sv_breakpoint = sv_breakpoint
                continue

            previous_sv_breakpoint = sv_breakpoint
            # convert to bed format
            # chrom,start,end,svtype,id,svlen,re,info
            if sv.chrom1 == sv.chrom2:
                if int(sv.pos1) > 1:
                    sv.pos1 = str(int(sv.pos1)-1)
                bed_line = "\t".join([sv.chrom1,sv.pos1,sv.pos2,sv.svtype,
                    sv.id,sv.svlen,sv.re])+"\n"
            else: #TRA
                if int(sv.pos1) > 1:
                    pos1_1 = str(int(sv.pos1)-1)
                    pos1_2 = sv.pos1
                elif int(sv.pos1) == 1:
                    pos1_1 = "1"
                    pos1_2 = "1"
                else:
                    continue # invalid position
                if int(sv.pos2) > 1:
                    pos2_1 = str(int(sv.pos2)-1)
                    pos2_2 = sv.pos2
                elif int(sv.pos2) == 1:
                    pos2_1 = "1"
                    pos2_2 = "1"
                else:
                    continue # invalid position
                bed_line1 = "\t".join([sv.chrom1,pos1_1,pos1_2,sv.svtype,
                    sv.id+"_1",sv.svlen,sv.re])+"\n"
                bed_line2 = "\t".join([sv.chrom2,pos2_1,pos2_2,sv.svtype,
                    sv.id+"_2",sv.svlen,sv.re])+"\n"
                bed_line = bed_line1+bed_line2
            bed_lines.append(bed_line)
        
        out_bed_path = os.path.join(out_dir,
                "{}.sv.database.bed".format(sv_id_prefix))
        with open(out_bed_path, "w") as out_hd:
            out_hd.writelines(bed_lines)


def get_args():
    parser = argparse.ArgumentParser(
        description="Prepare Local SV Database",
        usage="usage: %(prog)s [options]")
    parser.add_argument("--vcf_dir",
        help="vcf file directory [default %(default)s]", metavar="STR")
    parser.add_argument("--db_dir",
        help="database directory [default %(default)s]", metavar="STR")
    parser.add_argument("--min_re",
        help="minimum support reads number [default %(default)s]", type=float,
        default=2, metavar="INT")
    parser.add_argument("--threads",
        help="number of threads [default %(default)s]", type=int,
        default=4, metavar="INT")
    
    if len(sys.argv) <= 1:
        parser.print_help()
        exit()
    return parser.parse_args()


def main():
    args = get_args()
    # vcf file list
    vcfs = iglob(os.path.join(args.vcf_dir, "*.vcf"))

    if not os.path.exists(args.db_dir):
        os.mkdir(args.db_dir)

    work_args_list = [(i, args.min_re, args.db_dir,
        os.path.basename(i)[:-4]) for i in vcfs]
    with Pool(processes=args.threads) as pool:
        pool.map(vcf_to_db_bed, work_args_list)

if __name__ == "__main__":
    main()

