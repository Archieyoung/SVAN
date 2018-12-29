#!/usr/bin/env python3
"""
author: archieyoung<yangqi2@grandomics.com>
SVAN: A Struture Variation Annotation Tool
"""
import sys
import os
import subprocess
import operator
import argparse
import logging
from multiprocessing import Pool
from glob import iglob


import sv_vcf
import pubdb_prepare
from range_compare import RangeCompare
from insertion_compare import InsCompare
from traslocation_compare import TraCompare 



def queryPrepare_vcf(tmp, prefix, vcf_in, bed_out):
    """
    Convert VCF to BED
    Remove query chrom "chr" prefix if it exists, because database SV
    chrom do not have "chr" prefix.
    Sort bed file by "chrom" and "start", the same as "sort -k1,1 -k2,2n"
    """
    bed_tmp1 = os.path.join(tmp, prefix+".query.bed.tmp1")
    # vcf to bed
    sv_vcf.vcf_to_bed(vcf_in, bed_tmp1)
    # sort bed
    bed_list = []
    query_list = []
    with open(bed_tmp1, "r") as io:
        for line in io:
            line = line.strip()
            fields = line.split("\t")
            query_list.append(line)
            # Remove query chrom "chr" prefix if it exists
            fields[0] = fields[0].replace("chr", "")
            fields[1] = int(fields[1])
            bed_list.append(fields)
    bed_list.sort(key = operator.itemgetter(0, 1))

    with open(bed_out, "w") as io:
        for i in bed_list:
            i[1] = str(i[1])
            io.write("\t".join(i)+"\n")
    return query_list


def queryPrepare_bed(tmp, prefix, bed_in, bed_out):
    """
    Remove query chrom "chr" prefix if it exists, because database SV
    chrom do not have "chr" prefix.
    Sort bed file by "chrom" and "start", the same as "sort -k1,1 -k2,2n"
    """
    # sort bed
    bed_list = []
    query_list = []
    with open(bed_in, "r") as io:
        for line in io:
            line = line.strip()
            fields = line.split("\t")
            query_list.append(line)
            # Remove query chrom "chr" prefix if it exists
            fields[0] = fields[0].replace("chr","")
            fields[1] = int(fields[1])
            bed_list.append(fields)
    bed_list.sort(key = operator.itemgetter(0, 1))
    with open(bed_out, "w") as io:
        for i in bed_list:
            i[1] = str(i[1])
            io.write("\t".join(i)+"\n")
    return query_list


# db fields "chrom","start","end","svtype","svid","annot1","annot2","annot3"...
# format: "annot1;annot1","annot2,annot2";"annot3;annot3","chrom:start-end,svtype;..."
# if multiple result is founded in database, write each feature in one column,
# and seperate by semicolon
def format_pub_result(result, field_names):
    for i in result:
        variations = []
        # init annotations
        annotations = [[] for p in range(len(result[i][0])-5)]
        if (len(field_names)-1) != len(annotations):
            # print(field_names)
            # print(result[i][0][5:])
            raise RuntimeError("fields number error.")
        for j in result[i]:
            variations.append("{}:{}-{},{}".format(j[0],j[1],j[2],j[3]))
            for k in range(len(annotations)):
                annotations[k].append(j[k+5])
        variations_str = ";".join(variations)
        annotations_str = [";".join(l) for l in annotations]
        
        new_fields = [variations_str] + annotations_str
        result[i] = dict() # change the value inplace

        for m in range(len(field_names)):
            result[i][field_names[m]] = new_fields[m] 
    # return result


def format_local_result(result, local_db_size):
    for i in result:
        variations = []
        for j in result[i]:
            # variations.append("{}:{}-{},{}".format(j[0],j[1],j[2],j[3]))
            variations.append(j[4])
        variations_str = ";".join(variations)
        frequency = len(variations)/local_db_size
        if frequency > 1:
            frequency = 1 # set max frequency to 1
        frequency = "{:4.2f}".format(frequency)
        result[i] = {"GrandSV_Frequency":frequency,
                "GrandSV_Variant":variations_str} 
    # return result


def update_result(updated_result, intersect):
    for i in intersect:
        if i not in updated_result:
            updated_result[i] = intersect[i]
        else:
            updated_result[i].update(intersect[i])


def get_sv_by_type(bed_in, bed_out, svtype):
    result_list = []
    with open(bed_in, "r") as io:
        for line in io:
            fields = line.strip().split("\t")
            if fields[3] in svtype:
                result_list.append(line)
    with open(bed_out, "w") as io:
        io.writelines(result_list)


def local_sv_annot(_args):
    (_bedtools, _bed_query, _bed_query_ins, _bed_query_tra, _db_bed,
            _min_overlap, _max_dist, _tmp_dir, _prefix, _id) = _args

    _intersect = dict()
    
    # range compare
    range_intersect = RangeCompare(_bedtools, _bed_query, _db_bed,
            _min_overlap, _tmp_dir, _prefix, "local{}".format(_id))
    
    # insertion compare
    ins_intersect = InsCompare(_bedtools, _bed_query_ins, _db_bed, _max_dist,
            _tmp_dir, _prefix, "local{}".format(_id))

    # tra compare
    tra_intersect = TraCompare(_bedtools, _bed_query_tra, _db_bed, _max_dist,
            _tmp_dir, _prefix, "local{}".format(_id))

    # update intersect result
    _intersect.update(range_intersect)
    _intersect.update(ins_intersect)
    _intersect.update(tra_intersect)

    return _intersect


def table_maker(query_list, result, result_field_names):
    table = []
    query_field_names = ["Chr", "Start", "End", "SVTYPE", "SVID",
            "SVLEN", "RE", "INFO"]
    header = query_field_names + result_field_names
    for i in query_list:
        fields = i.split("\t")
        record_dict = dict([p for p in zip(query_field_names, fields)])
        svid = record_dict["SVID"]
        if svid in result:
            record_dict.update(result[svid])
            for j in result_field_names:
                if j not in record_dict:
                    record_dict[j] = "NA"
        else:
            for j in result_field_names:
                if j not in record_dict:
                    record_dict[j] = "NA"
        # print(record_dict)
        record = "\t".join([record_dict[k] for k in header]) + "\n"
        table.append(record)
    header_str = "\t".join(header) + "\n"
    return header_str, table


def get_args():
    parser = argparse.ArgumentParser(
        description="SVAN: Structure variation annotation",
        usage="usage: %(prog)s [options]")
    parser.add_argument("--vcf",
        help="vcf file [default %(default)s]", metavar="FILE")
    parser.add_argument("--bed",
        help="bed file [default %(default)s]", metavar="FILE")
    parser.add_argument("--pub_dbdir",
        help="public SV database directory [default %(default)s]", metavar="DIR")
    parser.add_argument("--local_dbdir",
        help="local SV database directory [default %(default)s]", metavar="DIR")
    parser.add_argument("--min_overlap",
        help="minimum reciprocal overlap fraction for \"Range\" SV A and"
        " \"Rang\" SV B [default %(default)s]", type=float,
        default=0.5, metavar="FLOAT")
    parser.add_argument("--max_dist",
        help="maximum distance for \"Point\" SV A and"
        " \"Point\" SV B [default %(default)s]", type=int,
        default=1000, metavar="FLOAT")
    parser.add_argument("--threads", default=4,
            help="number of threads [default %(default)s]", type=int, metavar="STR")
    parser.add_argument("--prefix",
        help="output file name [default %(default)s]", metavar="STR")
    parser.add_argument("--outdir",
        help="output directory [default %(default)s]", metavar="DIR")
    parser.add_argument("--tmp",
        help="temporary directory [default %(default)s]",
        default="tmp", metavar="DIR")

    if len(sys.argv) <= 1:
        parser.print_help()
        exit()
    return parser.parse_args()


def main():
    
    args = get_args()
    
    # make dirs
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)

    if not os.path.exists(args.tmp):
        os.mkdir(args.tmp)

    
    # prepare query SV
    bed_query = os.path.join(args.tmp, args.prefix+".query.bed")
    if args.vcf: # vcf query
        query_list = queryPrepare_vcf(args.tmp, args.prefix,
                args.vcf, bed_query)
    elif args.bed: # bed query
        query_list = queryPrepare_bed(args.tmp, args.prefix, 
                args.bed, bed_query)
    else:
        raise RuntimeError("Must provide either a bed or vcf a file!!!")

    
    # prepare database 
    one_thousand_sv_raw = os.path.join(args.pub_dbdir,
        "ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.1000genome.vcf.gz")
    dgv_raw = os.path.join(args.pub_dbdir,
        "DGV.GS.March2016.50percent.GainLossSep.Final.hg19.gff3")
    dbVar_raw = os.path.join(args.pub_dbdir,
        "nstd37.GRCh37.variant_call.vcf.gz")
    decipher_hi_raw = os.path.join(args.pub_dbdir,
        "decipher_HI_Predictions_Version3.bed.gz")

    for i in [one_thousand_sv_raw,dgv_raw,dbVar_raw,decipher_hi_raw]:
        if not os.path.exists(i):
            raise RuntimeError("Error: {} was not founded!".format(i))
    
    one_thousand_sv_db = os.path.join(args.tmp,
        "ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.1000genome.db.bed")
    dgv_db = os.path.join(args.tmp,
        "DGV.GS.March2016.50percent.GainLossSep.Final.hg19.db.bed")
    dbVar_db = os.path.join(args.tmp,
        "nstd37.GRCh37.variant_call.db.bed")
    decipher_hi_db = os.path.join(args.tmp,
        "decipher_HI_Predictions_Version3.db.bed")

    pubdb_prepare.one_thousand_sv.print_bed(one_thousand_sv_raw, one_thousand_sv_db)
    pubdb_prepare.dgv_gold_cnv.print_bed(dgv_raw, dgv_db)
    pubdb_prepare.dbVar_nstd37_sv.print_bed(dbVar_raw, dbVar_db)
    pubdb_prepare.decipher_HI.print_bed(decipher_hi_raw, decipher_hi_db)

    
    bedtools = "bedtools" # temp use

    # annotation field names
    one_thousand_field_names = ["1000g_SV", "1000g_subtype", "1000g_AF",
            "1000g_EAS_AF", "1000g_EUR_AF","1000g_AFR_AF","1000g_AMR_AF", 
            "1000g_SAS_AF"]
    dgv_field_names = ["dgv_SV", "dgv_AF","dgv_sample_size"]
    dbVar_field_names = ["dbVar_SV", "dbVar_CLNSIG","dbVar_PHENO"]
    decipher_hi_field_names = ["Decipher_region", "Decipher_HI"]

    
    # public database SV annotate, 1000g, dgv, dbVar, decipher_hi
    # range annotate
    # 1000g
    one_thousand_range_intersect = RangeCompare(bedtools, bed_query,
            one_thousand_sv_db, args.min_overlap,
            args.tmp, args.prefix, "1000genome")
    format_pub_result(one_thousand_range_intersect, one_thousand_field_names)

    # dgv
    dgv_range_intersect = RangeCompare(bedtools, bed_query,
            dgv_db, args.min_overlap, args.tmp, args.prefix, "dgv")
    format_pub_result(dgv_range_intersect, dgv_field_names)

    # dbVar
    dbVar_range_intersect = RangeCompare(bedtools, bed_query,
            dbVar_db, args.min_overlap, args.tmp, args.prefix, "dbVar")
    format_pub_result(dbVar_range_intersect, dbVar_field_names)

    # decipher_hi
    decipher_hi_range_intersect = RangeCompare(bedtools, bed_query,
            decipher_hi_db, 0, args.tmp, args.prefix,
            "decipher_hi")
    format_pub_result(decipher_hi_range_intersect, decipher_hi_field_names)

    
    # insertion annotate, 1000g and decipher_hi
    # get insertion
    bed_query_ins = os.path.join(args.tmp, args.prefix + ".query.ins.bed")
    get_sv_by_type(bed_query, bed_query_ins, ["INS"])

    # 1000g
    one_thousand_ins_intersect = InsCompare(bedtools, bed_query_ins,
            one_thousand_sv_db, args.max_dist,
            args.tmp, args.prefix, "1000genome")
    format_pub_result(one_thousand_ins_intersect, one_thousand_field_names)

    # decipher hi
    decipher_hi_ins_intersect = InsCompare(bedtools, bed_query_ins,
            decipher_hi_db, args.max_dist, args.tmp, args.prefix,
            "decipher_hi")
    format_pub_result(decipher_hi_ins_intersect, decipher_hi_field_names)

    # translocation annotate, decipher hi
    # get translocation
    bed_query_tra = os.path.join(args.tmp, args.prefix + ".query.tra.bed")
    get_sv_by_type(bed_query, bed_query_tra, ["TRA"])

    # decipher hi
    decipher_hi_tra_intersect = TraCompare(bedtools, bed_query_tra,
            decipher_hi_db, args.max_dist, args.tmp, args.prefix,
            "decipher_hi")
    format_pub_result(decipher_hi_tra_intersect, decipher_hi_field_names)

    # merge public result
    public_results = dict()
    public_intersect_list = [one_thousand_range_intersect,
            dgv_range_intersect,
            dbVar_range_intersect,
            decipher_hi_range_intersect,
            one_thousand_ins_intersect,
            decipher_hi_ins_intersect,
            decipher_hi_tra_intersect
            ]
    for i in public_intersect_list:
        update_result(public_results, i)
    

    # local SV database annotate
    # local SV beds
    local_db_beds = iglob(os.path.join(args.local_dbdir, "*.sv.database.bed"))


    pool = Pool(processes = args.threads)
    intersect_work_list = []

    
    n = 1 # used for number the temp files and calculate local database size
    for local_db_bed in local_db_beds:
        local_annot_args = ([bedtools, bed_query, bed_query_ins, bed_query_tra,
                local_db_bed, args.min_overlap, args.max_dist, args.tmp, args.prefix, n],)
        intersect_work_list.append(pool.apply_async(local_sv_annot,
            local_annot_args))
        n += 1

    pool.close()
    pool.join()

    # merge intersect results
    local_intersect_list = []
    for i in intersect_work_list:
        intersect = i.get()
        local_intersect_list.append(intersect)

    local_intersect_merged = dict()
    for i in local_intersect_list:
        for k in i:
            if k not in local_intersect_merged:
                local_intersect_merged[k] = i[k]
            else:
                local_intersect_merged[k].extend(i[k])
    format_local_result(local_intersect_merged, n)

    # OMG, got the END...
    annot_result = public_results
    update_result(annot_result, local_intersect_merged)
    # for k in annot_result:
    #    print(k, annot_result[k])

    # make table
    # fields
    result_field_names = (["GrandSV_Frequency", "GrandSV_Variant"] +
            one_thousand_field_names +
            dgv_field_names +
            dbVar_field_names +
            decipher_hi_field_names)

    header_str, table = table_maker(query_list, annot_result, result_field_names)

    prefix_with_dir = os.path.join(args.outdir, args.prefix)
    with open(prefix_with_dir + ".svan.tsv", "w") as io:
        io.write(header_str)
        io.writelines(table)


if __name__ == "__main__":
    main()
