"""
Author: ArchieYoung <yangqi2@grandomics.com>
Time: Thu Jul  5 09:24:07 CST 2018
"""
import os
import subprocess
import sys

from bed_intersect import intersect_f
from table_check import TableCheck

def TraCompare(bedtools, bedA, bedB, max_dist, tmp_dir, prefix, db_id):
    # result sv dict, key: query sv id, value: db fields
    tra_compare_result = {}

    # calculate query and database fields number
    query_field_num = TableCheck(bedA)
    # db_field_num = table_check(bedB)

    # padding
    padding = max_dist
    padding_tmp_bed = os.path.join(tmp_dir, 
            "{}.{}.tra.padding{}.tmp.bed".format(prefix, db_id, padding))
    
    with open(padding_tmp_bed, "w") as padding_io:
        with open(bedA, "r") as bedA_io:
            for line in bedA_io:    
                line = line.strip()
                if line[0] == "#":
                    print(line, file = padding_io)
                else:
                    fields = line.split("\t")
                    bed_start = int(fields[1])
                    bed_end = int(fields[2])
                    if bed_start > padding: # make sure start not less than 1
                        bed_start = bed_start - padding
                    else:
                        bed_start = 1
                    bed_end = bed_end + padding

                    fields[1] = str(bed_start)
                    fields[2] = str(bed_end)

                    print("\t".join(fields), file = padding_io)

    intersect_tmp_bed = os.path.join(tmp_dir,
            "{}.{}.tra.intersect.tmp.bed".format(prefix, db_id))

    # bedtools intersect, padding + intersect to get all SVs near the brk
    with open(intersect_tmp_bed,"w") as io:
        subprocess.run([bedtools, "intersect", "-a", padding_tmp_bed, "-b",
            bedB,"-wo"], stdout = io)

    # read intersect file
    # chrom start end svtype svid svlen re info;
    # svid is a unique identifier for the sv
    with open(intersect_tmp_bed, "r") as io:
        tra_pair = {}
        lines = io.readlines()
        for line in lines:
            fields = line.strip().split("\t")
            query_fields = fields[:query_field_num]
            db_fields = fields[query_field_num:-1]
            intersect_field = intersect_f(fields, query_field_num)

            # 'WILD' type database SVs match any query SV, do breakpoint annotation
            if intersect_field.db_svtype == "WILD":
                tra_compare_result.setdefault(intersect_field.query_svid,
                        []).append(db_fields)
            elif intersect_field.db_svtype == "TRA":
                tra_main_id = "_".join(
                        intersect_field.query_svid.split("_")[:-1])
                tra_pair.setdefault(tra_main_id,[]).append(fields)

    if tra_pair:
        for i in tra_pair:

            # print(tra_pair[i])
            # get paired database traslocation
            db_tra_pair = {}
            for r in tra_pair[i]:
                query_fields = r[:query_field_num]
                db_fields = r[query_field_num:-1]
                intersect_field = intersect_f(r, query_field_num)

                # print(intersect_field.db_svid)
                db_tra_main_id = "_".join(intersect_field.db_svid.split("_")[:-1])
                db_tra_pair.setdefault(db_tra_main_id,[]).append(r)
            # print(db_tra_pair)

            for p in db_tra_pair:
                if len(db_tra_pair[p]) == 2:
                    if (db_tra_pair[p][0][4] != db_tra_pair[p][1][4]):
                        # two query sv ids are not the same
                        tra_compare_result.setdefault(db_tra_pair[p][0][4],
                                []).append(db_tra_pair[p][0][query_field_num:-1])
                        tra_compare_result.setdefault(db_tra_pair[p][1][4],
                                []).append(db_tra_pair[p][1][query_field_num:-1])
    return tra_compare_result


def main():
    # test
    result = TraCompare("bedtools", sys.argv[1], sys.argv[2], 1000, "tmp1", "test_tra1", "TEST1")
    for i in result:
        print(i, ":", result[i])

if __name__ == "__main__":
    main()
