"""
Author: ArchieYoung <yangqi2@grandomics.com>
Time: Thu Jul  5 09:24:07 CST 2018
"""
import os
import subprocess


from bed_intersect import intersect_f
from table_check import TableCheck


def InsCompare(bedtools, bedA, bedB, max_dist, tmp_dir, prefix, db_id):
    # result sv dict, key: query sv id, value: db fields
    intersect_result = {}

    # calculate query and database fields number
    query_field_num = TableCheck(bedA)
    # db_field_num = table_check(bedB)

    # padding
    padding = max_dist
    padding_tmp_bed = os.path.join(tmp_dir, 
            "{}.{}.ins.padding{}.tmp.bed".format(prefix, db_id, padding))
    
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
            "{}.{}.ins.intersect.tmp.bed".format(prefix, db_id))

    # bedtools intersect
    with open(intersect_tmp_bed,"w") as io:
        subprocess.run([bedtools, "intersect", "-a", padding_tmp_bed, "-b",
            bedB,"-wo"], stdout = io)

    # read intersect file
    # chrom start end svtype svid svlen re info;
    # svid is a unique identifier for the sv
    with open(intersect_tmp_bed, "r") as io:
        lines = io.readlines()
        for line in lines:
            fields = line.strip().split("\t")
            query_fields = fields[:query_field_num]
            db_fields = fields[query_field_num:-1]

            intersect_field = intersect_f(fields, query_field_num)
            
            _db_svtype = intersect_field.db_svtype

            _query_svid = intersect_field.query_svid

            # WILD database SV type matchs any query SV type
            if (_db_svtype == "INS" or _db_svtype == "WILD"):
                intersect_result.setdefault(_query_svid, []).append(db_fields)

    return intersect_result
