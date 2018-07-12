"""
Author: ArchieYoung <yangqi2@grandomics.com>
Time: Thu Jul  5 09:24:07 CST 2018
"""
import os
import subprocess

from bed_intersect import intersect_f
from table_check import TableCheck

def RangeCompare(bedtools, bedA, bedB, min_overlap, tmp_dir, prefix,db_id):
    # result sv dict, key: query sv id, value: db fields
    intersect_result = {}

    # calculate query and database fields number
    query_field_num = TableCheck(bedA)
    # db_field_num = table_check(bedB)

    range_intersect_tmp_bed = os.path.join(tmp_dir
            ,prefix+"."+db_id+".intersect.tmp.bed")

    if min_overlap > 0:
        with open(range_intersect_tmp_bed,"w") as io:
            subprocess.run([bedtools,"intersect","-a",bedA,"-b",bedB,"-wo",
                "-f",str(min_overlap),"-r"],stdout=io)
    else:
        with open(range_intersect_tmp_bed,"w") as io:
            subprocess.run([bedtools,"intersect","-a",bedA,"-b",bedB,"-wo"],
                    stdout=io)

    # read intersect file
    # bedA query bed
    # chrom start end svtype svid svlen re info;
    # svid is a unique identifier for the sv
    with open(range_intersect_tmp_bed,"r") as io:
        lines = io.readlines()
        for line in lines:
            fields = line.strip().split("\t")
            query_fields = fields[:query_field_num]
            db_fields = fields[query_field_num:-1]
            
            intersect_field = intersect_f(fields, query_field_num)
            
            _query_svtype = intersect_field.query_svtype
            _db_svtype = intersect_field.db_svtype

            _query_svid = intersect_field.query_svid
            
            # DEL match DEL or CNV
            if (_query_svtype == "DEL" and (_db_svtype == "DEL" or
                _db_svtype == "CNV")):
                intersect_result.setdefault(_query_svid, []).append(db_fields)
            
                # DUP match DUP or CNV
            if (_query_svtype == "DUP" and (_db_svtype == "DUP" or
                _db_svtype == "CNV")):
                intersect_result.setdefault(_query_svid, []).append(db_fields)

            if _query_svtype == "INV" and _db_svtype == "INV":
                intersect_result.setdefault(_query_svid, []).append(db_fields)
            
            # WILD database SV type matchs any query SV type
            if _db_svtype == "WILD":
                intersect_result.setdefault(_query_svid, []).append(db_fields)

    return intersect_result

