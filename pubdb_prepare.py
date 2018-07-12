#!/usr/bin/env python3
"""
prepare SV database for annotation
convert 1000genome, DGV, dbVar SV files into bed files
"""
import sys
import gzip
import logging
import operator
import os
from glob import iglob
from datetime import date


from sv_vcf import SV


# 1000genome
class one_thousand_sv(object):
    def __init__(self,record):
        # 1000genome vcf file parse
        self.record = record
        fields = record.strip().split("\t")
        (self.chrom,self.pos1,self.id,self.ref,self.alt,self.qual,self.filter,
                self.info,self.format) = fields[:9]
        self.samples = fields[9:]
        # info dict
        self.info_dict = {}
        info_list = self.info.split(";")
        for i in info_list:
            if "=" in i:
                info_id,info_value = i.split("=")
                self.info_dict[info_id] = info_value
            else:
                self.info_dict[i] = i
        # end
        if "END" in self.info_dict:
            self.pos2 = self.info_dict["END"]
        else:
            # if can not find end in info, end = start(eg. insertion)
            self.pos2 = self.pos1

        # SVLEN
        if "SVLEN" in self.info_dict:
            self.svlen = self.info_dict["SVLEN"]
        else:
            self.svlen = "NA"

        # SVTYPE
        self.sub_svtype = self.info_dict["SVTYPE"]
        if self.sub_svtype in ["SVA","LINE1","ALU","INS"]:
            self.svtype = "INS"
        elif self.sub_svtype in ["DEL","DEL_ALU","DEL_HERV","DEL_LINE1",
                "DEL_SVA"]:
            self.svtype = "DEL"
        else:
            self.svtype = self.sub_svtype

        # allele frequency
        # multi-alleles(CNVs,0,1,2...) frequency is not considered here,
        # treated as bi-alleles(0,1) frequency
        af_populations = ["AF","EAS_AF","EUR_AF","AFR_AF","AMR_AF","SAS_AF"]
        self.AFs = [self._get_af(i) for i in af_populations]

    def _get_af(self,af_population):
        # af_population: AF=0.00698882;EAS_AF=0.0069;EUR_AF=0.0189;
        # AFR_AF=0.0;AMR_AF=0.0072;SAS_AF=0.0041;
        try:
            af = sum([float(i) for i in self.info_dict[af_population].split(
                ",")])
            af = "{:.6}".format(af)
        except:
            af = "NA"
            logging.warning('Can not find "{}" in INFO of record: {}'.format(
                af_population,self.record))
        return af

    @classmethod
    def print_bed(cls,vcf_gz,out_name):
        bed_list = []
        with gzip.open(vcf_gz,"r") as io:
            n = 0
            for line in io:
                line = line.decode("utf-8")
                if line[0] == "#":
                    continue
                
                db_svid = "1000genome{}".format(n) # make 1000genome SV id
                n += 1

                sv = one_thousand_sv(line)
                sv.pos1 = int(sv.pos1)
                bed = [sv.chrom, sv.pos1, sv.pos2, sv.svtype, db_svid,
                        sv.sub_svtype]+sv.AFs
                bed_list.append(bed)
        bed_list.sort(key = operator.itemgetter(0, 1))
        bed_lines = []
        for i in bed_list:
            i[1] = str(i[1])
            bed_lines.append("\t".join(i)+"\n")
        with open(out_name,"w") as io:
            io.writelines(bed_lines)


class dgv_gold_cnv(object):
    # dgv gff3 file parse
    def __init__(self,record):
        self.record = record
        fields = record.strip().split("\t")
        # remove "chr" prefix in chrom if it exists
        self.chrom = fields[0].replace("chr","")
        self.pos1 = fields[3]
        self.pos2 = fields[4]
        self.info_dict = {}
        for i in fields[-1].split(";"):
            if "=" in i:
                info_id,info_value = i.split("=")
                self.info_dict[info_id] = info_value
            else:
                self.info_dict[i] = i
        if self.info_dict["variant_sub_type"] == "Gain":
            self.svtype = "DUP"
        elif self.info_dict["variant_sub_type"] == "Loss":
            self.svtype = "DEL"
        else:
            raise RuntimeError('variant_sub_type can either be "Gain" or "Loss"')
        self.af = self.info_dict["Frequency"]
        self.af = str(float(self.af.replace("%",""))*0.01)
        self.sample_size = self.info_dict["num_unique_samples_tested"]

    @classmethod
    def print_bed(cls,gff3,out_name):
        bed_list = []
        with open(gff3,"r") as io:
            n = 0
            for line in io:
                if line[0] == "#":
                    continue
                sv = dgv_gold_cnv(line)

                db_svid = "dgv{}".format(n)
                n += 1

                sv.pos1 = int(sv.pos1)
                bed = [sv.chrom, sv.pos1, sv.pos2, sv.svtype, db_svid,
                        sv.af, sv.sample_size]
                bed_list.append(bed)
        bed_list.sort(key = operator.itemgetter(0, 1))
        bed_lines = []
        for i in bed_list:
            i[1] = str(i[1])
            bed_lines.append("\t".join(i)+"\n")
        with open(out_name,"w") as io:
            io.writelines(bed_lines)


class dbVar_nstd37_sv(object):
    # dbvar vcf file parse
    def __init__(self,record):
        self.record = record
        fields = record.strip().split("\t")
        (self.chrom,self.pos1,self.id,self.ref,self.alt,self.qual,self.filter,
                self.info) = fields[:8]
        # info dict
        self.info_dict = {}
        info_list = self.info.split(";")
        for i in info_list:
            if "=" in i:
                info_id,info_value = i.split("=")
                self.info_dict[info_id] = info_value
            else:
                self.info_dict[i] = i
        self.pos2 = self.info_dict["END"]
        self.svtype = self.info_dict["SVTYPE"]
        try:
            self.clnsig = self.info_dict["CLNSIG"]
        except KeyError:
            self.clnsig = "NA"
        try:
            self.pheno = self.info_dict["PHENO"]
        except KeyError:
            self.pheno = "NA"

    @classmethod
    def print_bed(cls,vcf_gz,out_name):
        bed_list = []
        with gzip.open(vcf_gz,"r") as io:
            n = 0
            for line in io:
                line = line.decode("utf-8")
                if line[0] == "#":
                    continue
                sv = dbVar_nstd37_sv(line)

                db_svid = "dbvar{}".format(n)
                n += 1

                sv.pos1 = int(sv.pos1)
                bed = [sv.chrom, sv.pos1, sv.pos2, sv.svtype, db_svid,
                        sv.clnsig, sv.pheno]
                bed_list.append(bed)
        bed_list.sort(key = operator.itemgetter(0, 1))
        bed_lines = []
        for i in bed_list:
            i[1] = str(i[1])
            bed_lines.append("\t".join(i)+"\n")
        with open(out_name,"w") as io:
            io.writelines(bed_lines)


class decipher_HI(object):
    """
    Convert decipher_HI_Predictions_Version3.bed.gz to database bed
    Huang N, Lee I, Marcotte EM, Hurles ME (2010) Characterising and Predicting Haploinsufficiency in the Human Genome. PLOS Genetics 6(10): e1001154.
    """
    def __init__(self,record):
        fields = record.strip().split("\t")
        self.chrom,self.pos1,self.pos2,self.gene_hi = fields[:4]
        # remove "chr"
        self.chrom = self.chrom.replace("chr","")
        self.svtype = "WILD" # wild means that it can match any SV type, for doing svtye-insensity annotation

    @classmethod
    def print_bed(cls,input_gz,out_name):
        bed_list = []
        with gzip.open(input_gz,"r") as io:
            io.readline() # remove header
            n = 0
            for line in io:
                line = line.decode("utf-8")
                sv = decipher_HI(line)
                sv.pos1 = int(sv.pos1)
                
                db_svid = "decipherHI{}".format(n)
                n += 1

                bed = [sv.chrom, sv.pos1, sv.pos2, sv.svtype, db_svid,
                        sv.gene_hi]
                bed_list.append(bed)
        bed_list.sort(key = operator.itemgetter(0, 1))
        bed_lines = []
        for i in bed_list:
            i[1] = str(i[1])
            bed_lines.append("\t".join(i)+"\n")
        with open(out_name,"w") as io:
            io.writelines(bed_lines)


class cosmic_cnv(object):
    """
    Convert CosmicCompleteCNA.tsv.gz(CNV) into database bed
    too many records 31723168, need refine for annotation, beta!!!
    """
    def __init__(self,record):
        fields = record.strip().split("\t")
        self.CNV_ID = fields[0]
        self.Primary_site = fields[5]
        self.Primary_histology = fields[9]
        self.svtype = fields[-4]
        if self.svtype == "gain":
            self.svtype = "DUP"
        if self.svtype == "loss":
            self.svtype = "DEL"
        sv_positions = fields[-1] # chrom:start..end
        if ":" and ".." in sv_positions:
            sp1 = sv_positions.split(":")
            sp2 = sp1[1].split("..")
            self.chrom = sp1
            self.pos1 = sp2[0]
            self.pos2 = sp2[1]
        else:
            raise RuntimeError("{} not match 'chrom:start..end'".format(
                sv_positions))

    @classmethod
    def print_bed(cls,input_gz,out_name):
        bed_list = []
        cnv_ids = []
        with gzip.open(input_gz,"r") as io:
            io.readline() # remove header
            n = 0
            for line in io:
                line = line.decode("utf-8")
                sv = cosmic_cnv(line)
                if sv.CNV_ID in cnv_ids:
                    continue # remove 'Duplicated' record. CosmicCNA store CNV considering gene informations which is not necessary here
                else:
                    cnv_ids.append(sv.CNV_ID)
                    sv.pos1 = int(sv.pos1)
                    
                    db_svid = "cosmic{}".format(n)
                    n += 1

                    bed = [sv.chrom, sv.pos1, sv.pos2, sv.svtype, db_svid,
                            sv.Primary_site, sv.Primary_histology]
                    bed_list.append(bed)
        bed_list.sort(key = operator.itemgetter(0, 1))
        bed_lines = []
        for i in bed_list:
            i[1] = str(i[1])
            bed_lines.append("\t".join(i)+"\n")
        with open(out_name,"w") as io:
            io.writelines(bed_lines)


#class cosmic_sv(object):
#    """
#    convert cosmic CosmicStructExport.tsv.gz into database bed
#    """
#    def __init__(self,record):
#        fileds = record.strip().split("\t")


def main():
    #one_thousand_sv.print_bed(sys.argv[1],sys.argv[2])
    #dgv_gold_cnv.print_bed(sys.argv[1],sys.argv[2])
    #dbVar_nstd37_sv.print_bed(sys.argv[1],sys.argv[2])
    #decipher_HI.print_bed(sys.argv[1],sys.argv[2])
    #cosmic_cnv.print_bed(sys.argv[1],sys.argv[2])
    #make_grand_sv_db(sys.argv[1], "tmp")
    pass

if __name__ == "__main__":
    main()

