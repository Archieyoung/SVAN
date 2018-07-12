"""
A Universal Stucture Variant VCF parsing module
tested vcf: sniffles vcf, nanosv vcf, picky vcf
shared INFO ID are: SVTYPE, END, SVLEN
RE(reads evidence): sniffles, picky; nano SV: RT(2d,template,complement)
BND shared format: N[ref:pos2[

BND format:
    N]chr6:25647927]
STAT  REF  ALT   Meaning
s1    s    t[p[  piece extending to the right of p is joined after t
s2    s    t]p]  reverse comp piece extending left of p is joined after t
s3    s    ]p]t  piece extending to the left of p is joined before t
s4    s    [p[t  reverse comp piece extending right of p is joined before t
"""
import sys
import logging


class bnd(object):
    # for nano sv BND only, must be adjusted when use in other cases
    def __init__(self,bnd_string):
        self.bnd_string = bnd_string

        # fix bnd_string, no matter what is the ref character('N','A','T','C','G','t'), replace it with 'N'
        if "[" in self.bnd_string:
            if self.bnd_string[0] == "[":
                index_2 = self.bnd_string[1:].index("[")
                self.bnd_string = self.bnd_string[:index_2+2]+"N"
            else:
                index_1 = self.bnd_string.index("[")
                self.bnd_string = "N"+self.bnd_string[index_1:]
        else:
            if self.bnd_string[0] == "]":
                index_2 = self.bnd_string[1:].index("]")
                self.bnd_string = self.bnd_string[:index_2+2]+"N"
            else:
                index_1 = self.bnd_string.index("]")
                self.bnd_string = "N"+self.bnd_string[index_1:]

        # N[chr1:123456[, N]chr1:123456], ]chr1:123456]N, [chr1:123456[N
        if self.bnd_string[:2] == "N[":
            self.stat = "s1"
            self.pos = self.bnd_string[2:-1]
        if self.bnd_string[:2] == "N]":
            self.stat = "s2"
            self.pos = self.bnd_string[2:-1]
        if self.bnd_string[0] == "]":
            self.stat = "s3"
            self.pos = self.bnd_string[1:-2]
        if self.bnd_string[0] == "[":
            self.stat = "s4"
            self.pos = self.bnd_string[1:-2]
        self.chrom = self.pos.split(":")[0]
        self.pos_num = self.pos.split(":")[1]


class SV(object):
    def __init__(self,record):
        self.record = record
        fields = record.strip().split("\t")
        # checked if ID is unique for sniffles, nanosv, picky vcf
        (self.chrom1,self.pos1,self.id,self.ref,self.alt,self.qual,self.filter,
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

        self.svtype = self.info_dict["SVTYPE"]
        # BND, svtype, pos2
        if self.svtype == "BND":
            bnd_pos = bnd(self.alt) # BND position
            self.chrom2 = bnd_pos.chrom
            self.pos2 = bnd_pos.pos_num
            if (self.chrom1 == self.chrom2 and (bnd_pos.stat == "s2"
                    or bnd_pos.stat == "s4")): # INV
                self.svtype = "INV"
            elif self.chrom1 != self.chrom2:
                self.svtype = "TRA"
            else:
                raise RuntimeError("bad line {}".format(record))
        elif self.svtype == "TRA": # sniffles TRA (BND not specified)
            self.chrom2 = self.info_dict["CHR2"]
            self.pos2 = self.info_dict["END"]
        else:
            self.chrom2 = self.chrom1
            self.pos2 = self.info_dict["END"]
            # exchange pos1 and pos2, if pos1 > pos2
            if int(self.pos1) > int(self.pos2):
                tmp = self.pos1
                self.pos1 = self.pos2
                self.pos2 = tmp
        # svlen
        try:
            self.svlen = self.info_dict["SVLEN"]
        except KeyError:
            # INS, TRA do not have SVLEN attribute
            self.svlen = "NA"
        # RE(number of read evidence)
        if "RE" in self.info_dict: # sniffles and picky
            self.re = self.info_dict["RE"]
        elif "RT" in self.info_dict: # nanosv
            # RT=2d,template,compementary; no matter what kind of reads they
            # are, add them up
            self.re = sum([int(i) for i in self.info_dict["RT"].split(",")])
            self.re = str(self.re)
        else:
            self.re = "NA"
            logging.warning("Can not get RE(support reads num) "
            "for {}".format(record))

def vcf_to_bed(vcf, out_name, filter=False, min_support_reads=0,
        chr_remove=False, sv_id_prefix=None):
    with open(vcf,"r") as io:
        lines = io.readlines()

        # chromosomes to be kept
        main_chr = (list(range(1,23)) + ["X", "Y", "MT", "M", "chrX", "chrY",
            "chrM", "chrMT"] + ["chr"+str(i) for i in range(1, 23)])

        # output bedlines
        bed_lines = []

        # check if id is unique
        id_dict = {}
        #chrom1,pos1,chrom2,pos2
        previous_sv_breakpoint = ["NA","NA","NA","NA"]
        for line in lines:
            #skip comment lines
            if line.strip()[0] == "#":
                continue
            sv = SV(line)

            # filter
            if filter:
                if sv.chrom1 not in main_chr or sv.chrom2 not in main_chr:
                    continue

                if int(sv.re) < min_support_reads:
                    continue

            # remove 'chr' in chromosome id
            if chr_remove:
                sv.chrom1 = sv.chrom1.replace("chr", "")
                sv.chrom2 = sv.chrom2.replace("chr", "")

            # rename sv id
            if sv_id_prefix:
                sv.id = sv_id_prefix + sv.id

            if sv.id not in id_dict:
                id_dict[sv.id] = 1
            else:
                raise RuntimeError("Duplicated SV ID in you VCF "
                "file {}".format(sv.id))
            sv_breakpoint = [sv.chrom1,sv.pos1,sv.chrom2,sv.pos2]
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
                    sv.id,sv.svlen,sv.re,sv.info])+"\n"
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
                    sv.id+"_1",sv.svlen,sv.re,sv.info])+"\n"
                bed_line2 = "\t".join([sv.chrom2,pos2_1,pos2_2,sv.svtype,
                    sv.id+"_2",sv.svlen,sv.re,sv.info])+"\n"
                bed_line = bed_line1+bed_line2
            bed_lines.append(bed_line)
    with open(out_name,"w") as io:
        io.writelines(bed_lines)

def main():
    #test
    vcf_to_bed(sys.argv[1],sys.argv[2])

if __name__ == "__main__":
    main()

