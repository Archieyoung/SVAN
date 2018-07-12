"""
Author: ArchieYoung <yangqi2@grandomics.com>
Time: Thu Jul  5 09:24:07 CST 2018
"""

class intersect_f(object):
    def __init__(self, fields, query_field_num):
        (self.query_chrom, self.query_start, self.query_end, self.query_svtype,
                self.query_svid) = fields[:5]
        (self.db_chrom, self.db_start, self.db_end, self.db_svtype,
                self.db_svid) = fields[query_field_num:query_field_num+5]
