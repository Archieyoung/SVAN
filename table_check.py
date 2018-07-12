"""
Author: ArchieYoung <yangqi2@grandomics.com>
Time: Thu Jul  5 09:24:07 CST 2018
"""


def TableCheck(table):
    """
    check if the table have consistent field number of each row
    return field number of the table
    """
    with open(table,"r") as io:
        # skip comment lines
        check_momment_line = io.readline()
        while check_momment_line[0] == "#":
            # comment lines start with "#"
            check_momment_line = io.readline()
            pass
        first_record_line = check_momment_line # first record line
        first_record_field_num = len(first_record_line.split("\t"))
        lines = io.readlines()
        for line in lines:
            field_num = len(line.split("\t"))
            if field_num != first_record_field_num:
                raise RuntimeError("field number not consistent: "
                        "first record field num is {}, but {} "
                        "field num is {}".format(first_record_field_num,
                            field_num,line))
    return first_record_field_num
