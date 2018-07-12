"""
Author: ArchieYoung <yangqi2@grandomics.com>
Time: Thu Jul  5 09:24:07 CST 2018
"""


"""
db fields "chrom","start","end","svtype","annot1","annot2","annot3"...
format: "annot1;annot1","annot2,annot2";"annot3;annot3","chrom:start-end,svtype;..."
if multiple result is founded in database, write all feature in one column,
and seperate them by semicolon
"""
def format_result_pub_db(result):
    for i in result:
        variations = []
        #init annotations
        annotations = [[] for p in range(len(result[i][0])-4)]
        for j in result[i]:
            variations.append("{}:{}-{},{}".format(j[0],j[1],j[2],j[3]))
            for k in range(len(annotations)):
                annotations[k].append(j[k+4])
        variations_str = ";".join(variations)
        annotations_str = [";".join(l) for l in annotations]
        result[i] = annotations_str+[variations_str]
    return result

