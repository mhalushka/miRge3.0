#!/usr/bin/python
from difflib import unified_diff, Differ
d = Differ()
mirDict=dict()
with open("human_mature_miRBase.fa") as mir:
    for mil in mir:
        mil = mil.strip()
        if '>' in mil:
            headmil = mil.replace(">","")
        else:
            mirDict[headmil] = mil

with open("bad_isomir.csv") as gi:
#with open("good_isomir.csv") as gi:
    for glin in gi:
        glin = glin.strip()
        glines = glin.split("\t")
        new_string=""
        if "." in glines[1]:
            key_miR = glines[1].split(".")[0]
        else:
            key_miR = glines[1]
        try: 
            master_seq = mirDict[key_miR]
            query_seq = glines[0]
            result = list(d.compare(master_seq, query_seq))
            re_len = len(master_seq)
            variant="0"
            for idx, bases in enumerate(result):
                if not bases.startswith("-"):
                    bases = bases.replace(" ","")
                    if bases.startswith("+"):
                        if idx < re_len:
                            variant = "1"
                        new_string += "["+ str(bases) +"]"
                    else:
                        new_string += bases
                else:
                    try:
                        if bases[idx+1].startswith("+"):
                            pass 
                        else:
                            new_string += "[-]"
                    except IndexError:
                        pass
            #print(result)
            #print(master_seq +"\n"+ new_string)
            print(glin+"\t"+ master_seq +"\t"+ new_string+"\t"+variant)
        except KeyError:
            print(glines[1]+"\t-\t-\t-") 

#bad_isomir.csv
