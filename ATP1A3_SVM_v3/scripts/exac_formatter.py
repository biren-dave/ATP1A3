import csv

single_letter = {"gly": "G", "ala": "A", "leu": "L", "met": "M", "phe": "F",
                 "trp": "W", "lys": "K", "gln": "Q", "glu": "E", "ser": "S",
                 "pro": "P", "val": "V", "ile": "I", "cys": "C", "tyr": "Y",
                 "his": "H", "arg": "R", "asn": "N", "asp": "D", "thr": "T"}

def csv_dreader(infile):
    data = []
    with open(infile, "r") as fh:
        for row in csv.DictReader(fh):
            data.append(row)
    return data

def p_var_converter(in_var):

    ref_res = ((in_var.split(".")[-1])[:3]).lower()
    alt_res = ((in_var.split(".")[-1])[-3:]).lower()
    res_pos = (in_var.split(".")[-1])[3:-3]

    return [int(res_pos), single_letter[ref_res], single_letter[alt_res]]

def coord_conveter(exac_coord):
    return (-1 * exac_coord) + 42501650

def complement(nuc):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return comp[nuc]

benign_vars = csv_dreader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/raw_data/benign_variants.csv")
b_var_data = []
total = 0

'''
for var in benign_vars:
    if var["Annotation"] in ("missense", "synonymous"):
        total += int(var["Allele Count"])
'''

for var in benign_vars:
    if var["Annotation"] in ("missense", "synonymous"):
        temp = []
        temp.extend(p_var_converter(var["Protein Consequence"]))
        temp.extend([coord_conveter(int(var["Position"])), complement(var["Reference"]), complement(var["Alternate"])])
        temp.extend([int(var["Allele Count"])])
        b_var_data.append(temp)

with open("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/raw_data/benign_variants_formatted.csv", "w") as fh:
    writer = csv.writer(fh)
    writer.writerow(["AA_pos", "ref_AA", "alt_AA", "nuc_pos", "ref_nuc", "alt_nuc", "count"])
    writer.writerows(b_var_data)
