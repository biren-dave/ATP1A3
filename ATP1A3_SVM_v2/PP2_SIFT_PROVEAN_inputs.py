'''
This script reads benign_variants.csv (from ExAc) and pathogenic_variants.csv.
An output CSV is generated for both types of variants, which can be used to
submit as inputs for batch PROVEAN/SIFT and PolyPhen-2 queries.
'''

import csv
import statistics as stats

single_letter = {"gly": "G",
                 "ala": "A",
                 "leu": "L",
                 "met": "M",
                 "phe": "F",
                 "trp": "W",
                 "lys": "K",
                 "gln": "Q",
                 "glu": "E",
                 "ser": "S",
                 "pro": "P",
                 "val": "V",
                 "ile": "I",
                 "cys": "C",
                 "tyr": "Y",
                 "his": "H",
                 "arg": "R",
                 "asn": "N",
                 "asp": "D",
                 "thr": "T"}

def csv_reader(infile):
    data = []
    with open(infile, "r") as fh:
        for row in csv.reader(fh):
            if len(row) == 1:
                data.append(row[0])
            else:
                data.append(row)
    return data

def csv_dreader(infile):
    data = []
    with open(infile, "r") as fh:
        for row in csv.DictReader(fh):
            data.append(row)
    return data

def notation_converter(in_var):
    '''
    input: amino acid change in the form p.Aaa[pos]Bbb
    output: amino acid change in the form A[pos]B
    '''
    ref_res = ((in_var.split(".")[-1])[:3]).lower()
    alt_res = ((in_var.split(".")[-1])[-3:]).lower()
    res_pos = (in_var.split(".")[-1])[3:-3]
    out_var = single_letter[ref_res] + res_pos + single_letter[alt_res]
    return out_var

def set_maker(d):
    '''
    input: dictionary with keys mapped to numerical values
    output: a tuple of lists: the first list contains the [keys, values] of the
    dictionary, where values were greater than the median value, and the second
    list contains key, value lists where values were less than the median value
    '''
    l_1, l_2 = [], []
    for key in d.keys():
        if d[key] > stats.median(d.values()):
            l_1.append([key, d[key]])
        else:
            l_2.append([key, d[key]])
    return (l_1, l_2)

def output_writer(outfile, lists, delim):
    with open(outfile, "w") as fh:
        writer = csv.writer(fh, delimiter = delim)
        for l in lists:
            writer.writerows(l)

def formatter(l, i):
    '''
    input: list containing protein variants in the form A[pos]B, and column index
    of the variant
    output: list in the form [p_id, pos, A, B] for PROVEAN, SIFT and PolyPhen-2
    online batch queries
    '''
    out_l = []
    NCBI_id = "NP_689509.1"
    for v in l:
        ref_res = v[i][0]
        alt_res = v[i][-1]
        pos = v[i][1:-1]
        out_l.append([NCBI_id, pos, ref_res, alt_res])
    return out_l

p_var = csv_reader("pathogenic_variants.csv")
b_var = csv_dreader("benign_variants.csv")

p_var_counts = {}
b_var_counts = {}

for variant in p_var:
    p_var_counts[variant] = p_var.count(variant)

for variant in b_var:
    if variant["Annotation"] in ("missense", "synonymous"):
        b_var_counts[notation_converter(variant["Protein Consequence"])] = int(variant["Allele Count"])

p_vars_train, p_vars_test = set_maker(p_var_counts)[0], set_maker(p_var_counts)[1]
b_vars_train, b_vars_test = set_maker(b_var_counts)[0], set_maker(b_var_counts)[1]

output_writer("path_vars.csv", [p_vars_train, p_vars_test], ",")
output_writer("benign_vars.csv", [b_vars_train, b_vars_test], ",")
output_writer("path_vars_train.csv", [formatter(p_vars_train, 0)], " ")
output_writer("path_vars_test.csv", [formatter(p_vars_test, 0)], " ")
output_writer("benign_vars_train.csv", [formatter(b_vars_train, 0)], " ")
output_writer("benign_vars_test.csv", [formatter(b_vars_test, 0)], " ")
