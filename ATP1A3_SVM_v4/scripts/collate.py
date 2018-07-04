'''
- Takes pathogenic and benign variant data and combine it with corresponding
conservation scores
- Splits pathogenic and benign variant data randomly into training and testing
sets (4 CSVs are generated)
'''
import csv
import random
import math

def csv_dreader(infile):
    data = []
    with open(infile, "r") as fh:
        for row in csv.DictReader(fh):
            data.append(row)
    return data

p_cons = csv_dreader("/home/biren_dave/Documents/ATP1A3/alignments/protein_conservation.csv")
n_cons = csv_dreader("/home/biren_dave/Documents/ATP1A3/alignments/nucleotide_conservation.csv")

p_vars = csv_dreader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/raw_data/de_novos.csv")
b_vars = csv_dreader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/raw_data/benign_variants_formatted.csv")

def collate(vars, pcons=p_cons, ncons=n_cons):
    combined = []
    for var in vars:
        temp = []
        p_pos = int(var["AA_pos"])
        n_pos = int(var["benchling_pos"])
        temp.extend([p_pos, var["ref_AA"], var["alt_AA"]])
        temp.extend([n_pos, var["ref_nuc"], var["alt_nuc"]])

        for aa in pcons:
            if p_pos == int(aa["position"]):
                temp.append(aa["conservation"])

        for nuc in ncons:
            if n_pos == int(nuc["position"]):
                temp.append(nuc["conservation"])

        temp.append(int(var["count"]))
        combined.append(temp)

    return combined

def list_split(l):
    m = int(math.ceil(len(l)/2))
    l1 = l[:m]
    l2 = l[m:]
    return (l1, l2)

def csv_output(outfile, l, headers):
    with open(outfile, "w") as fh:
        writer = csv.writer(fh)
        writer.writerow(headers)
        writer.writerows(l)

headers = ["aa_pos", "ref_aa", "alt_aa", "n_pos", "ref_n", "alt_n", "aa_cons", "n_cons", "count"]

path_vars = collate(p_vars)
benign_vars = collate(b_vars)

random.seed(1)
random.shuffle(path_vars)
random.shuffle(benign_vars)

p_train, p_test = list_split(path_vars)[0], list_split(path_vars)[1]
b_train, b_test = list_split(benign_vars)[1], list_split(benign_vars)[1]

csv_output("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v4/processed_data/pathogenic_train.csv", p_train, headers)
csv_output("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v4/processed_data/pathogenic_test.csv", p_test, headers)
csv_output("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v4/processed_data/benign_train.csv", b_train, headers)
csv_output("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v4/processed_data/benign_test.csv", b_test, headers)
