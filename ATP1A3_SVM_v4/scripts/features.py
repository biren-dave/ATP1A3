import csv
import json

pfam_domains = {range(33,102): 0, # cation_ATPase_N
                range(153,345): 1, #"E1-E2_ATPase
                range(416,512): 2, # cation_ATPase
                range(789,999): 3} # cation_ATPase_C

exons = {range(1,161): 1, range(161,248): 2, range(248,308): 3,
         range(308,512): 4, range(512,625): 5, range(625,761): 6,
         range(761,879): 7, range(879,1148): 8, range(1148,1347): 9,
         range(1347,1457): 10, range(1457,1592): 11, range(1592,1785): 12,
         range(1785,1961): 13, range(1961,2098): 14, range(2098,2249): 15,
         range(2249,2418): 16, range(2418,2573): 17, range(2573,2679): 18,
         range(2679,2843): 19, range(2843,2974): 20, range(2974,3076): 21,
         range(3076,3168): 22, range(3168,3552): 23}

bases = {"A": 1, "G": 1, # purines
         "C": 2, "T": 2} # pyramidines

properties = {"R": 0, "K": 0, "D": 0, "E": 0, # charged
              "Q": 1, "N": 1, "H": 1, "S": 1, # polar
              "T": 1, "Y": 1, "C": 1, "W": 1,
              "A": 2, "I": 2, "L": 2, "M": 2, # hydrophobic
              "F": 2, "V": 2, "P": 2, "G": 2}

def csv_dreader(infile):
    data = []
    with open(infile, "r") as fh:
        for row in csv.DictReader(fh):
            data.append(row)
    return data

def get_domain(pos, d=pfam_domains):
    for key in d.keys():
        if pos == key:
            return d[key]
        elif (type(key) == range) and (pos in key):
            return d[key]
    return 4

def get_chem_properties(residue):
    return properties[residue]

def get_base(base):
    return bases[base]

def get_exon(pos, d=exons):
    for key in d.keys():
        if pos in key:
            return d[key]
    return 0

def features(l, target):
    features = []
    for var in l:
        p_pos = int(var["aa_pos"])
        n_pos = int(var["n_pos"])
        p_con = float(var["aa_cons"])
        n_con = float(var["n_cons"])
        ref_aa = get_chem_properties(var["ref_aa"])
        alt_aa = get_chem_properties(var["alt_aa"])
        ref_n = get_base(var["ref_n"])
        alt_n = get_base(var["alt_n"])
        domain = get_domain(p_pos)
        exon = get_exon(n_pos)
        count = int(var["count"])
        if ref_aa != alt_aa:
            features.append([p_pos, n_pos, p_con, n_con, ref_aa, alt_aa, ref_n, alt_n, domain, exon, count, target])
    return features

def csv_output(outfile, l):
    with open(outfile, "w") as fh:
        writer = csv.writer(fh)
        writer.writerow(["p_pos", "n_pos", "p_con", "n_con", "ref_aa", "alt_aa", "ref_n", "alt_n", "domain", "exon", "count", "target"])
        writer.writerows(l)

path_train = csv_dreader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v4/processed_data/pathogenic_train.csv")
benign_train = csv_dreader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v4/processed_data/benign_train.csv")

path_test = csv_dreader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v4/processed_data/pathogenic_test.csv")
benign_test = csv_dreader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v4/processed_data/benign_test.csv")

training_data = features(path_train, "pathogenic") + features(benign_train, "benign")
testing_data = features(path_test, "pathogenic") + features(benign_test, "benign")

csv_output("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v4/features/training_missense.csv", training_data)
csv_output("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v4/features/testing_missense.csv", testing_data)
