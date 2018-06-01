import csv
from math import ceil

ncbi_domains = {37: 0, # phosphorylation site
                56: 0,
                range(72,75): 1, # interaction with phosphoinositide-3 kinase
                range(78,99): 2, # transmembrane region
                range(122,143): 2,
                218: 0,
                265: 0,
                range(279,299): 2,
                range(311,329): 2,
                442: 0,
                548: 0,
                range(763,783): 2,
                range(834,857): 2,
                range(909,929): 2,
                933: 3, # phosphorylation site by PKA
                range(942,961): 2,
                range(976,997): 2}

pfam_domains = {range(33,102): 0, # cation_ATPase_N
                range(153,345): 1, #"E1-E2_ATPase
                range(416,512): 2, # cation_ATPase
                range(789,999): 3} # cation_ATPase_C

properties = {"R": 0, "K": 0, "D": 0, "E": 0, # charged
              "Q": 1, "N": 1, "H": 1, "S": 1, # polar
              "T": 1, "Y": 1, "C": 1, "W": 1,
              "A": 2, "I": 2, "L": 2, "M": 2, # hydrophobic
              "F": 2, "V": 2, "P": 2, "G": 2}

bases = {"A": 0, "C": 1, "G": 2, "T": 3}

dna_seq = "/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/ATP1A3.fasta"

def csv_dreader(infile):
    data = []
    with open(infile, "r") as fh:
        for row in csv.DictReader(fh):
            data.append(row)
    return data

def get_seq(path):
    data = []
    with open(path, "r") as fh:
        data = fh.readlines()
    return data[1].strip()

def get_domain(pos, d):
    for key in d.keys():
        if pos == key:
            return d[key]
        elif (type(key) == range) and (pos in key):
            return d[key]
    return 4

def get_gc_content(seq, position, range):
    try:
        start = int(ceil(position - (0.5 * range)))
        if start > 0:
            end = start + range
            region = (seq[start:end]).lower()
            return (region.count("c") + region.count("g")) / len(region)
        elif start + (0.5 * range) > len(seq):
            region = seq[start:].lower()
            return (region.count("c") + region.count("g")) / len(region)
    except ZeroDivisionError:
        return "div by 0 error"

def get_chem_properties(residue):
    return properties[residue]

def get_base(base):
    return bases[base]

def feature_builder(l, target):
    features = []
    for var in l:
        p_pos = int(var["AA_pos"])
        ref_prop = get_chem_properties(var["ref_AA"])
        alt_prop = get_chem_properties(var["alt_AA"])
        ncbi_dom = get_domain(p_pos, ncbi_domains)
        pfam_dom = get_domain(p_pos, pfam_domains)
        n_pos = int(var["nuc_pos"])
        ref_nuc = get_base(var["ref_nuc"])
        alt_nuc = get_base(var["alt_nuc"])
        gc_content = get_gc_content(seq, n_pos, 100)
        count = int(var["count"])
        t = target
        features.append([p_pos, ref_prop, alt_prop, ncbi_dom, pfam_dom, n_pos, ref_nuc, alt_nuc, gc_content, count, t])
    return features

def csv_output(outfile, l):
    with open(outfile, "w") as fh:
        writer = csv.writer(fh)
        writer.writerow(["p_pos", "ref_prop", "alt_prop", "ncbi_dom", "pfam_dom", "n_pos", "ref_nuc", "alt_nuc", "gc_content", "count", "target"])
        writer.writerows(l)

seq = get_seq(dna_seq)

path_train = csv_dreader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/training_and_testing_data/pathogenic_train.csv")
benign_train = csv_dreader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/training_and_testing_data/benign_train.csv")

path_test = csv_dreader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/training_and_testing_data/pathogenic_test.csv")
benign_test = csv_dreader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/training_and_testing_data/benign_test.csv")

training_data = feature_builder(path_train, "pathogenic") + feature_builder(benign_train, "benign")
testing_data = feature_builder(path_test, "pathogenic") + feature_builder(benign_test, "benign")

csv_output("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/feature_vectors/training_features.csv", training_data)
csv_output("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/feature_vectors/testing_features.csv", testing_data)
