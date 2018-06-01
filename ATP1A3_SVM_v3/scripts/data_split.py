import csv
import random
import math

def csv_reader(infile):
    data = []
    with open(infile, "r") as fh:
        reader = csv.reader(fh)
        next(reader, None)
        for row in reader:
            data.append(row)
    return data

def list_split(l):
    midpoint = int(math.ceil(len(l) / 2))
    l1 = l[:midpoint]
    l2 = l[midpoint:]
    return (l1, l2)

def csv_output(outfile, l):
    with open(outfile, "w") as fh:
        writer = csv.writer(fh)
        writer.writerow(["AA_pos", "ref_AA", "alt_AA", "nuc_pos", "ref_nuc", "alt_nuc", "count"])
        writer.writerows(l)

pathogenic_vars = csv_reader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/raw_data/de_novos.csv")
benign_vars = csv_reader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/raw_data/benign_variants_formatted.csv")

random.seed(1)
random.shuffle(pathogenic_vars)
random.shuffle(benign_vars)

pathogenic_train = list_split(pathogenic_vars)[0]
pathogenic_test = list_split(pathogenic_vars)[1]

benign_train = list_split(benign_vars)[0]
benign_test = list_split(benign_vars)[0]

csv_output("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/training_and_testing_data/benign_train.csv", benign_train)
csv_output("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/training_and_testing_data/benign_test.csv", benign_test)
csv_output("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/training_and_testing_data/pathogenic_train.csv", pathogenic_train)
csv_output("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/training_and_testing_data/pathogenic_test.csv", pathogenic_test)
