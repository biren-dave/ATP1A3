import csv
import json
import random
from math import ceil

NCBI_id = "NP_689509.1"

def csv_reader(infile):
    data = []
    with open(infile, "r") as fh:
        for row in csv.reader(fh):
            data.append(row)
    return data

def csv_dreader(infile):
    data = []
    with open(infile, "r") as fh:
        for row in csv.DictReader(fh):
            data.append(row)
    return data

d = csv_reader("missense_de_novos.csv")
benigns = csv_reader("missense_benign.csv")
b_data = csv_reader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/raw_data/benign_variants_formatted.csv")
p_data = csv_reader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/raw_data/de_novos.csv")

def get_info(pos, ref, alt):
    for i in b_data[1:]:
        if (int(i[0]) == int(pos)) and (i[1] == ref) and (i[2] == alt):
            return (int(i[4]), int(i[11]))


def get_pos(pos, ref, alt):
    for i in p_data[1:]:
        if (int(i[0]) == pos) and (i[1] == ref) and (i[2] == alt):
            return int(i[4])
    return "error"

de_novos = []
a = []

for i in d:
    a.append(i[0])
    if (i not in de_novos) and (i[0] != "F266S") and (i[0] != "G734E"):
        de_novos.append(i)

random.seed(1)
random.shuffle(de_novos)
random.shuffle(benigns)

p_train = de_novos[:ceil(len(de_novos)/2)]
p_test = de_novos[ceil(len(de_novos)/2):]
b_train = benigns[:ceil(len(benigns)/2)]
b_test = benigns[ceil(len(benigns)/2):]

print(len(de_novos), len(p_train), len(p_test))
print(len(benigns), len(b_train), len(b_test))

p_tr, p_ts = {}, {}
b_tr, b_ts = {}, {}

p_ts_csv = []
b_ts_csv = []

count = 1

for i in p_train:
    p_tr[i[0]] = {"id": count, "pos": int(i[0][1:-1]), "ref": i[0][0], "alt": i[0][-1], "target": "pathogenic", "count": a.count(i[0]), "n_pos": get_pos(int(i[0][1:-1]), i[0][0], i[0][-1])}
    count += 1

for i in p_test:
    p_ts[i[0]] = {"id": count, "pos": int(i[0][1:-1]), "ref": i[0][0], "alt": i[0][-1], "target": "pathogenic", "count": a.count(i[0]), "n_pos": get_pos(int(i[0][1:-1]), i[0][0], i[0][-1])}
    p_ts_csv.append([NCBI_id, int(i[0][1:-1]), i[0][0], i[0][-1]])
    count += 1

for i in b_test:
    b_tr[i[1] + i[0] + i[2]] = {"id": count, "pos": int(i[0]), "ref": i[1], "alt": i[2], "target": "benign", "count": get_info(int(i[0]), i[1], i[2])[1], "n_pos": get_info(int(i[0]), i[1], i[2])[0]}
    count += 1

for i in b_train:
    b_ts[i[1] + i[0] + i[2]] = {"id": count, "pos": int(i[0]), "ref": i[1], "alt": i[2], "target": "benign", "count": get_info(int(i[0]), i[1], i[2])[1], "n_pos": get_info(int(i[0]), i[1], i[2])[0]}
    b_ts_csv.append([NCBI_id, int(i[0]), i[1], i[2]])
    count += 1

with open("pathogenic_train.json", "w") as fh:
    json.dump(p_tr, fh, indent=4)

with open("pathogenic_test.json", "w") as fh:
    json.dump(p_ts, fh, indent=4)

with open("benign_train.json", "w") as fh:
    json.dump(b_tr, fh, indent=4)

with open("benign_test.json", "w") as fh:
    json.dump(b_ts, fh, indent=4)

with open("pathogenic_test.csv", "w") as fh:
    csv.writer(fh, delimiter=" ").writerows(p_ts_csv)

with open("benign_test.csv", "w") as fh:
    csv.writer(fh, delimiter=" ").writerows(b_ts_csv)

print(a.count("D801N"))
