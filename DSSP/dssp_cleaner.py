import csv
import json

with open("dssp_raw.csv", "r") as fh:
    data = fh.readlines()

clean = []
clean_d = {}

for i in data:
    tokens = i.split()
    try:
        clean.append([int(tokens[1]) - 3, tokens[3], tokens[-2]])
        clean_d[int(tokens[1])] = int(tokens[-2])
    except:
        pass

with open("dssp_clean.csv", "w") as fh:
    writer = csv.writer(fh)
    writer.writerow(["pos", "res", "accessibility"])
    writer.writerows(clean)

with open("solvent_acc.json", "w") as fh:
    json.dump(clean_d, fh, indent=4)
