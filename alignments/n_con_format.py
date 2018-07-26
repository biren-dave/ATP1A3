import json

with open("n_conservation.json", "r") as fh:
    n_cons = json.load(fh)

def average(l):
    return (sum(l) / len(l))

codon_cons = {}

aa = 1

for i in range(155,3194,3):
    s = [n_cons[str(i)], n_cons[str(i + 1)], n_cons[str(i + 2)]]
    codon_cons[aa] = average(s)
    aa += 1

with open("n_con_codon.json", "w") as fh:
    json.dump(codon_cons, fh, indent=4)
