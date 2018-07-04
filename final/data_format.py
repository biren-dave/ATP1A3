import json
import csv
from tqdm import tqdm

def json_reader(f_path):
    with open(f_path, "r") as fh:
        return json.load(fh)

def distance(i, j):
    for pair in distances.keys():
        i_p = int(pair.split()[0])
        j_p = int(pair.split()[1])
        if ((i, j) == (i_p, j_p)) or ((i, j) == (j_p, i_p)):
            return distances[pair]
    return 0

def json_writer(f_path, d):
    with open(f_path, "w") as fh:
        json.dump(d, fh, indent=4)

distances = json_reader("/home/biren_dave/Documents/ATP1A3/final/distances.json")
consurf = json_reader("/home/biren_dave/Documents/ATP1A3/final/consurf.json")

d_to_ATP_site = {}
consurf_std = {}

for pos in tqdm(range(1, 1014)):
    d_to_ATP_site[pos] = distance(pos, 480)
    consurf_std[pos] = consurf[str(pos)]["normalized"]

json_writer("/home/biren_dave/Documents/ATP1A3/final/d_to_ATP_site.json", d_to_ATP_site)
json_writer("/home/biren_dave/Documents/ATP1A3/final/consurf_std.json", consurf_std)
