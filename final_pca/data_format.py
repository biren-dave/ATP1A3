import json
import csv
from tqdm import tqdm

exons = {range(1,3): 1, range(3,32): 2, range(32,52): 3,
         range(52,120): 4, range(120,158): 5, range(158,203): 6,
         range(203,242): 7, range(242,332): 8, range(332,398): 9,
         range(398,435): 10, range(435,480): 11, range(480,544): 12,
         range(544,603): 13, range(603,649): 14, range(649,699): 15,
         range(699,755): 16, range(755,807): 17, range(807,848): 18,
         range(848,897): 19, range(897,941): 20, range(941,975): 21,
         range(975,1005): 22, range(1005,1014): 23}

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

def get_domain(pos):
    for key in pfam_domains.keys():
        if pos in key:
            return pfam_domains[key]
    return 0

def get_exon(pos):
    for key in exons.keys():
        if pos in key:
            return exons[key]
    return 0

pfam_domains = {range(33,102): 1, # cation_ATPase_N
                range(153,345): 1, #"E1-E2_ATPase
                range(416,512): 1, # cation_ATPase
                range(789,999): 1} # cation_ATPase_C

d_to_mg_site_1, d_to_mg_site_2, d_to_active_site = {}, {}, {}
pfam = {}
exon = {}
# pos  = {}

distances = json_reader("/home/biren_dave/Documents/ATP1A3/final/distances.json")

for pos in tqdm(range(1, 1014)):
    # d_to_mg_site_1[pos] = distance(pos, 710)
    # d_to_mg_site_2[pos] = distance(pos, 714)
    # pfam[pos] = get_domain(pos)
    exon[pos] = get_exon(pos)
    d_to_active_site[pos] = distance(pos, 369)




# json_writer("d_to_mg_site_1.json", d_to_mg_site_1)
# json_writer("d_to_mg_site_2.json", d_to_mg_site_2)
# json_writer("pfam_domains.json", pfam)

json_writer("d_to_active_site.json", d_to_active_site)
json_writer("exons.json", exon)
