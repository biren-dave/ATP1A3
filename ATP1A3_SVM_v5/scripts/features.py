import json
from tqdm import tqdm

def json_reader(file_path):
    with open(file_path, "r") as fh:
        return json.load(fh)

def json_writer(d, file_path):
    with open(file_path, "w") as fh:
        json.dump(d, fh, indent=4)

p_train = json_reader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v5/pathogenic_train.json")
p_test = json_reader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v5/pathogenic_test.json")
b_train = json_reader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v5/benign_train.json")
b_test = json_reader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v5/benign_test.json")

distances = json_reader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v5/distances.json")
consurf = json_reader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v5/consurf.json")
p_cons = json_reader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v5/p_conservation.json")
n_cons = json_reader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v5/n_conservation.json")
sol_acc = json_reader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v5/solvent_acc.json")

pfam_domains = {range(33,102): 1, # cation_ATPase_N
                range(153,345): 1, #"E1-E2_ATPase
                range(416,512): 1, # cation_ATPase
                range(789,999): 1} # cation_ATPase_C

properties = {"R": 0, "K": 0, "D": 0, "E": 0, # charged
              "Q": 1, "N": 1, "H": 1, "S": 1, # polar
              "T": 1, "Y": 1, "C": 1, "W": 1,
              "A": 2, "I": 2, "L": 2, "M": 2, # hydrophobic
              "F": 2, "V": 2, "P": 2, "G": 2}

exons = {range(1,161): 1, range(161,248): 2, range(248,308): 3,
         range(308,512): 4, range(512,625): 5, range(625,761): 6,
         range(761,879): 7, range(879,1148): 8, range(1148,1347): 9,
         range(1347,1457): 10, range(1457,1592): 11, range(1592,1785): 12,
         range(1785,1961): 13, range(1961,2098): 14, range(2098,2249): 15,
         range(2249,2418): 16, range(2418,2573): 17, range(2573,2679): 18,
         range(2679,2843): 19, range(2843,2974): 20, range(2974,3076): 21,
         range(3076,3168): 22, range(3168,3552): 23}

def distance(i, j):
    for pair in distances.keys():
        i_p = int(pair.split()[0])
        j_p = int(pair.split()[1])
        if ((i, j) == (i_p, j_p)) or ((i, j) == (j_p, i_p)):
            return distances[pair]
    return 0

def consurf_score(p_pos):
    for pos in consurf.keys():
        if p_pos == int(pos):
            return consurf[pos]["normalized"]

def solvent_accessibility(p_pos):
    for pos in sol_acc.keys():
        if p_pos == int(pos):
            return sol_acc[pos]
    return 0

def p_conservation(p_pos):
    for pos in p_cons.keys():
        if p_pos == int(pos):
            return p_cons[pos]

def n_conservation(n_pos):
    for pos in n_cons.keys():
        if n_pos == int(pos):
            return n_cons[pos]

def get_domain(pos, d=pfam_domains):
    for key in d.keys():
        if pos == key:
            return d[key]
        elif (type(key) == range) and (pos in key):
            return d[key]
    return 0

def chem_prop(residue):
    return properties[residue]

def get_exon(pos, d=exons):
    for key in d.keys():
        if pos in key:
            return d[key]
    return 0

p_combined = dict(p_train)
p_combined.update(p_test)

p_total = 0

for key in p_combined.keys():
    p_total += p_combined[key]["count"]

b_combined = dict(b_train)
b_combined.update(b_test)

b_total = 0

for key in b_combined.keys():
    b_total += b_combined[key]["count"]

def features(var_json):
    for var in tqdm(var_json.keys()):

        pos = var_json[var]["pos"]

        var_json[var]["d_to_active_site"] = distance(369, pos)
        var_json[var]["d_to_ATP_site"] = distance(480, pos)
        var_json[var]["d_to_mg_site_1"] = distance(710, pos)
        var_json[var]["d_to_mg_site_2"] = distance(714, pos)
        var_json[var]["consurf"] = consurf_score(pos)
        var_json[var]["p_con"] = p_conservation(pos)
        var_json[var]["n_con"] = n_conservation(var_json[var]["n_pos"])
        var_json[var]["sol_acc"] = solvent_accessibility(pos)
        var_json[var]["ref_prop"] = chem_prop(var_json[var]["ref"])
        var_json[var]["alt_prop"] = chem_prop(var_json[var]["alt"])
        var_json[var]["exon"] = get_exon(var_json[var]["n_pos"])
        var_json[var]["pfam_domain"] = get_domain(pos)

        if var_json[var] == "pathogenic":
            var_json[var]["rel_freq"] = var_json[var]["count"] / p_total
        else:
            var_json[var]["rel_freq"] = var_json[var]["count"] / b_total

features(p_train)
features(p_test)
features(b_train)
features(b_test)

json_writer(p_train, "/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v5/features/p_train_features.json")
json_writer(p_test, "/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v5/features/p_test_features.json")
json_writer(b_train, "/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v5/features/b_train_features.json")
json_writer(b_test, "/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v5/features/b_test_features.json")
