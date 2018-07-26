import json
import csv

def json_reader(f_path):
    with open(f_path, "r") as fh:
        return json.load(fh)

NCBI_id = "NP_689509.1"

p_test = json_reader("p_test_features.json")
b_test = json_reader("b_test_features.json")

var_list = []

for var in p_test:
    # print(NCBI_id, p_test[var]["pos"], p_test[var]["ref"], p_test[var]["alt"])
    print(p_test[var]["pos"])

print("")

for var in b_test:
    # print(NCBI_id, b_test[var]["pos"], b_test[var]["ref"], b_test[var]["alt"])
    print(b_test[var]["pos"])

print("")

print(len(p_test), len(b_test), len(p_test) + len(b_test))


# for var in b_test:
    #print(var["pos"], var["ref"], var["alt"])
