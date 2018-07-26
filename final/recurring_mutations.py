import pandas as pd

p_train = pd.read_json("p_train_features.json", orient="records")
b_train = pd.read_json("b_train_features.json", orient="records")
p_test = pd.read_json("p_test_features.json", orient="records")
b_test = pd.read_json("b_test_features.json", orient="records")
all_path = pd.concat([p_train, p_test], axis=1)
all_benign = pd.concat([b_train, b_test], axis=1)

NCBI_id = "NP_689509.1"

recur_p = []
recur_b = []

# for key in all_path:
#     if all_path[key]["count"] > 1:
#         print(NCBI_id, all_path[key]["pos"], all_path[key]["ref"], all_path[key]["alt"])
#         if str(all_path[key]["pos"]) not in recur_p:
#             recur_p.append(str(all_path[key]["pos"]))
#
# print(NCBI_id, "810", "I", "N")
#
# recur_p.append("810")
# print(",".join(recur_p))

# for var in all_benign:
#     print(all_benign[var].keys())
        #print(all_benign[var]["pos"])
        #recur_b.append(str(all_path[var]["pos"]))

#print(",".join(recur_b))

for key in b_train.keys():
    if b_train[key]["count"] > 1:
        #print(NCBI_id, b_train[key]["pos"], b_train[key]["ref"], b_train[key]["alt"])
        recur_b.append(str(b_train[key]["pos"]))

for key in b_test.keys():
    if b_test[key]["count"] > 1:
        #print(NCBI_id, b_test[key]["pos"], b_test[key]["ref"], b_test[key]["alt"])
        recur_b.append(str(b_test[key]["pos"]))


print(len(recur_b))
print(",".join(recur_b))
