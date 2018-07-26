import csv
import numpy as np
import matplotlib.pyplot as plt

data = []

def reader(fp):
    with open(fp) as fh:
        for row in csv.DictReader(fh):
            data.append(row)
    return data

performance = reader("performance.csv")
predictions = [[], [], [], [], []]
methods = ["Chinese", "PROVEAN", "SIFT", "PolyPhen-2", "6F PCA linear CLF"]
mutations = []

for mut in performance:

    predictions[0].append(int(mut["Chinese"]))
    predictions[1].append(int(mut["PROVEN"]))
    predictions[2].append(int(mut["SIFT"]))
    predictions[3].append(int(mut["PolyPhen-2"]))
    predictions[4].append(int(mut["6F PCA linear CLF"]))

    mutations.append(mut["ref"] + mut["pos"] + mut["alt"])

predictions = np.array(predictions)
#np.swapaxes(predictions, 0, 1)

fig, ax = plt.subplots(figsize=(20,30))
im = ax.imshow(predictions, aspect=5)

ax.set_xticks(np.arange(len(mutations)))
ax.set_yticks(np.arange(len(methods)))
ax.set_xticklabels(mutations, fontsize=10, rotation=45, ha="right", rotation_mode="anchor", va="top")
ax.set_yticklabels(methods, fontsize=14)

plt.xlabel("mutation", fontsize=20)
plt.ylabel("prediction method", fontsize=20)

plt.savefig("performance.svg")

# predictions = np.array([[1,1,1,1,1],
#                         [0,0,1,1,0],
#                         [0,1,1,1,1],
#                         [1,0,1,0,0],
#                         [1,1,1,1,1],
#                         [1,1,1,1,0],
#                         [1,1,1,1,0],
#                         [0,1,1,1,0],
#                         [1,1,1,1,1],
#                         [0,1,1,0,1],
#                         [0,1,1,0,1],
#                         [1,1,1,1,1],
#                         [1,1,1,1,1],
#                         [0,1,1,1,1],
#                         [1,1,1,1,1],
#                         [1,1,1,1,0],
#                         [1,1,1,1,0],
#                         [1,1,1,1,1],
#                         [1,1,1,1,1],
#                         [1,1,1,1,1],
#                         [1,1,1,1,1],
#                         [1,1,1,1,0],
#                         [1,0,0,1,1],
#                         [1,1,1,1,0],
#                         [1,1,1,1,1]])
#
#
# print(len(predictions))
