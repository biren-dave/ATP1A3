import csv
from matplotlib import pyplot as plt
from sklearn import svm

def csv_dreader(infile):
    data = []
    with open(infile, "r") as fh:
        for row in csv.DictReader(fh):
            data.append(row)
    return data

# def f_importances(coef, names):
#     imp = coef
#     imp,names = zip(*sorted(zip(imp,names)))
#     plt.barh(range(len(names)), imp, align='center')
#     plt.yticks(range(len(names)), names)
#     plt.show()

training_data = csv_dreader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v4/features/training.csv")

X, Xp, Xb = [], [], []
y, yp, yb = [], [], []

for v in training_data:
    X.append([int(v["p_pos"]), int(v["n_pos"]), float(v["p_con"]), float(v["n_con"]),int(v["ref_aa"]), int(v["alt_aa"]), int(v["ref_n"]), int(v["alt_n"]), int(v["domain"]), int(v["exon"])])
    #X.append([float(v["p_con"]), float(v["n_con"])])
    if v["target"] == "pathogenic":
        y.append(0)
        Xp.append(float(v["p_con"]))
        yp.append(float(v["n_con"]))
    elif v["target"] == "benign":
        y.append(1)
        Xb.append(float(v["p_con"]))
        yb.append(float(v["n_con"]))

features_names = ["p_pos", "n_pos", "p_con", "n_con", "ref_aa", "alt_aa", "ref_n", "alt_n", "domain", "exon"]
svm = svm.SVC(kernel='linear')
svm.fit(X, y)

# f_importances(svm.coef_, features_names)

#print(svm.coef_)

plt.scatter(Xp, yp, color="black", zorder=10)
plt.scatter(Xb, yb, color="white", edgecolors="black")
plt.show()
