import json
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn import svm
import matplotlib.pyplot as plt
import numpy as np

p_train = pd.read_json("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v5/features/p_train_features.json", orient="records")
b_train = pd.read_json("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v5/features/b_train_features.json", orient="records")

p_test = pd.read_json("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v5/features/p_test_features.json", orient="records")
b_test = pd.read_json("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v5/features/b_test_features.json", orient="records")

training_data = pd.concat([p_train, b_train], axis=1)
testing_data = pd.concat([p_test, b_test], axis=1)

features = ["pos", "d_to_active_site", "d_to_ATP_site", "d_to_mg_site_1", "d_to_mg_site_2", "consurf", "p_con", "n_con", "sol_acc", "ref_prop", "alt_prop", "exon", "pfam_domain"]

def format_df(input_df):

    headers = ["pos", "d_to_active_site", "d_to_ATP_site", "d_to_mg_site_1", "d_to_mg_site_2", "consurf", "p_con", "n_con", "sol_acc", "ref_prop", "alt_prop", "exon", "pfam_domain", "rel_freq", "target"]
    d = []

    for i in [j for j in input_df.keys()]:
        temp = []
        for h in headers:
            temp.append(input_df[i][h])
        d.append(temp)

    return pd.DataFrame(data=d, columns=headers)

def standardize(input_df):
    x = input_df.loc[:, features].values
    x = StandardScaler().fit_transform(x)
    return x

def two_feature_std(input_df):
    x = input_df.loc[:, ["d_to_ATP_site", "consurf"]]
    x = StandardScaler().fit_transform(x)
    return x

def predictions(classifier, x_test, y_test):

    preds = []

    for x, y in zip(x_test, y_test):
        pred = classifier.predict([x])[0]
        if pred == "benign":
            prob = classifier.predict_proba([x])[0][0]
        elif pred == "pathogenic":
            prob = classifier.predict_proba([x])[0][1]
        preds.append({"x": x, "actual": y, "predicted": pred, "probability": prob})

    return preds

def stats(predictions, t=0.5):
    tp, tn, fp, fn = 0, 0, 0, 0
    unsure = 0
    tp_x, tn_x, fp_x, fn_x, u_x = [], [], [], [], []
    for i in predictions:
        if i["probability"] >= t:
            if (i["actual"] == "pathogenic") and (i["predicted"] == "pathogenic"):
                tp += 1
                tp_x.append(i["x"])
            elif (i["actual"] == "benign") and (i["predicted"] == "benign"):
                tn += 1
                tn_x.append(i["x"])
            elif (i["actual"] == "benign") and (i["predicted"] == "pathogenic"):
                fp += 1
                fp_x.append(i["x"])
            elif (i["actual"] == "pathogenic") and (i["predicted"] == "benign"):
                fn += 1
                fn_x.append(i["x"])
        elif i["probability"] < t:
            unsure += 1
            u_x.append(i["x"])

    tp_x = np.array(tp_x)
    fp_x = np.array(fp_x)
    tn_x = np.array(tn_x)
    fn_x = np.array(fn_x)
    u_x = np.array(u_x)

    return {"tp": tp_x, "tn": tn_x, "fp": fp_x, "fn": fn_x, "unsure": u_x}

X = two_feature_std(format_df(training_data))
y = np.array([i[0] for i in format_df(training_data).loc[:, ["target"]].values])
t = y
weights = [i[0] for i in format_df(training_data).loc[:, ["rel_freq"]].values]
pos_train = [int(i[0]) for i in format_df(training_data).loc[:, ["pos"]].values]

colors, z = [], []

for i in y:
    if i == "pathogenic":
        colors.append("black")
        z.append(10)
    elif i == "benign":
        colors.append("white")
        z.append(5)

colors = np.array(colors)

clf_2_w = svm.SVC(kernel="linear")
clf_2_w.fit(X, y, weights)

clf_2_u = svm.SVC(kernel="linear", probability=True)
clf_2_u.fit(X, y)

x_test = two_feature_std(format_df(testing_data))
y_test = [i[0] for i in format_df(testing_data).loc[:, ["target"]].values]
pos_test = [int(i[0]) for i in format_df(testing_data).loc[:, ["pos"]].values]

for x, p in zip(X, pos_train):
    if p == 801:
        print(p, clf_2_u.predict([x]), clf_2_u.predict_proba([x]))

for x, p in zip(x_test, pos_test):
    if p == 810:
        print(p, clf_2_u.predict([x]), clf_2_u.predict_proba([x]))

pred_y = predictions(clf_2_u, x_test, y_test)
test_stats = stats(pred_y, 0.50)

# xx, yy = np.meshgrid(np.linspace(-2, 3, 500), np.linspace(-1, 3, 500))
# Z = clf_2_u.decision_function(np.c_[xx.ravel(), yy.ravel()])
# Z = Z.reshape(xx.shape)
#
# w = clf_2_u.coef_[0]
# e = -w[0] / w[1]
# x2 = np.linspace(-2, 3)  # make sure the line is long enough
# y2 = e * x2 - (clf_2_u.intercept_[0]) / w[1]
# plt.plot(x2, y2, lw=5, c="black", ls="--")
#
# plt.contourf(xx, yy, Z, alpha=0.75, cmap=plt.cm.bone_r)
# a = plt.scatter(test_stats["tp"][:, 0], test_stats["tp"][:, 1], c="green", s=100, edgecolors='black', label="true positive")
# b = plt.scatter(test_stats["tn"][:, 0], test_stats["tn"][:, 1], c="green", s=100, marker="X", edgecolors='black', label="true negative")
# c = plt.scatter(test_stats["fp"][:, 0], test_stats["fp"][:, 1], c="red", s=100, edgecolors='black', label="false positive")
# d = plt.scatter(test_stats["fn"][:, 0], test_stats["fn"][:, 1], c="red", s=100, marker="X", edgecolors='black', label="false negative")
# if len(test_stats["unsure"]) > 0:
#     e = plt.scatter(test_stats["unsure"][:, 0], test_stats["unsure"][:, 1], c="yellow", s=100, marker="s", edgecolors='black', label="unsure")
#     plt.legend(handles=[a, b, c, d, e], loc=1)
# else:
#     plt.legend(handles=[a, b, c, d], loc=1)
# plt.xlabel("distance from alpha-carbon to ATP binding site (normalized)")
# plt.ylabel("ConSurf score (normalized)")
# plt.ylim(-1, 3)
# plt.xlim(-2, 3)
# plt.savefig("clf_test_0.90.png")
# plt.show()
