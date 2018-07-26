import argparse
import json
import csv
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn import svm
import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate

def json_reader(f_path):
    with open(f_path, "r") as fh:
        return json.load(fh)

def make_list(dict):
    l = []
    for key in dict.keys():
        l.append([key, dict[key]])
    return l

d_to_ATP_site =  pd.DataFrame(make_list(json_reader("d_to_ATP_site.json")), columns=["pos", "dist_ATP"])
consurf = pd.DataFrame(make_list(json_reader("consurf_std.json")), columns=["pos", "consurf"])
sol_acc = pd.DataFrame(make_list(json_reader("solvent_acc.json")), columns=["pos", "sol_acc"])
d_to_mg_site_1 = pd.DataFrame(make_list(json_reader("d_to_mg_site_1.json")), columns=["pos", "dist_1"])
p_cons = pd.DataFrame(make_list(json_reader("p_conservation.json")), columns=["pos", "p_cons"])
d_to_mg_site_2 = pd.DataFrame(make_list(json_reader("d_to_mg_site_2.json")), columns=["pos", "dist_2"])
pfam_domain = pd.DataFrame(make_list(json_reader("pfam_domains.json")), columns=["pos", "pfam"])
exon = pd.DataFrame(make_list(json_reader("exons.json")), columns=["pos", "exon"])
n_cons = pd.DataFrame(make_list(json_reader("n_con_codon.json")), columns=["pos", "n_cons"])
positions = pd.DataFrame([i for i in range(1, 1014)], columns=["pos"])
d_to_active_site = pd.DataFrame(make_list(json_reader("d_to_active_site.json")), columns=["pos", "dist_active"])

# master = pd.concat([d_to_ATP_site.loc[:, ["dist_ATP"]], consurf.loc[:, ["consurf"]],
#          sol_acc.loc[:, ["sol_acc"]], d_to_mg_site_1.loc[:, ["dist_1"]], p_cons.loc[:, ["p_cons"]],
#          d_to_mg_site_2.loc[:, ["dist_2"]], pfam_domain.loc[:, ["pfam"]], exon.loc[:, ["exon"]],
#          n_cons.loc[:, ["n_cons"]], positions, d_to_active_site.loc[:, ["dist_active"]]], axis=1)

master = pd.concat([d_to_ATP_site.loc[:, ["dist_ATP"]], consurf.loc[:, ["consurf"]],
         sol_acc.loc[:, ["sol_acc"]], d_to_mg_site_1.loc[:, ["dist_1"]], p_cons.loc[:, ["p_cons"]],
         d_to_mg_site_2.loc[:, ["dist_2"]]], axis=1)

x = StandardScaler().fit_transform(master.values)
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x)
std_master_pca = pd.DataFrame(principalComponents, columns=["pc_1", "pc_2"])

p_train = pd.read_json("p_train_features.json", orient="records")
b_train = pd.read_json("b_train_features.json", orient="records")
training_data = pd.concat([p_train, b_train], axis=1)

x_train, y_train = [], []

for var in training_data:
    p = training_data[var]["pos"]
    target = training_data[var]["target"]
    x_train.append([std_master_pca.iloc[p - 1]["pc_1"], std_master_pca.iloc[p - 1]["pc_2"]])
    y_train.append(target)

p_test = pd.read_json("p_test_features.json", orient="records")
b_test = pd.read_json("b_test_features.json", orient="records")
testing_data = pd.concat([p_test, b_test], axis=1)

x_test, y_test = [], []
p_test_pos = []

# for v in b_test:
#     p = str(b_test[v]["pos"])
#     p_test_pos.append(p)
#
# print(len(p_test_pos))
# print(",".join(p_test_pos))

for var in testing_data:
    p = testing_data[var]["pos"]
    target = testing_data[var]["target"]
    x_test.append([std_master_pca.iloc[p - 1]["pc_1"], std_master_pca.iloc[p - 1]["pc_2"]])
    y_test.append(target)

#clf = svm.LinearSVC()
clf = svm.SVC(kernel="linear", probability=True)
clf.fit(x_train, y_train)

def predictions(classifier, x_test, y_test):

    preds = []

    for x, y in zip(x_test, y_test):
        pred = classifier.predict([x])[0]
        prob_b = classifier.predict_proba([x])[0][0]
        prob_p = classifier.predict_proba([x])[0][1]
        # preds.append({"x": x, "actual": y, "prediction": pred, "prob_b": prob_b, "prob_b": prob_p})
        preds.append({"x": x, "actual": y, "predicted": pred, "prob_b": prob_b, "prob_p": prob_p})

    # for x, y in zip(x_test, y_test):
    #     pred = classifier.predict([x])[0]
    #     if pred == "benign":
    #         prob = classifier.predict_proba([x])[0][0]
    #     elif pred == "pathogenic":
    #         prob = classifier.predict_proba([x])[0][1]
    #     preds.append({"x": x, "actual": y, "predicted": pred, "probability": prob})

    return preds

def stats(predictions):
    tp, tn, fp, fn = 0, 0, 0, 0
    tp_x, tn_x, fp_x, fn_x= [], [], [], []
    for i in predictions:
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

    tp_x = np.array(tp_x)
    fp_x = np.array(fp_x)
    tn_x = np.array(tn_x)
    fn_x = np.array(fn_x)

    return {"tp": tp_x, "tn": tn_x, "fp": fp_x, "fn": fn_x}
    # return {"tp": len(tp_x), "tn": len(tn_x), "fp": len(fp_x), "fn": len(fn_x), "unsure": len(u_x)}


# def stats(predictions, t=0.4):
#     tp, tn, fp, fn = 0, 0, 0, 0
#     unsure = 0
#     tp_x, tn_x, fp_x, fn_x, u_x = [], [], [], [], []
#     for i in predictions:
#         if i["probability"] >= t:
#             if (i["actual"] == "pathogenic") and (i["predicted"] == "pathogenic"):
#                 tp += 1
#                 tp_x.append(i["x"])
#             elif (i["actual"] == "benign") and (i["predicted"] == "benign"):
#                 tn += 1
#                 tn_x.append(i["x"])
#             elif (i["actual"] == "benign") and (i["predicted"] == "pathogenic"):
#                 fp += 1
#                 fp_x.append(i["x"])
#             elif (i["actual"] == "pathogenic") and (i["predicted"] == "benign"):
#                 fn += 1
#                 fn_x.append(i["x"])
#         elif i["probability"] < t:
#             unsure += 1
#             u_x.append(i["x"])
#
#     tp_x = np.array(tp_x)
#     fp_x = np.array(fp_x)
#     tn_x = np.array(tn_x)
#     fn_x = np.array(fn_x)
#     u_x = np.array(u_x)
#
#     return {"tp": tp_x, "tn": tn_x, "fp": fp_x, "fn": fn_x, "unsure": u_x}
#     # return {"tp": len(tp_x), "tn": len(tn_x), "fp": len(fp_x), "fn": len(fn_x), "unsure": len(u_x)}

def parse_input(positions):
    p = positions.split(",")
    return p

def get_x(p):
    x = []
    for i in p:
        j = int(i)
        x.append([std_master_pca.iloc[j - 1]["pc_1"], std_master_pca.iloc[j - 1]["pc_2"]])
    return x

def predict(clf, x, pos):
    result = []
    for i, p in zip(x, pos):
        pred = clf.predict([i])[0]
        d = clf.decision_function([i])
        result.append([p, pred, d])
        # prob = clf.predict_proba([i])
        # if pred == "benign":
        #     prob_b, prob_p = prob[0][0], prob[0][1]
        # elif pred == "pathogenic":
        #     prob_p, prob_b = prob[0][0], prob[0][1]

        # prob_b = prob[0][0]
        # prob_p = prob[0][1]

        #result.append([p, pred, prob_b, prob_p])
    return result

# pred_y = predictions(clf, x_test, y_test)
# test_stats = stats(pred_y)
#
# print("tp: ", len(test_stats["tp"]))
# print("fp: ",len(test_stats["fp"]))
# print("tn: ",len(test_stats["tn"]))
# print("fn: ",len(test_stats["fn"]))
# print("unsure: ",len(test_stats["unsure"]))

parser = argparse.ArgumentParser()
parser.add_argument("positions", help="Input the amino acid position(s) of hATP1A3 (expected range: 1 - 1013) to predict the pathogenicity of if mutated. If inputting more than one position, separate each with a comma (e.g. 45,67,124)", type=str)
#parser.add_argument("-p", "--probabilities", help="Display probabilities associated with each prediction. WARNING: probabilites are meaningless for small inputs.", type=str)
parser.add_argument("-o", "--output_file", help="Enter the file name of the output file.", type=str)
#parser.add_argument("-c", "--conf_cutoff", help="Specify the confidence threshold below which predictions are called unsure (expected range: 0.5 - 1.0).", type=float)


args = parser.parse_args()

var_set = args.positions
#cutoff = args.conf_cutoff
out_file = args.output_file

# if cutoff:
#     pass
# else:
#     cutoff = 0.5

parsed = parse_input(var_set)
x = get_x(parsed)
results = predict(clf, x, parsed)

p, b = 0, 0

for i in results:
    if i[1] == "pathogenic":
        p += 1
    elif i[1] == "benign":
        b += 1

if out_file:
    with open(out_file, "w") as fh:
        writer = csv.writer(fh)
        #writer.writerow(["position", "prediction", "probability_benign", "probability_pathogenic"])
        writer.writerow(["position", "prediction", "distance_to_decision_boundary"])
        writer.writerows(results)
else:
    #print(tabulate(results, headers=["pos", "prediction", "probability_benign", "probability_pathogenic"]))
    print(tabulate(results, headers=["position", "prediction", "distance_to_decision_boundary"]))
    print("\npathogenic: {}".format(p))
    print("benign: {}".format(b))
    print("total: {}".format(p + b))

xx, yy = np.meshgrid(np.linspace(-5, 5, 500), np.linspace(-4, 4, 500))
Z = clf.decision_function(np.c_[xx.ravel(), yy.ravel()])
Z = Z.reshape(xx.shape)

w = clf.coef_[0]
e = -w[0] / w[1]
x2 = np.linspace(-5, 5)  # make sure the line is long enough
y2 = e * x2 - (clf.intercept_[0]) / w[1]
plt.plot(x2, y2, lw=5, c="black", ls="--")

plt.contourf(xx, yy, Z, alpha=0.75, cmap=plt.cm.bone_r)
c = plt.scatter([std_master_pca.iloc[i - 1]["pc_1"] for i in [681,220,583,277,89,371,597,756]], [std_master_pca.iloc[i - 1]["pc_2"] for i in [681,220,583,277,89,371,597,756]], c="red", s=100, edgecolors='black', label="false negative")
d = plt.scatter([std_master_pca.iloc[i - 1]["pc_1"] for i in [817,988,941,115]], [std_master_pca.iloc[i - 1]["pc_2"] for i in [817,988,941,115]], c="blue", s=100, marker="X", edgecolors='black', label="false positive")
# a = plt.scatter(test_stats["tp"][:, 0], test_stats["tp"][:, 1], c="green", s=100, edgecolors='black', label="true positive")
# b = plt.scatter(test_stats["tn"][:, 0], test_stats["tn"][:, 1], c="green", s=100, marker="X", edgecolors='black', label="true negative")
# c = plt.scatter(test_stats["fp"][:, 0], test_stats["fp"][:, 1], c="red", s=100, edgecolors='black', label="false positive")
# d = plt.scatter(test_stats["fn"][:, 0], test_stats["fn"][:, 1], c="red", s=100, marker="X", edgecolors='black', label="false negative")
# if len(test_stats["unsure"]) > 0:
#     e = plt.scatter(test_stats["unsure"][:, 0], test_stats["unsure"][:, 1], c="yellow", s=100, marker="s", edgecolors='black', label="unsure")
#     plt.legend(handles=[a, b, c, d, e], loc=1)
# else:
#     plt.legend(handles=[a, b, c, d], loc=1)
plt.legend(handles=[c, d], loc=1)
plt.xlabel("principal component 1")
plt.ylabel("principal component 2")
plt.ylim(-4, 4)
plt.xlim(-5, 5)
plt.savefig("false_calls.png")
#plt.show()

# pred = predictions(clf, x_test, y_test)
# print(stats(pred, 0.80))


# x1 = [i[0] for i in x]
# x2 = [i[1] for i in x]
# c = []
#
# for i in y:
#     if i == "pathogenic":
#         c.append("red")
#     elif i == "benign":
#         c.append("green")
#
# plt.scatter(x1, x2, c=c)
# plt.show()
