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
import random

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

master = pd.concat([d_to_ATP_site.loc[:, ["dist_ATP"]], consurf.loc[:, ["consurf"]],
         sol_acc.loc[:, ["sol_acc"]], d_to_mg_site_1.loc[:, ["dist_1"]], p_cons.loc[:, ["p_cons"]],
         d_to_mg_site_2.loc[:, ["dist_2"]]], axis=1)

x = StandardScaler().fit_transform(master.values)
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x)
std_master_pca = pd.DataFrame(principalComponents, columns=["pc_1", "pc_2"])

#print(pca.explained_variance_ratio_)

p_train = [123,755,316,923,801,811,927,706,947,772,137,597,322,839,277,771,771,815,
           810,755,274,810,773,137]
p_test = [992,220,893,154,802,371,613,742,955,888,804,333,818,801,358,715,808,773,
          140,805,755,801,370,923]

b_train = pd.read_json("b_train_features.json", orient="records")
b_test = pd.read_json("b_test_features.json", orient="records")

x_train, y_train = [], []

for var in b_train:
    p = b_train[var]["pos"]
    target = b_train[var]["target"]
    x_train.append([std_master_pca.iloc[p - 1]["pc_1"], std_master_pca.iloc[p - 1]["pc_2"]])
    y_train.append(target)

for p in p_train:
    x_train.append([std_master_pca.iloc[p - 1]["pc_1"], std_master_pca.iloc[p - 1]["pc_2"]])
    y_train.append("pathogenic")

clf = svm.SVC(kernel="linear")
clf.fit(x_train, y_train)

x_test, y_test = [], []
test_pos = []

for var in b_test:
    p = b_test[var]["pos"]
    target = b_test[var]["target"]
    x_test.append([std_master_pca.iloc[p - 1]["pc_1"], std_master_pca.iloc[p - 1]["pc_2"]])
    y_test.append(target)
    test_pos.append(p)

for p in p_test:
    x_test.append([std_master_pca.iloc[p - 1]["pc_1"], std_master_pca.iloc[p - 1]["pc_2"]])
    y_test.append("pathogenic")
    test_pos.append(p)

# for p, x, y in zip(test_pos, x_test, y_test):
#     pred = clf.predict([x])
#     d = clf.decision_function([x])
#     print(p, y, pred, d)

parser = argparse.ArgumentParser()
parser.add_argument("positions", help="Input the amino acid position(s) of hATP1A3 (expected range: 1 - 1013) to predict the pathogenicity of if mutated. If inputting more than one position, separate each with a comma (e.g. 45,67,124)", type=str)
parser.add_argument("-o", "--output_file", help="Enter the file name of the output file containing tabulated predictions and distances. Note: benign and pathogenic predictions are associated with negative and positive distance vectors, respectively.", type=str)
parser.add_argument("-f", "--output_figure", help="Enter the file name of the output graph.", type=str)
args = parser.parse_args()

queries = args.positions.split(",")
out_file = args.output_file
fig_name = args.output_figure

plot_x, plot_y = [], []
labels = []
colors = []
distances = []
predictions = []

for q in queries:
    q = int(q)
    x = [std_master_pca.iloc[q - 1]["pc_1"], std_master_pca.iloc[q - 1]["pc_2"]]
    plot_x.append(x[0])
    plot_y.append(x[1])
    labels.append(q)
    p = clf.predict([x])[0]
    d = clf.decision_function([x])
    distances.append(d)

    predictions.append([q, p, d])

    if "benign" in p:
        colors.append("green")
    elif "pathogenic" in p:
        colors.append("red")

if out_file:
    with open(out_file, "w") as fh:
        writer = csv.writer(fh)
        writer.writerow(["position", "prediction", "distance_to_decision_boundary"])
        writer.writerows(predictions)
else:
    print(tabulate(predictions, headers=["position", "prediction", "distance_to_decision_boundary"]))


def make_figure():
    fig, ax = plt.subplots()
    plt.Figure(figsize=(8, 8))

    xx, yy = np.meshgrid(np.linspace(-6, 6, 500), np.linspace(-6, 6, 500))
    Z = clf.decision_function(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)

    w = clf.coef_[0]
    e = -w[0] / w[1]
    x2 = np.linspace(-6, 6)
    y2 = e * x2 - (clf.intercept_[0]) / w[1]
    ax.plot(x2, y2, lw=5, c="black", ls="--")

    ax.contourf(xx, yy, Z, alpha=0.75, cmap=plt.cm.bone_r)
    ax.scatter(plot_x, plot_y, c=colors, s=100)
    for i, txt in enumerate(labels):
        s = str(txt) + ", d = " + str(distances[i][0])
        ax.annotate(s, (plot_x[i], plot_y[i]))
        plt.xlabel("principal component 1")
        plt.ylabel("principal component 2")
        plt.ylim(-6, 6)
        plt.xlim(-6, 6)
        plt.savefig(fig_name)
        #plt.show()

if fig_name:
    make_figure()

# # for v in b_test:
# #     p = str(b_test[v]["pos"])
# #     p_test_pos.append(p)
# #
# # print(len(p_test_pos))
# # print(",".join(p_test_pos))
#
#
# def predictions(classifier, x_test, y_test):
#
#     preds = []
#
#     for x, y in zip(x_test, y_test):
#         pred = classifier.predict([x])[0]
#         prob_b = classifier.predict_proba([x])[0][0]
#         prob_p = classifier.predict_proba([x])[0][1]
#         # preds.append({"x": x, "actual": y, "prediction": pred, "prob_b": prob_b, "prob_b": prob_p})
#         preds.append({"x": x, "actual": y, "predicted": pred, "prob_b": prob_b, "prob_p": prob_p})
#
#     # for x, y in zip(x_test, y_test):
#     #     pred = classifier.predict([x])[0]
#     #     if pred == "benign":
#     #         prob = classifier.predict_proba([x])[0][0]
#     #     elif pred == "pathogenic":
#     #         prob = classifier.predict_proba([x])[0][1]
#     #     preds.append({"x": x, "actual": y, "predicted": pred, "probability": prob})
#
#     return preds
#
# def stats(predictions):
#     tp, tn, fp, fn = 0, 0, 0, 0
#     tp_x, tn_x, fp_x, fn_x= [], [], [], []
#     for i in predictions:
#         if (i["actual"] == "pathogenic") and (i["predicted"] == "pathogenic"):
#             tp += 1
#             tp_x.append(i["x"])
#         elif (i["actual"] == "benign") and (i["predicted"] == "benign"):
#             tn += 1
#             tn_x.append(i["x"])
#         elif (i["actual"] == "benign") and (i["predicted"] == "pathogenic"):
#             fp += 1
#             fp_x.append(i["x"])
#         elif (i["actual"] == "pathogenic") and (i["predicted"] == "benign"):
#             fn += 1
#             fn_x.append(i["x"])
#
#     tp_x = np.array(tp_x)
#     fp_x = np.array(fp_x)
#     tn_x = np.array(tn_x)
#     fn_x = np.array(fn_x)
#
#     return {"tp": tp_x, "tn": tn_x, "fp": fp_x, "fn": fn_x}
#     # return {"tp": len(tp_x), "tn": len(tn_x), "fp": len(fp_x), "fn": len(fn_x), "unsure": len(u_x)}
#
#
# # def stats(predictions, t=0.4):
# #     tp, tn, fp, fn = 0, 0, 0, 0
# #     unsure = 0
# #     tp_x, tn_x, fp_x, fn_x, u_x = [], [], [], [], []
# #     for i in predictions:
# #         if i["probability"] >= t:
# #             if (i["actual"] == "pathogenic") and (i["predicted"] == "pathogenic"):
# #                 tp += 1
# #                 tp_x.append(i["x"])
# #             elif (i["actual"] == "benign") and (i["predicted"] == "benign"):
# #                 tn += 1
# #                 tn_x.append(i["x"])
# #             elif (i["actual"] == "benign") and (i["predicted"] == "pathogenic"):
# #                 fp += 1
# #                 fp_x.append(i["x"])
# #             elif (i["actual"] == "pathogenic") and (i["predicted"] == "benign"):
# #                 fn += 1
# #                 fn_x.append(i["x"])
# #         elif i["probability"] < t:
# #             unsure += 1
# #             u_x.append(i["x"])
# #
# #     tp_x = np.array(tp_x)
# #     fp_x = np.array(fp_x)
# #     tn_x = np.array(tn_x)
# #     fn_x = np.array(fn_x)
# #     u_x = np.array(u_x)
# #
# #     return {"tp": tp_x, "tn": tn_x, "fp": fp_x, "fn": fn_x, "unsure": u_x}
# #     # return {"tp": len(tp_x), "tn": len(tn_x), "fp": len(fp_x), "fn": len(fn_x), "unsure": len(u_x)}
#
# def parse_input(positions):
#     p = positions.split(",")
#     return p
#
# def get_x(p):
#     x = []
#     for i in p:
#         j = int(i)
#         x.append([std_master_pca.iloc[j - 1]["pc_1"], std_master_pca.iloc[j - 1]["pc_2"]])
#     return x
#
# def predict(clf, x, pos):
#     result = []
#     for i, p in zip(x, pos):
#         pred = clf.predict([i])[0]
#         d = clf.decision_function([i])
#         result.append([p, pred, d])
#         # prob = clf.predict_proba([i])
#         # if pred == "benign":
#         #     prob_b, prob_p = prob[0][0], prob[0][1]
#         # elif pred == "pathogenic":
#         #     prob_p, prob_b = prob[0][0], prob[0][1]
#
#         # prob_b = prob[0][0]
#         # prob_p = prob[0][1]
#
#         #result.append([p, pred, prob_b, prob_p])
#     return result
#
# # pred_y = predictions(clf, x_test, y_test)
# # test_stats = stats(pred_y)
# #
# # print("tp: ", len(test_stats["tp"]))
# # print("fp: ",len(test_stats["fp"]))
# # print("tn: ",len(test_stats["tn"]))
# # print("fn: ",len(test_stats["fn"]))
# # print("unsure: ",len(test_stats["unsure"]))
#
# parser = argparse.ArgumentParser()
# parser.add_argument("positions", help="Input the amino acid position(s) of hATP1A3 (expected range: 1 - 1013) to predict the pathogenicity of if mutated. If inputting more than one position, separate each with a comma (e.g. 45,67,124)", type=str)
# parser.add_argument("-p", "--probabilities", help="Display probabilities associated with each prediction. WARNING: probabilites are meaningless for small inputs.", type=str)
# parser.add_argument("-o", "--output_file", help="Enter the file name of the output file.", type=str)
# #parser.add_argument("-c", "--conf_cutoff", help="Specify the confidence threshold below which predictions are called unsure (expected range: 0.5 - 1.0).", type=float)
#
#
# args = parser.parse_args()
#
# var_set = args.positions
# #cutoff = args.conf_cutoff
# out_file = args.output_file
#
# # if cutoff:
# #     pass
# # else:
# #     cutoff = 0.5
#
# parsed = parse_input(var_set)
# x = get_x(parsed)
# results = predict(clf, x, parsed)
#
# p, b = 0, 0
#
# for i in results:
#     if i[1] == "pathogenic":
#         p += 1
#     elif i[1] == "benign":
#         b += 1
#
# if out_file:
#     with open(out_file, "w") as fh:
#         writer = csv.writer(fh)
#         #writer.writerow(["position", "prediction", "probability_benign", "probability_pathogenic"])
#         writer.writerow(["position", "prediction", "distance_to_decision_boundary"])
#         writer.writerows(results)
# else:
#     #print(tabulate(results, headers=["pos", "prediction", "probability_benign", "probability_pathogenic"]))
#     print(tabulate(results, headers=["position", "prediction", "distance_to_decision_boundary"]))
#     print("\npathogenic: {}".format(p))
#     print("benign: {}".format(b))
#     print("total: {}".format(p + b))
#

#
# # pred = predictions(clf, x_test, y_test)
# # print(stats(pred, 0.80))
#
#
# # x1 = [i[0] for i in x]
# # x2 = [i[1] for i in x]
# # c = []
# #
# # for i in y:
# #     if i == "pathogenic":
# #         c.append("red")
# #     elif i == "benign":
# #         c.append("green")
# #
# # plt.scatter(x1, x2, c=c)
# # plt.show()
