import argparse
import json
import csv
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn import svm
import matplotlib.pyplot as plt
import numpy as np

def json_reader(f_path):
    with open(f_path, "r") as fh:
        return json.load(fh)

def make_list(dict):
    l = []
    for key in dict.keys():
        l.append([key, dict[key]])
    return l

def normalize(df, col):
    x = df.loc[:, [col]].values
    x = StandardScaler().fit_transform(x)
    return x

def parse_input(positions):
    p = positions.split(",")
    return p

def get_x(p):
    x = []
    for i in p:
        x.append([norm_dist[int(i) - 1][0], consurf[i]])
    return x

def predict(clf, x, pos, c):
    result = []
    for i, p in zip(x, pos):
        prob = clf.predict_proba([i])
        prob_b = prob[0][0]
        prob_p = prob[0][1]
        if (prob_b > c) or (prob_p > c):
            pred = clf.predict([i])[0]
        else:
            pred = "unsure"
        result.append([p, pred, prob_b, prob_p])
    return result

p_train = pd.read_json("p_train_features.json", orient="records")
b_train = pd.read_json("b_train_features.json", orient="records")
training_data = pd.concat([p_train, b_train], axis=1)

p_test = pd.read_json("p_test_features.json", orient="records")
b_test = pd.read_json("b_test_features.json", orient="records")
testing_data = pd.concat([p_test, b_test], axis=1)

distances = pd.DataFrame(make_list(json_reader("d_to_ATP_site.json")), columns=["pos", "dist"])
norm_dist = normalize(distances, "dist")
consurf = json_reader("consurf_std.json")

x_train, y_train = [], []

for var in training_data:
    p = training_data[var]["pos"]
    x_train.append([norm_dist[p - 1][0], consurf[str(p)]])
    y_train.append(training_data[var]["target"])

x_test, y_test = [], []

for var in testing_data:
    p = testing_data[var]["pos"]
    x_test.append([norm_dist[p - 1][0], consurf[str(p)]])
    y_test.append(testing_data[var]["target"])

clf = svm.SVC(kernel="linear", probability=True)
clf.fit(x_train, y_train)

parser = argparse.ArgumentParser()
parser.add_argument("positions", help="Input the amino acid position(s) of hATP1A3 (expected range: 1 - 1013) to predict the pathogenicity of if mutated. If inputting more than one position, separate each with a comma (e.g. 45,67,124)", type=str)
parser.add_argument("-c", "--conf_cutoff", help="Specify the confidence threshold below which predictions are called unsure (expected range: 0.5 - 1.0).", type=float)
parser.add_argument("-o", "--output_file", help="Enter the file name of the output file.", type=str)
args = parser.parse_args()

var_set = args.positions
cutoff = args.conf_cutoff
out_file = args.output_file

parsed = parse_input(var_set)
x = get_x(parsed)

with open(out_file, "w") as fh:
    writer = csv.writer(fh)
    writer.writerow(["position", "prediction", "probability benign", "probability pathogenic"])
    writer.writerows(predict(clf, x, parsed, cutoff))
