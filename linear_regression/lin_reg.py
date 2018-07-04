from sklearn import linear_model
import csv
import numpy as np

def csv_dreader(infile):
    data = []
    with open(infile, "r") as fh:
        for row in csv.DictReader(fh):
            data.append(row)
    return data

training_data = csv_dreader("/home/biren_dave/Documents/ATP1A3/linear_regression/data/training.csv")

X = []
y = []

for v in training_data:
    #X.append([int(v["p_pos"]), int(v["n_pos"]), float(v["p_con"]), float(v["n_con"]),int(v["ref_aa"]), int(v["alt_aa"]), int(v["ref_n"]), int(v["alt_n"]), int(v["domain"]), int(v["exon"])])
    X.append([float(v["p_con"]), float(v["n_con"])])
    if v["target"] == "pathogenic":
        y.append(0)
    elif v["target"] == "benign":
        y.append(1)

lm = linear_model.LinearRegression()
lm.fit(X, y)

print(lm.score(X, y))
