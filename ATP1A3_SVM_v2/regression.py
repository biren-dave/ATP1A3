import csv
import numpy as np
from sklearn import linear_model, datasets

def csv_dreader(infile):
    data = []
    with open(infile, "r") as fh:
        for row in csv.DictReader(fh):
            data.append(row)
    return data

b_training_data = csv_dreader("features/b_train_features.csv")
p_training_data = csv_dreader("features/p_train_features.csv")

training_data = []
weights = []

for b_var in b_training_data:
    training_data.append(([float(b_var["PROVEAN"]), float(b_var["SIFT"]), float(b_var["PolyPhen-2"]), int(b_var["Position"])], 0))
    weights.append(int(b_var["Count"]))

for p_var in p_training_data:
    training_data.append(([float(p_var["PROVEAN"]), float(p_var["SIFT"]), float(p_var["PolyPhen-2"]), int(p_var["Position"])], 1))
    weights.append(int(p_var["Count"]))

X = [v[0] for v in training_data]
y = [v[1] for v in training_data]

reg = linear_model.LinearRegression()
reg.fit(X, y, weights)

print(reg.coef_)
print(reg.get_params)
