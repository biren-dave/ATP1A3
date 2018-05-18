import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn import svm
import numpy as np
import math

# read training CSVs as pandas dataframes, normalize and perform PCA
def principal_ca(filepath):

    df = pd.read_csv(filepath)

    features = ["PROVEAN", "SIFT", "PP2", "position"]
    x = df.loc[:, features].values
    y = df.loc[:, ["target"]].values

    x = StandardScaler().fit_transform(x)

    pca = PCA(n_components = 2)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents, columns = ["principal component 1", "principal component 2"])
    finalDf = pd.concat([principalDf, df[["target"]], df[["count"]]], axis = 1)

    return finalDf

df_train = principal_ca("features/PCA_training_set.csv")
df_test = principal_ca("features/PCA_testing_set.csv")

#print(pca.explained_variance_ratio_)

# generates PCA plot for training data
fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1)
ax.set_xlabel("Principal Component 1", fontsize = 15)
ax.set_ylabel("Principal Component 2", fontsize = 15)
ax.set_title("2 component PCA", fontsize = 20)

targets = ["pathogenic", "benign"]
colors = ["black", "white"]
for target, color in zip(targets,colors):
    indicesToKeep = df_train["target"] == target
    ax.scatter(df_train.loc[indicesToKeep, "principal component 1"], df_train.loc[indicesToKeep, "principal component 2"], c = color, s = 50, edgecolors = "black")
ax.legend(targets)
ax.grid()
plt.savefig("2-component PCA.jpg")

# turn normalized, PCA training data into a list
training_data = df_train.values.tolist()

X = []
y = []
weights = []
p_weights = []

# parse training data to extract features, outcome and weights
for variant in training_data:
    X.append([variant[0], variant[1]])
    if variant[2] == "benign":
        y.append(1)
    elif variant[2] == "pathogenic":
        y.append(-1)
    weights.append(variant[3])
    p_weights.append(5)

X = np.array(X)
y = np.array(y)
weights = np.array(weights)
p_weights = np.array(p_weights)

# generate decision function
def plot_decision_function(classifier, sample_weight, axis, title):
    xx, yy = np.meshgrid(np.linspace(-4, 5, 500), np.linspace(-4, 5, 500))

    Z = classifier.decision_function(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)

    axis.contourf(xx, yy, Z, alpha=0.75, cmap=plt.cm.bone)
    axis.scatter(X[:, 0], X[:, 1], c = y, s = sample_weight * 10, alpha=0.9, cmap=plt.cm.bone, edgecolors="black")

    axis.axis('off')
    axis.set_title(title)

#train classifiers with weighted and unweighted data
clf_weighted = svm.SVC()
clf_weighted.fit(X, y, weights)

clf_unweighted = svm.SVC()
clf_unweighted.fit(X, y)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
plot_decision_function(clf_unweighted, p_weights, axes[0],
                       "unweighted")
plot_decision_function(clf_weighted, weights, axes[1],
                       "weighted")
plt.savefig("decision_regions.jpg")

#testing the classifiers
testing_data = df_test.values.tolist()

X_test = []
y_test = []
w_test = []

for variant in testing_data:
    X_test.append([variant[0], variant[1]])
    if variant[2] == "benign":
        y_test.append(1)
    elif variant[2] == "pathogenic":
        y_test.append(-1)
    w_test.append(variant[3])


#print(clf_weighted.score(X_test, y_test, w_test))
#print(clf_unweighted.score(X_test, y_test, w_test))

tp, tn, fp, fn = 0, 0, 0, 0

for x, y in zip(X_test, y_test):
    pred = clf_weighted.predict([x])

    if (y == -1) and (pred == -1):
        tp += 1
    elif (y == -1) and (pred == 1):
        fp += 1
    elif (y == 1) and (pred == 1):
        tn += 1
    elif (y == 1) and (pred == -1):
        fn += 1

precision = tp / (tp + fp)
recall = tp / (tp + fn)
mcc = ((tp * tn) - (fp * fn)) / math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))

print("accuracy: {}".format((tp + tn) / (tp + tn + fp + fn)))
print("precision: {}".format(precision))
print("recall/sensitivity: {}".format(recall))
print("specificity: {}".format(tn / (fp + tn)))
print("F-score: {}".format(2 * ((precision * recall) / (precision + recall))))
print("negative predictive value: {}".format(tn / (tn + fn)))
print("matthews correlation coefficient: {}".format(mcc))
