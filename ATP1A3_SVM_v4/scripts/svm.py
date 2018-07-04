import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn import svm
import numpy as np
import math

def principal_ca(filepath):

    df = pd.read_csv(filepath)

    features = ["p_pos", "n_pos", "p_con", "n_con", "ref_aa", "alt_aa", "ref_n", "alt_n", "domain", "exon"]
    #features = ["p_con", "n_con"]
    x = df.loc[:, features].values
    y = df.loc[:, ["target"]].values

    x = StandardScaler().fit_transform(x)

    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents, columns = ["principal component 1", "principal component 2"])
    finalDf = pd.concat([principalDf, df[["target"]], df[["count"]]], axis = 1)
    #print(pca.explained_variance_ratio_)
    return finalDf

df_train = principal_ca("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v4/features/training_missense.csv")
df_test = principal_ca("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v4/features/testing_missense.csv")

fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1)
ax.set_xlabel("Principal Component 1", fontsize = 15)
ax.set_ylabel("Principal Component 2", fontsize = 15)
ax.set_title("2 component PCA", fontsize = 20)

targets = ["pathogenic", "benign"]
colors = ["black", "white"]
zorders = [10, 0]
for target, color, z in zip(targets,colors,zorders):
    indicesToKeep = df_train["target"] == target
    ax.scatter(df_train.loc[indicesToKeep, "principal component 1"], df_train.loc[indicesToKeep, "principal component 2"], c = color, s = 50, zorder=z, edgecolors = "black")
ax.legend(targets)
ax.grid()
plt.savefig("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v4/figures/pca.jpg")

training_data = df_train.values.tolist()
testing_data = df_test.values.tolist()

x_train, y_train, w_train = [], [], []
x_test, y_test, w_test = [], [], []

path_count = 0
benign_count = 0

for v in (training_data + testing_data):
    if v[2] == "pathogenic":
        path_count += 1
    elif v[2] == "benign":
        benign_count += 1

for v in training_data:
    x_train.append([v[0], v[1]])
    if v[2] == "pathogenic":
        y_train.append(-1)
        w_train.append(v[3] / path_count)
    elif v[2] == "benign":
        y_train.append(1)
        w_train.append(v[3] / benign_count)

for v in testing_data:
    x_test.append([v[0], v[1]])
    if v[2] == "pathogenic":
        y_test.append(-1)
        w_test.append(v[3] / path_count)
    elif v[2] == "benign":
        y_test.append(1)
        w_test.append(v[3] / benign_count)

x_train = np.array(x_train)
y_train = np.array(y_train)
w_train = np.array(w_train)

def plot_decision_function(classifier, sample_weight, axis, title):
    xx, yy = np.meshgrid(np.linspace(-4, 5, 500), np.linspace(-4, 5, 500))

    Z = classifier.decision_function(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)

    axis.contourf(xx, yy, Z, alpha=0.75, cmap=plt.cm.bone)

    if title == "weighted":
        axis.scatter(x_train[:, 0], x_train[:, 1], c = y_train, s = sample_weight * 1000, alpha=0.9, cmap=plt.cm.bone, edgecolors="black")
    else:
        axis.scatter(x_train[:, 0], x_train[:, 1], c = y_train, s = sample_weight * 75, alpha=0.9, cmap=plt.cm.bone, edgecolors="black")

    axis.axis('off')
    axis.set_title(title)

clf_w = svm.SVC()
clf_w.fit(x_train, y_train, w_train)

clf = svm.SVC()
clf.fit(x_train, y_train)

p_weights = np.ones(len(x_train))

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
plot_decision_function(clf, p_weights, axes[0],
                       "unweighted")
plot_decision_function(clf_w, w_train, axes[1],
                       "weighted")
plt.savefig("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v4/features/decision_regions.jpg")
#plt.show()


#print(len(x_train))
tp, tn, fp, fn = 0, 0, 0, 0

for x, y in zip(x_train, y_train):
    pred = clf.predict([x])
    print(pred, y)
    if (y == -1) and (pred == -1):
        tp += 1
    elif (y == -1) and (pred == 1):
        fp += 1
    elif (y == 1) and (pred == 1):
        tn += 1
    elif (y == 1) and (pred == -1):
        fn += 1

print("tp: {}, tn: {}, fp: {}, fn: {}".format(tp, tn, fp, fn))
print("score: {}".format(clf.score(x_test, y_test)))

try:
    precision = tp / (tp + fp)
except:
    pass
try:
    recall = tp / (tp + fn)
except:
    pass
try:
    mcc = ((tp * tn) - (fp * fn)) / math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
except:
    pass
try:
    print("accuracy: {}".format((tp + tn) / (tp + tn + fp + fn)))
except:
    pass
try:
    print("precision: {}".format(precision))
except:
    pass
try:
    print("recall/sensitivity: {}".format(recall))
except:
    pass
try:
    print("specificity: {}".format(tn / (fp + tn)))
except:
    pass
try:
    print("F-score: {}".format(2 * ((precision * recall) / (precision + recall))))
except:
    pass
try:
    print("negative predictive value: {}".format(tn / (tn + fn)))
except:
    pass
try:
    print("matthews correlation coefficient: {}".format(mcc))
except:
    pass
