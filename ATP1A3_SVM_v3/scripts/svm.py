
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn import svm
import numpy as np
import math

def principal_ca(filepath):

    df = pd.read_csv(filepath)

    features = ["p_pos", "ref_prop", "alt_prop", "ncbi_dom", "pfam_dom", "n_pos", "ref_nuc", "alt_nuc", "gc_content"]
    x = df.loc[:, features].values
    y = df.loc[:, ["target"]].values

    x = StandardScaler().fit_transform(x)

    pca = PCA(n_components = 3)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents, columns = ["principal component 1", "principal component 2", "pc 3"])
    finalDf = pd.concat([principalDf, df[["target"]], df[["count"]]], axis = 1)

    print(pca.explained_variance_ratio_)

    return finalDf

df_train = principal_ca("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/feature_vectors/training_features.csv")
#df_test = principal_ca("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/feature_vectors/testing_features.csv")

'''
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
plt.show()
#plt.savefig("2-component PCA.jpg")
'''
