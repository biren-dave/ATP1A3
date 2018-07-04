import json
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
from sklearn import svm

p_train = pd.read_json("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v5/features/p_train_features.json", orient="records")
b_train = pd.read_json("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v5/features/b_train_features.json", orient="records")

p_test = pd.read_json("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v5/features/p_test_features.json", orient="records")
b_test = pd.read_json("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v5/features/b_test_features.json", orient="records")


training_data = pd.concat([p_train, b_train], axis=1)
testing_data = pd.concat([p_test, b_test], axis=1)

def format_df(input_df):

    headers = ["pos", "d_to_active_site", "d_to_ATP_site", "d_to_mg_site_1", "d_to_mg_site_2", "consurf", "p_con", "n_con", "sol_acc", "ref_prop", "alt_prop", "exon", "pfam_domain", "rel_freq", "target"]
    d = []

    for i in [j for j in input_df.keys()]:
        temp = []
        for h in headers:
            temp.append(input_df[i][h])
        d.append(temp)

    return pd.DataFrame(data=d, columns=headers)

def standardized_pca(input_df):
    features = ["pos", "d_to_active_site", "d_to_ATP_site", "d_to_mg_site_1", "d_to_mg_site_2", "consurf", "p_con", "n_con", "sol_acc", "ref_prop", "alt_prop", "exon", "pfam_domain"]

    x = input_df.loc[:, features].values
    rel_freq = input_df.loc[:, ["rel_freq"]].values
    y = input_df.loc[:, ["target"]].values

    x = StandardScaler().fit_transform(x)

    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents, columns = ['pc_1', 'pc_2'])

    return pd.concat([principalDf, input_df[["rel_freq"]], input_df[["target"]]], axis = 1)

train_s_df = standardized_pca(format_df(training_data))
test_s_df = standardized_pca(format_df(testing_data))

def raw_pca(input_df):
    features = ["pos", "d_to_active_site", "d_to_ATP_site", "d_to_mg_site_1", "d_to_mg_site_2", "consurf", "p_con", "n_con", "sol_acc", "ref_prop", "alt_prop", "exon", "pfam_domain"]

    x = input_df.loc[:, features].values
    rel_freq = input_df.loc[:, ["rel_freq"]].values
    y = input_df.loc[:, ["target"]].values

    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents, columns = ['pc_1', 'pc_2'])

    return pd.concat([principalDf, input_df[["rel_freq"]], input_df[["target"]]], axis = 1)

train_r_df = raw_pca(format_df(training_data))
test_r_df = raw_pca(format_df(testing_data))

x_train_r = train_r_df.loc[:, ["pc_1", "pc_2"]].values
y_train_r = [i[0] for i in train_r_df.loc[:, ["target"]].values]
weights = [i[0] for i in train_r_df.loc[:, ["rel_freq"]].values]

x_test_r = test_r_df.loc[:, ["pc_1", "pc_2"]].values
y_test_r = [i[0] for i in test_r_df.loc[:, ["target"]].values]

print("learning")
'''
clf_w = svm.SVC(kernel="poly")
clf_w.fit(x_train_r, y_train_r, weights)

clf = svm.SVC(kernel="poly")
clf.fit(x_train_r, y_train_r)
'''
# CLF trained on normalized data
x_train_s = train_s_df.loc[:, ["pc_1", "pc_2"]].values
y_train_s = [i[0] for i in train_s_df.loc[:, ["target"]].values]
weights = [i[0] for i in train_s_df.loc[:, ["rel_freq"]].values]

x_test_s = test_s_df.loc[:, ["pc_1", "pc_2"]].values
y_test_s = [i[0] for i in test_s_df.loc[:, ["target"]].values]

clf_w = svm.SVC(kernel="sigmoid")
clf_w.fit(x_train_s, y_train_s, weights)

clf = svm.SVC(kernel="sigmoid")
clf.fit(x_train_s, y_train_s)


def stats(classifier, x_test, y_test):
    tp, tn, fp, fn = 0, 0, 0, 0
    for x, y in zip(x_test, y_test):
        pred = classifier.predict([x])[0]
        if (y == "pathogenic") and (pred == "pathogenic"):
            tp += 1
        elif (y == "benign") and (pred == "pathogenic"):
            fp += 1
        elif (y == "benign") and (pred == "benign"):
            tn += 1
        elif (y == "pathogenic") and (pred == "benign"):
            fn += 1
    return {"tp": tp, "fp": fp, "tn": tn, "fn": fn}

print("testing")

print("weighted ", stats(clf_w, x_test_s, y_test_s))
print("unweighted ", stats(clf, x_test_s, y_test_s))



#print(test_s_df)

# fig = plt.figure(figsize = (8,8))
# ax = fig.add_subplot(1,1,1)
# ax.set_xlabel('pc_1', fontsize = 15)
# ax.set_ylabel('pc_2', fontsize = 15)
# ax.set_title('2 component PCA', fontsize = 20)
#
# targets = ["pathogenic", "benign"]
# colors = ["r", "b"]
# for target, color in zip(targets,colors):
#     indicesToKeep = train_s_df['target'] == target
#     ax.scatter(train_s_df.loc[indicesToKeep, 'pc_1']
#                , train_s_df.loc[indicesToKeep, 'pc_2']
#                , c = color
#                , s = 50)
# ax.legend(targets)
# ax.grid()
#
# plt.show()

# fig = plt.figure(figsize = (8,8))
# ax = fig.add_subplot(1,1,1)
# ax.set_xlabel('pc_1', fontsize = 15)
# ax.set_ylabel('pc_2', fontsize = 15)
# ax.set_title('2 component PCA', fontsize = 20)
#
# targets = ["pathogenic", "benign"]
# colors = ["r", "b"]
# for target, color in zip(targets,colors):
#     indicesToKeep = train_r_df['target'] == target
#     ax.scatter(train_r_df.loc[indicesToKeep, 'pc_1']
#                , train_r_df.loc[indicesToKeep, 'pc_2']
#                , c = color
#                , s = 50)
# ax.legend(targets)
# ax.grid()
#
# plt.show()
