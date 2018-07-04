import csv
import matplotlib.pyplot as plt
import json

def csv_dreader(infile):
    data = []
    with open(infile, "r") as fh:
        for row in csv.DictReader(fh):
            data.append(row)
    return data

def csv_output():
    pass

consurf_data = csv_dreader("/home/biren_dave/Documents/ATP1A3/consurf/clean_output.csv")
path_vars = csv_dreader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/raw_data/de_novos.csv")
b_vars = csv_dreader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/raw_data/benign_missense.csv")
silent = csv_dreader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/raw_data/benign_silent.csv")

px_de_novo, px_b = [], []
nx_de_novo, nx_b = [], []
n_silent = []

path = []
benign = []

for var in path_vars:
    px_de_novo.append(int(var["AA_pos"]))

for var in b_vars:
    px_b.append(int(var["AA_pos"]))

x_c, y_c, z_c = [], [], []
scores = {}

for i in consurf_data:
    x_c.append(int(i["pos"]))
    y_c.append(float(i["normalized_score"]))
    z_c.append(int(i["consurf_grade"]))
    scores[int(i["pos"])] = {"score": float(i["score"]), "normalized": float(i["normalized_score"]), "grade": float(i["consurf_grade"])}

with open("consurf.json", "w") as fh:
    json.dump(scores, fh, indent=4)

def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

plt.plot(x_c, y_c, color=lighten_color("black"))
p =  plt.scatter(px_de_novo, [scores[i]["normalized"] for i in px_de_novo], color="black", zorder=10, label="pathogenic missense")
b = plt.scatter(px_b, [scores[i]["normalized"] for i in px_b], color="white", edgecolors="black", zorder=5, label="benign missense")
plt.xlim(0, 1013)
plt.ylabel("normalized consurf score")
plt.xlabel("amino acid position")
plt.legend(handles=[p, b], loc=1)
plt.savefig("consurf_normalized.png")
