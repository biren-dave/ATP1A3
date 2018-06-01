from skbio import TabularMSA, DNA, Protein
import csv
import matplotlib.pyplot as plt

n_aln = "ATP1A3_muscle_250_n.aln"
n_query_seq = "NM_152296.5" #human ATP1A3 transcript 1

msa = TabularMSA.read(n_aln, constructor=DNA)
con_scores = msa.conservation(gap_mode="ignore")
n_seqs = msa.to_dict()

nuc_scores = {}
count = 1

for c, score in zip(str(n_seqs[n_query_seq]), con_scores):
    if c != "-":
        nuc_scores[count] = (c, score)
        count += 1

x_n = [key for key in nuc_scores.keys()]
y_n = [nuc_scores[key][1] for key in nuc_scores.keys()]

p_aln = "ATP1A3_muscle_250_p.aln"
p_query_seq = "NP_689509.1"

msa = TabularMSA.read(p_aln, constructor=Protein)
con_scores = msa.conservation(gap_mode="ignore")
p_seqs = msa.to_dict()

aa_scores = {}
p = ""
count = 1

for r, score in zip(str(p_seqs[p_query_seq]), con_scores):
    if r != "-":
        aa_scores[count] = (r, score)
        count += 1

x = [key for key in aa_scores.keys()]
y = [aa_scores[key][1] for key in aa_scores.keys()]

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

def csv_dreader(infile):
    data = []
    with open(infile, "r") as fh:
        for row in csv.DictReader(fh):
            data.append(row)
    return data

path_vars = csv_dreader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/raw_data/de_novos.csv")
b_vars = csv_dreader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/raw_data/benign_missense.csv")
silent = csv_dreader("/home/biren_dave/Documents/ATP1A3/ATP1A3_SVM_v3/raw_data/benign_silent.csv")

px_de_novo, px_b = [], []
nx_de_novo, nx_b = [], []
n_silent = []

for var in path_vars:
    px_de_novo.append(int(var["AA_pos"]))
    nx_de_novo.append(int(var["benchling_pos"]))

for var in b_vars:
    px_b.append(int(var["AA_pos"]))
    nx_b.append(int(var["benchling_pos"]))

for var in silent:
    n_silent.append(int(var["benchling_pos"]))

with open("protein_conservation.csv", "w") as fh:
    writer = csv.writer(fh)
    writer.writerow(["position", "ref", "conservation"])
    for key in aa_scores.keys():
        writer.writerow([key, aa_scores[key][0], aa_scores[key][1]])

with open("nucleotide_conservation.csv", "w") as fh:
    writer = csv.writer(fh)
    writer.writerow(["position", "ref", "conservation"])
    for key in nuc_scores.keys():
        writer.writerow([key, nuc_scores[key][0], nuc_scores[key][1]])

# Generating plots
'''
plt.subplot(2, 1, 1)
plt.plot(x, y, color=lighten_color("black"), zorder=25)
a = plt.fill_between(list(range(33,102)), 0, 1.05, color=lighten_color("green",0.6), zorder=0, label="cation_ATPase_N")
b = plt.fill_between(list(range(154,345)), 0, 1.05, color=lighten_color("red"), zorder=0, label="e1-e2_ATPase")
c = plt.fill_between(list(range(416,512)), 0, 1.05, color=lighten_color("blue"), zorder=0, label="cation_ATPase")
d = plt.fill_between(list(range(789,999)), 0, 1.05, color=lighten_color("yellow"), zorder=0, label="cation_ATPase_C")
pp = plt.scatter(px_de_novo, [aa_scores[i][1] for i in px_de_novo], color="black", zorder=35, label="pathogenic de novo nonsynonymous")
pb = plt.scatter(px_b, [aa_scores[i][1] for i in px_b], color="white", edgecolors="black", zorder=30, label="benign nonsynonymous (ExAc)")
plt.xlim(0, 1013)
plt.xlabel("amino acid position")
plt.ylim(0.35, 1.05)
plt.ylabel("evolutionary conservation")
first_leg = plt.legend(handles=[pp, pb], loc=4, facecolor="white")
ax = plt.gca().add_artist(first_leg)
plt.legend(handles=[a, b, c, d], loc=3, ncol=2)

plt.subplot(2, 1, 2)
plt.plot(x_n, y_n, color=lighten_color("black"), lw=0.5)
np = plt.scatter(nx_de_novo, [nuc_scores[i][1] for i in nx_de_novo], color="black", zorder=15, label="pathogenic de novo SNV")
nb = plt.scatter(nx_b, [nuc_scores[i][1] for i in nx_b], color="white", edgecolors="black",zorder=10, label="benign nonsynonymous SNV (ExAc)")
ns = plt.scatter(n_silent, [nuc_scores[i][1] for i in n_silent], color="blue", zorder=5, label="benign synonymous SNV (ExAc)")
plt.xlim(155, 3197)
plt.xlabel("nucleotide position")
plt.ylabel("evolutionary conservation")
plt.legend(handles=[np, nb, ns], loc=4, facecolor="white")

plt.show()
'''
