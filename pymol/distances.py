from pymol import cmd
import json
from tqdm import tqdm

distances = {}

for i in tqdm(range(23, 1017)):
    for j in range(i + 1, 1017):
        distances["{} {}".format(i, j)] = cmd.get_distance(atom1="{}/CA".format(i), atom2="{}/CA".format(j), state=-1)

with open("distances.json", "w") as fh:
    json.dump(distances, fh, indent=4)

with open("distances.json", "r") as fh:
    dist = json.load(fh)
