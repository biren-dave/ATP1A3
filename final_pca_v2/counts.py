import json

def json_reader(fp):
    with open(fp) as fh:
        return json.load(fh)

data = json_reader("b_test_features.json")

pos = []

for var in data:
    pos.append(str(data[var]["pos"]))

print(",".join(pos))
