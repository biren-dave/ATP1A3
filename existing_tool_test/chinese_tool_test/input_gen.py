import json
import csv

def json_reader(f_path):
    with open(f_path, "r") as fh:
        return json.load(fh)

def csv_writer(f_path, l):
    with open(f_path, "w") as fh:
        for i in l:
            fh.write("{}\n".format(i))

p_test = json_reader("/home/biren_dave/Documents/ATP1A3/existing_tool_test/pathogenic_test.json")
b_test = json_reader("/home/biren_dave/Documents/ATP1A3/existing_tool_test/benign_test.json")

p_pos, b_pos = [], []

for key in p_test.keys():
    p_pos.append(str(p_test[key]["pos"]))

for key in b_test.keys():
    b_pos.append(str(b_test[key]["pos"]))

csv_writer("/home/biren_dave/Documents/ATP1A3/existing_tool_test/chinese_tool_test/b_input.tsv", b_pos)
csv_writer("/home/biren_dave/Documents/ATP1A3/existing_tool_test/chinese_tool_test/p_input.tsv", p_pos)
