import csv
import random

all_path = []
NCBI_id = "NP_689509.1"

with open("all_path.csv") as fh:
    next(fh)
    for row in csv.reader(fh):
        row = [NCBI_id] + row
        all_path.append(row)

random.seed(1)
random.shuffle(all_path)

m = len(all_path) // 2

p_train = all_path[:m]
p_test = all_path[m:]

for i in p_test:
    print(" ".join(i))

# print(p_train)

# with open("path_train.csv", "w") as fh:
#     writer = csv.writer(fh)
#     writer.writerows(p_train)
#
# with open("path_test.csv", "w") as fh:
#         writer = csv.writer(fh)
#         writer.writerows(p_test)
