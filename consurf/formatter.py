import csv

def csv_reader(infile):
    data = []
    with open(infile, "r") as fh:
        for row in csv.reader(fh):
            data.append(row)
    return data

headers = ["pos", "score", "normalized_score", "consurf_grade"]

def csv_writer(outfile, l, headers=headers):
    with open(outfile, "w") as fh:
        writer = csv.writer(fh)
        writer.writerow(headers)
        writer.writerows(l)

consurf_out = csv_reader("output.csv")

clean = []

for row in consurf_out:
    score = row[1].split(" ")[1]
    clean.append([row[0], score, row[2], row[3]])

csv_writer("clean_output.csv", clean)
