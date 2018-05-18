import csv

def csv_reader(infile, delim):
    data = []
    with open(infile, "r") as fh:
        for row in csv.reader(fh, delimiter = delim):
            if len(row) == 1:
                data.append(row[0])
            else:
                data.append(row)
    return data

def csv_dreader(infile, delim):
    data = []
    with open(infile, "r") as fh:
        for row in csv.DictReader(fh, delimiter = delim):
            data.append(row)
    return data

def cleaner(l):
    '''
    removes whitespace from the elements of a 2D list
    '''
    clean = []
    for row in l:
        temp = []
        for e in row:
            temp.append(e.strip())
        clean.append(temp)
    return clean


b_train_pp2 = cleaner(csv_reader("scores/b_vars_train_pp.tsv", "\t"))
b_train_ps = cleaner(csv_reader("scores/b_vars_train_provean_sift.tsv", "\t"))

print(len(b_train_ps))
