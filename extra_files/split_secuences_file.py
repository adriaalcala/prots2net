import csv
with open('data/secuences/aa_M7_borja.fasta') as f:
    rows = list(csv.reader(f))
print(len(rows))
count = 0
for row in rows:
    if '>' in row[0]:
        count += 1
        if count % 10000 == 0:
            print(count)
            print(count//10000)
    line = row[0]+'\n'
    with open(f'data/secuences/aaborjaM7/{count//10000}.faa', 'a') as f:
        f.writelines([line])
