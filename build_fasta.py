import csv
with open('data/others/salinibacter_all.csv') as f:
    reader = csv.reader(f)
    rows = list(reader)

with open('data/secuences/sal.faa', 'w') as f:
    for row in rows[1:]:
        f.write(f">{row[0]}\n")
        f.write(f"{row[1]}\n")
