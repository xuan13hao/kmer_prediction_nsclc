import sys
import pandas as pd

helper = "call by ~ <metadata> <dir>\n"
helper += "\tget the OTU-like table using samples listed in <metadata> from files in <dir>"
helper = "output the OTU-like tbale as TSV table"
if len(sys.argv) == 3:
    dname = sys.argv[1]
    fdir = sys.argv[2]
    outname = f'{fdir}.count.tsv'
elif len(sys.argv) == 4:
    dname = sys.argv[1]
    fdir = sys.argv[2]
    outname = sys.argv[3]
else:
    print(helper)
    exit()

# UproC output is CSV table for each sample
sep = ','

metadata = pd.read_csv(dname)
samples = metadata.iloc[:,0]

table = []
for sid in samples:
    fname = f"{fdir}/{sid}.csv"
    tmp_data = {'Sample':sid}
    with open(fname) as fin:
        for l in fin:
            if l[0] == '#':
                continue
            cont = l.strip().split(sep)
            tmp_data[cont[0]] = cont[1]
        table.append(tmp_data)

df = pd.DataFrame(table).fillna(0)
df.T.to_csv(outname, sep='\t', header=None)
