import sys
import pandas as pd

helper = "call by ~ <metadata> <dir> (outname)\n"
helper += "\tget the OTU-like table using samples listed in <metadata> from files in <dir>"
helper = "output the OTU-like tbale as CSV table"

def filterClade(name):
    cont = name.split('|')
    return cont[-1][0] == 's'

def renameClade(name):
    cont = name.split('|')
    return ";".join(cont)

if len(sys.argv) == 3:
    dname = sys.argv[1]
    fdir = sys.argv[2]
    outname = f'{fdir}.count.csv'
elif len(sys.argv) == 4:
    dname = sys.argv[1]
    fdir = sys.argv[2]
    outname = sys.argv[3]
else:
    print(helper)
    exit()

sep = '\t'
metadata = pd.read_csv(dname)
samples = metadata.iloc[:,0]

table = []
for sid in samples:
    fname = f"{fdir}/{sid}.tsv"
    tmp_data = {'name':sid}
    with open(fname) as fin:
        for l in fin:
            if l[0] == '#':
                continue
            cont = l.split()
            if filterClade(cont[0]):
                rname = renameClade(cont[0])
                tmp_data[rname] = cont[4]
        table.append(tmp_data)

df = pd.DataFrame(table).fillna(0)
df.T.to_csv(outname, header=None)
