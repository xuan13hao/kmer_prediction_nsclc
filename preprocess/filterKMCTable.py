import os
import sys
os.environ['OPENBLAS_NUM_THREADS'] = '1'
sys.path.append("/home/ben0522/toolkit/pyLib")
import pandas as pd
import numpy as np

if len(sys.argv) == 3:
    inName = sys.argv[1]
    outName = sys.argv[2]
    ratio = 0.6
elif len(sys.argv) == 4:
    inName = sys.argv[1]
    outName = sys.argv[2]
    ratio = int(sys.argv[3]) / 100.0
else:
    print("call by ~ <inName> <outName> {<ratio>}")
    exit()

ID2sample = []
sample_group = []
CB_MAPPING = {"PR":1, "CR":1, "PD":0, "DD":0, "SD":1}
with open("/home/ben0522/Data/wgs/prediction/3dataset/metadata.csv") as fin:
    for line in fin:
        if not line[0] == '#':
            cont = line.strip().split(',')
            ID2sample.append(cont[0])
            sample_group.append(CB_MAPPING[cont[5]])

TOTAL_YES = np.sum(sample_group)
TOTAL_NO = len(sample_group) - TOTAL_YES
print(TOTAL_YES, TOTAL_NO)
def takeLine(line, ratio=0.6):
    cont = line.split()
    cnt_yes = 0
    cnt_no = 0
    for i, v in enumerate(cont[1:]):
        if i%2==0:
            sample_id = int(v)
            if sample_group[sample_id]:
                cnt_yes += 1
            else:
                cnt_no += 1
    if cnt_yes >= ratio*TOTAL_YES or cnt_no >= ratio*TOTAL_NO:
        return True
    return False

def getLine(line):
    ans = ['0']*417
    cont = line.split()
    for i, v in enumerate(cont[1:]):
        if i%2==0:
            sample_id = int(v)
        else:
            ans[sample_id] = v
    return [cont[0]] + ans

fout = open(outName+'.count.tsv', 'w')
fout.write("\t".join(['name']+ID2sample))
fout.write("\n")

with open(inName, 'r') as fin:
    for line in fin:
        if takeLine(line, ratio):
            fout.write("\t".join(getLine(line)))
            fout.write("\n")

fout.close()
