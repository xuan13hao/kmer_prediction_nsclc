import os
import sys
os.environ['OPENBLAS_NUM_THREADS'] = '16'
sys.path.append("/home/xuan/K-mer-prediction-main")

import pandas as pd
import numpy as np
import pickle

from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.feature_selection import RFECV
from sklearn.feature_selection import SelectKBest, mutual_info_classif

from sklearn.metrics import recall_score
from sklearn.metrics import precision_score
from sklearn.metrics import f1_score
from sklearn.metrics import accuracy_score

import warnings
warnings.filterwarnings('ignore')
filter = "bwa"
helper = "call by\n    ~ <pf> <ss_name> <ds> <code_stats> <profile ID>"
# pf = 'CB'

if len(sys.argv) > 5:
    pf = sys.argv[1]
    ss_name = sys.argv[2]
    ds = sys.argv[3]
    code_stats = sys.argv[4]
    selected = sys.argv[5:]
else:
    print(helper)
    exit()

pf_list = ["pfam", "kegg", "tax", "k30", "k16", "k60", "com"]
if pf not in pf_list:
    print(helper)
    print(f"you must choose <pf> from {' '.join(pf_list)}!")
    exit()

sp1 = 'CB'
sp2 = 'CB' # we are predicting response!!

if ds not in ['mix', 'dds3', 'ds1', 'ds2', 'ds3', 'ds12', 'pfs1', 'pfs2', 'pfs3']:
    print(helper)
    print("you must choose <ds> from {mix, dds3, ds1, ds2, ds3, ds12, pfs1, pfs2}!")
    print("    mix: mix all samples;")
    print("    dds3: train on ds1, ds2, and part of ds3, test on the rest of ds3;")
    print("    ds12: mix 12 together;")
    print("    ds1, ds2, ds3: only use the corresponding dataset.")
    print("    pfs1: train on 50 PFS samples, test on 79 response samples (ds12)")
    print("    pfs2: train on 50 PFS samples, test on 29 response samples (ds12)")
    print("    pfs3: train on 50 PFS samples, test on 50 PFS samples (ds12) - verify cancers results")
    exit()

if not ((not code_stats[1:] or code_stats[1:].isdigit()) and code_stats[0] in "@-0123456789" ):
    print(helper)
    print("code_stats must be in {<int>, @<int>, -}")
    print("    - for old computing;")
    print("    <int>, use <int>")
    print("    @<int>, use random int from 0 to 1000")
    exit()

if code_stats.isdigit():
    top = 1
else:
    top = int(code_stats[1:])

if pf == 'k30':
    for i,v in enumerate(selected):
        selected[i] = int(v)
elif pf == "com":
    for i,v in enumerate(selected):
        if v.isdigit():
            selected[i] = int(v)

if ds == 'pfs3':
    sp1 = 'pfs'
    sp2 = 'pfs'
"""== Parameters ============================================================"""
ff_skipNormalize = False
map_class = {
    "SR" : {"PR":1, "CR":1, "PD":0, "DD":0, "SD":0},
    "CB" : {"PR":1, "CR":1, "PD":0, "DD":0, "SD":1},
    "3C" : {"PR":2, "CR":2, "PD":0, "DD":0, "SD":1},
    "2C" : {"PR":1, "CR":1, "PD":0, "DD":0},
    "org": {"PR":2, "CR":3, "PD":0, "DD":0, "SD":1},
    "pfs": {">6":1, "<3":0},
}
base_dir   = f"uproc_{filter}/"
work_space = f"{base_dir}TestUsingGivenProfile/"
os.system("mkdir -p " + work_space) # create work_space if not exist
"""=========================================================================="""

"""== Utility ==============================================================="""
def getDS(labels):
    gt = np.zeros(len(labels))
    for (i, v) in enumerate(labels):
        if v[0] == 'J':
            gt[i] = 0
        elif v[0] == 'E':
            gt[i] = 1
        elif v[0] == 'S':
            gt[i] = 2
    return gt

def getGT(labels, tag):
    gt = np.zeros(len(labels))
    for (i, v) in enumerate(labels):
        gt[i] = map_class[tag][v[-2:]]
    return gt

def normalize_feature(df):
    ans = []
    df_std = df.std()
    df_mean = df.mean()
    for col in df.columns:
        ans.append((df[col] - df_mean[col]) / (df_std[col]+1e-17))
    return pd.concat(ans, axis=1)

def normalize_test(df, df_train):
    ans = []
    df_std = df_train.std()
    df_mean = df_train.mean()
    for col in df.columns:
        ans.append((df[col] - df_mean[col]) / (df_std[col]+1e-17))
    return pd.concat(ans, axis=1)

renameSample_test = {}
renameSample_pfs = {}
stable_list = []
with open("metadata.csv") as fin:
    for line in fin:
        if not line[0] == '#':
            cont = line.strip().split(',')
            renameSample_test[cont[0]] = cont[0] + '_' +cont[5]
            if len(cont[12]) > 0:
                renameSample_pfs[cont[0]] = cont[0] + '_' +cont[12]
            if cont[5] == 'SD':
                stable_list.append(cont[0] + '_' +cont[5])

training_sample = [renameSample_test[x] for x in renameSample_pfs]

acc_length_kegg = {}
with open('../kegg.length.csv', 'r') as fin:
    for line in fin:
        cont = line.strip().split(',')
        acc_length_kegg[cont[0]] = float(cont[1])*3/1000

acc_name_kegg = {}
with open('../kegg.name.tsv', 'r') as fin:
    for line in fin:
        cont = line.split()
        acc_name_kegg[cont[0]] = cont[1]

acc_length_pfam = {}
acc_name_pfam = {}
with open('../pfam.acclengname.csv', 'r') as fin:
    for line in fin:
        cont = line.strip().split(',')
        acc_length_pfam[cont[0]] = float(cont[1])/1000
        acc_name_pfam[cont[0]] = cont[2]

def loadPFAM(fname, ff_fast=False):
    df_pfam = pd.DataFrame(pd.read_csv(fname, sep='\t', index_col=0))
    if not ff_fast:
        for col in df_pfam.columns:
            df_pfam[col] = df_pfam[col]*1e6 / df_pfam[col].sum()
        for row in df_pfam.index:
            df_pfam.loc[row] = df_pfam.loc[row] / acc_length_pfam[row]
    df_pfam.sort_index(axis=1, inplace=True)
    return df_pfam

def loadKEGG(fname, ff_fast=False):
    df_kegg  = pd.DataFrame(pd.read_csv(fname, sep='\t', index_col=0))
    if not ff_fast:
        for col in df_kegg.columns:
            df_kegg[col] = df_kegg[col]*1e6 / df_kegg[col].sum()
        for row in df_kegg.index:
            df_kegg.loc[row] = df_kegg.loc[row] / acc_length_kegg[row]
    df_kegg.sort_index(axis=1, inplace=True)
    return df_kegg

def loadTAX(fname):
    df_tax  = pd.DataFrame(pd.read_csv(fname, sep=',', index_col=0))
    df_tax.sort_index(axis=1, inplace=True)
    return df_tax

def loadKMC(fname):
    df_kmc  = pd.DataFrame(pd.read_csv(fname, sep='\t', index_col=0))
    df_kmc.sort_index(axis=1, inplace=True)
    return df_kmc

map_data = {}
if pf == 'pfam':
    df_pfam = loadPFAM(f"{base_dir}/pfam.count.tsv", ff_skipNormalize)
    df_pfam.rename(columns=renameSample_test,inplace=True)
    map_data[pf] = df_pfam
elif pf == 'kegg':
    df_kegg = loadKEGG(f"{base_dir}/kegg.count.tsv", ff_skipNormalize)
    df_kegg.rename(columns=renameSample_test,inplace=True)
    map_data[pf] = df_kegg
elif pf == 'tax':
    df_tax = loadTAX(f"{base_dir}/otu.count.csv")
    df_tax.rename(columns=renameSample_test,inplace=True)
    map_data[pf] = df_tax
elif pf == 'k30':
    df_kmc = loadKMC(f"{base_dir}/k30.count.tsv")
    df_kmc.rename(columns=renameSample_test,inplace=True)
    map_data[pf] = df_kmc
elif pf == 'k16':
    df_kmc = loadKMC(f"{base_dir}/k16.count.tsv")
    df_kmc.rename(columns=renameSample_test,inplace=True)
    map_data[pf] = df_kmc
elif pf == 'k60':
    df_kmc = loadKMC(f"{base_dir}/k60.count.tsv")
    df_kmc.rename(columns=renameSample_test,inplace=True)
    map_data[pf] = df_kmc
elif pf == "com":
    df_tax = loadTAX(f"{base_dir}/otu.count.csv")
    df_tax.rename(columns=renameSample_test,inplace=True)
    df_kmc = loadKMC(f"{base_dir}/k30.count.tsv")
    df_kmc.rename(columns=renameSample_test,inplace=True)
    map_data[pf] = pd.concat([df_tax, df_kmc])
    print(df_tax.shape, df_kmc.shape)


print(f"{pf}, {map_data[pf].shape}")
"""== Test bestK ============================================================"""
def pickFeature(df, select):
    if len(set(df.columns) & set(select)) == 0:
        print("None of the selected feature was found in data")
        return df
    return df.T.reindex(select).fillna(0).T.copy()

def simpleName(name):
    if type(name) == int:
        return str(name)
    if name[:2] == 'PF' and name[2:].isdigit():
        return name[2:]
    elif name[0] == 'K' and name[1:].isdigit():
        return name[1:]
    cont = name.split(';')
    return cont[-1][3:]

def my_split1(df_join, Y, ds):
    ds_label = getDS(df_join.index)
    if ds == 'dds3':
        df_train = df_join[ds_label!=2]
        y_train = Y[ds_label!=2]
        df_join = df_join[ds_label==2]
        Y = Y[ds_label==2]
        N = len(ds_label)
    elif ds in ['ds1', 'ds2', 'ds3']:
        df_join = df_join[ds_label==(int(ds[-1])-1)]
        Y = Y[ds_label==(int(ds[-1])-1)]
        df_train = pd.DataFrame()
        y_train = np.array([])
        N = len(Y)
    elif ds == 'ds12':
        df_join = df_join[ds_label!=2]
        Y = Y[ds_label!=2]
        df_train = pd.DataFrame()
        y_train = np.array([])
        N = len(Y)
    else:
        df_train = pd.DataFrame()
        y_train = np.array([])
        N = len(Y)
    rs_sp = np.random.randint(9999999)
    kf = StratifiedShuffleSplit(n_splits=1, test_size=int(np.ceil(0.3*N)), random_state=rs_sp)
    for T_index, t_index in kf.split(df_join, Y):
        df_tmp , y_tmp   = df_join.iloc[T_index], Y[T_index]
        df_test,  y_test = df_join.iloc[t_index], Y[t_index]
    df_train = pd.concat([df_train, df_tmp])
    y_train = np.concatenate((y_train, y_tmp))
    return (rs_sp, df_train, y_train, df_test, y_test)

def my_training_testing_split(df_join, ds):
    if ds in ['mix', 'dds3', 'ds1', 'ds2', 'ds3', 'ds12']:
        Y = getGT(df_join.index, sp2)
        return my_split1(df_join, Y, ds)
    if ds == 'pfs3':
        df_all = df_join.loc[training_sample]
        pfs_index = [renameSample_pfs[x[:-3]] for x in df_all.index]
        Y = getGT(pfs_index, 'pfs')
        rs_sp = np.random.randint(9999999)
        kf = StratifiedShuffleSplit(n_splits=1, test_size=0.3, random_state=rs_sp)
        for T_index, t_index in kf.split(df_all, Y):
            df_train, y_train = df_all.iloc[T_index], Y[T_index]
            df_test,  y_test  = df_all.iloc[t_index], Y[t_index]
        return (rs_sp, df_train, y_train, df_test, y_test)
    ds_label = getDS(df_join.index)
    df_train = df_join.loc[training_sample]
    if sp1 == 'pfs':
        pfs_index = [renameSample_pfs[x[:-3]] for x in df_train.index]
        y_train  = getGT(pfs_index, sp1)
    else:
        y_train  = getGT(df_train.index, sp1)

    if ds == 'pfs1':
        df_test  = df_join[ds_label!=2]
    elif ds == 'pfs2':
        df_test  = df_join.loc[set(df_join.index[ds_label!=2])-set(training_sample)]
    y_test   = getGT(df_test.index, sp2)
    return (-1, df_train, y_train, df_test, y_test)

def my_training(XX, yy, X, y):
    if code_stats[0] == '-':
        max_mcc = -2
        mark = -1
        for i in range(0,1000):
            try:
                clf = MLPClassifier(hidden_layer_sizes=(100,100), max_iter=400, random_state=i, early_stopping=True)
                clf.fit(XX, yy)
                p_test = clf.predict(X)
                mcc = f1_score(y, p_test)
                if mcc > max_mcc:
                    max_mcc = mcc
                    mark = i
            except:
                print("bad random_state: {}".format(i))
        return mark
    elif code_stats[0] == '@':
        return np.random.randint(1000)
    else:
        return int(code_stats)

def crossValidateSingle():
    df_join = pickFeature(map_data[pf].T.copy(), selected)
    fout = open(f"{work_space}{ds}_({code_stats})_{pf}_train({sp1})_test({sp2})_{ss_name}_prob.csv", 'w')
    fout.close()
    rs_states  = np.zeros((top, 5))
    for index in range(top):
        (rs_sp, df_train, y_train, df_test, y_test) = my_training_testing_split(df_join, ds)
        print(df_train.shape, df_test.shape, set(y_train), set(y_test))
        X_train = normalize_feature(df_train)
        X_test  = normalize_test(df_test, df_train)
        rs_model = my_training(X_train, y_train, X_test, y_test)
        print(rs_model)
        clf = MLPClassifier(hidden_layer_sizes=(100,100), random_state=rs_model, max_iter=400, early_stopping=True)
        clf.fit(X_train, y_train)
        # random states
        p_test = clf.predict(X_test)
        rs_states[index,0] = rs_sp
        rs_states[index,1] = rs_model
        rs_states[index,2] = recall_score(    y_test, p_test)
        rs_states[index,3] = precision_score( y_test, p_test, zero_division=1)
        rs_states[index,4] = f1_score(        y_test, p_test)
        # Probability
        pred = clf.predict_proba(X_test)
        l1 = ""
        l2 = ""
        for it in sorted(zip(pred[:,1], y_test)):
            l1 += str(it[0]) + ','
            l2 += str(int(it[1])) + ','

        fout = open(f"{work_space}{ds}_({code_stats})_{pf}_train({sp1})_test({sp2})_{ss_name}_prob.csv", 'a')
        fout.write(l1[:-1]+'\n')
        fout.write(l2[:-1]+'\n')
        fout.close()

    df_scores = pd.DataFrame(rs_states)
    df_scores.columns = ["rs_sp", "rs_model", "recall", "precision", "f1_score"]
    df_scores.to_csv(f"{work_space}{ds}_({code_stats})_{pf}_train({sp1})_test({sp2})_{ss_name}.csv")

crossValidateSingle()
