import os
import sys
# os.environ['OPENBLAS_NUM_THREADS'] = '4'
# sys.path.append("/home/ben0522/toolkit/pyLib")

import pandas as pd
import numpy as np
import pickle
from sklearn.ensemble import RandomForestClassifier

from sklearn.model_selection import StratifiedShuffleSplit

from sklearn.metrics import recall_score
from sklearn.metrics import precision_score
from sklearn.metrics import f1_score

import warnings
warnings.filterwarnings('ignore')

helper =  "call by\n    ~ <profile> <sp_code>"
helper += "\t\t<profile>     - code for profile, must from {pfam, tax, k30}"
helper += "\t\t<sp_code>     - code for training-testing split, must from {mix, dds3, ds3, ds12}"
helper += "\t\t                 mix,  [1,2,3|1,2,3], mix all samples from all three datasets"
helper += "\t\t                 dds3, [1,2,3|3], using all DS1 and DS2 in addition to part of DS3, same number of training sample as [1,2,3|1,2,3]"
helper += "\t\t                 ds3,  [3|3]"
helper += "\t\t                 ds12, [1,2|1,2]"
helper += "\t\t                 default: mix"

if len(sys.argv) == 2:
    pf = sys.argv[1]
    ds = 'mix'
elif len(sys.argv) == 3:
    pf = sys.argv[1]
    ds = sys.argv[2]
else:
    print(helper)
    exit()

if ds not in ['mix', 'dds3', 'ds3', 'ds12']:
    print("you must choose <sp_code> from {mix, dds3, ds3, ds12}!")
    print(helper)
    exit()

if pf not in ['pfam', 'k30', 'tax']:
    print("you must choose <profile> from {pfam, tax, k30}!")
    exit()

top = 100

"""== Parameters ============================================================"""
ff_skipNormalize = False
map_class = {"PR":1, "CR":1, "PD":0, "DD":0, "SD":1}
base_dir   = f"../data"
work_space = f"{base_dir}/imp"
os.system("mkdir -p " + work_space) # create work_space if not exist
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

def getGT(labels):
    gt = np.zeros(len(labels))
    for (i, v) in enumerate(labels):
        gt[i] = map_class[v[-2:]]
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
with open(f"{base_dir}/metadata.csv") as fin:
    for line in fin:
        if not line[0] == '#':
            cont = line.strip().split(',')
            renameSample_test[cont[0]] = cont[0] + '_' +cont[5]
            if len(cont[12]) > 0:
                renameSample_pfs[cont[0]] = cont[0] + '_' +cont[12]
            if cont[5] == 'SD':
                stable_list.append(cont[0] + '_' +cont[5])

training_sample = [renameSample_test[x] for x in renameSample_pfs]

acc_length_pfam = {}
acc_name_pfam = {}
with open(f"{base_dir}/pfam.acclengname.csv", 'r') as fin:
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

def loadTAX(fname):
    df_tax  = pd.DataFrame(pd.read_csv(fname, sep=',', index_col=0))
    df_tax.sort_index(axis=1, inplace=True)
    return df_tax

def loadKMC(fname):
    df_kmc  = pd.DataFrame(pd.read_csv(fname, sep='\t', index_col=0))
    df_kmc.sort_index(axis=1, inplace=True)
    return df_kmc

if pf == 'pfam':
    df_pfam = loadPFAM(f"{base_dir}/pfam.count.tsv", ff_skipNormalize)
    df_pfam.rename(columns=renameSample_test,inplace=True)
    map_data = df_pfam
elif pf == 'tax':
    df_tax = loadTAX(f"{base_dir}/otu.count.csv")
    df_tax.rename(columns=renameSample_test,inplace=True)
    map_data = df_tax
elif pf == 'k30':
    df_kmc = loadKMC(f"{base_dir}/k30.count.tsv")
    df_kmc.rename(columns=renameSample_test,inplace=True)
    map_data = df_kmc

def my_training(XX, yy, X, y):
    max_mcc = -2
    mark = -1
    for i in range(0,1000):
        try:
            clf = RandomForestClassifier(n_estimators=100, random_state=i, n_jobs=1)
            clf.fit(XX, yy)
            p_test = clf.predict(X)
            mcc = f1_score(y, p_test)
            if mcc > max_mcc:
                max_mcc = mcc
                mark = i
        except:
            print("bad random_state: {}".format(i))
    return mark

def my_training_testing_split(df_join, ds, rs_sp=None):
    Y = getGT(df_join.index)
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
    if not rs_sp:
        rs_sp = np.random.randint(9999999)
    else:
        rs_sp = int(rs_sp)
    kf = StratifiedShuffleSplit(n_splits=1, test_size=int(np.ceil(0.3*N)), random_state=rs_sp)
    for T_index, t_index in kf.split(df_join, Y):
        df_tmp , y_tmp   = df_join.iloc[T_index], Y[T_index]
        df_test,  y_test = df_join.iloc[t_index], Y[t_index]
    df_train = pd.concat([df_train, df_tmp])
    y_train = np.concatenate((y_train, y_tmp))
    return (rs_sp, df_train, y_train, df_test, y_test)

def computeImp():
    df_join = map_data.T.copy()
    feature_names = np.array(list(df_join.columns))
    feature_importances = pd.DataFrame(np.zeros((len(feature_names),top)), index=feature_names)
    for index in range(top):
        rs_sp = np.random.randint(9999999)
        (_, df_train, y_train, df_test, y_test) = my_training_testing_split(df_join, ds, rs_sp)
        print(df_train.shape, df_test.shape, set(y_train), set(y_test))
        X_train = normalize_feature(df_train)
        X_test  = normalize_test(df_test, df_train)
        rs_model = my_training(X_train, y_train, X_test, y_test)
        # benchmark
        clf = RandomForestClassifier(n_estimators=100, random_state=rs_model, n_jobs=16)
        clf.fit(X_train, y_train)

        importances = clf.feature_importances_
        indices = np.argsort(importances)[::-1]
        for i in range(len(importances)):
            feature_importances.loc[feature_names[indices[i]]][index] = importances[indices[i]]
        print(index)

    feature_importances.to_csv(f"{work_space}/RF_{pf}_{ds}_importance.csv", header=False)

####
computeImp()
