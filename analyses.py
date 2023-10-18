# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 16:02:39 2023

@author: Gulben AVSAR
"""
import matplotlib as mpl
import pandas as pd
from scripts import pfuncs as pf
import matplotlib.pyplot as plt
mpl.rcdefaults()


#-----------------------------------------------------------------------------------------------------
# load significant pairs in all datasets
dNames = ['stA','stA2','stA3','stB','stB2','stB3']
pairs = {}
for i in dNames:
    df = pd.read_csv('./cellphonedb_results/{}/significant_means.txt'.format(i), sep='\t', index_col=0)
    df1 = df.iloc[:,12:].dropna(how='all')
    df = df.iloc[df1.index,:12].join(df1)
    pairs[i] = df.copy()
del df,df1


#-----------------------------------------------------------------------------------------------------
# CCI Heatmaps
for i in ['stA', 'stA2', 'stA3', 'stB', 'stB2', 'stB3']:
    df = pairs[i].copy()
    tmtx = pf._tot_itrc_mtx(df[df.columns[12:]])
    pf._pl_heatmap_indv(tmtx, i)


#-----------------------------------------------------------------------------------------------------
# get top most common 15 pairs after excluding EMC related partners
tot, top15 = pf.fnd_topN(pairs, ['stA','stA2','stA3','stB','stB2','stB3'], 15)

c_tot = tot['interacting_pair'].value_counts()
pair_counts = pd.DataFrame(index=c_tot.index, columns=['occurance'])
for x in c_tot.index:
    pair_counts.loc[x,'occurance'] = tot.loc[tot.interacting_pair == x,:].counts.values.sum()
pair_counts=pair_counts.sort_values(by=['occurance'], ascending=False)
del c_tot, x, tot

pair_counts.index[:15]

speclr = pair_counts.iloc[:15,:].index.tolist()

fdf = pd.DataFrame(index=speclr) # dataframe of the given L/R (LGALS9 here)
for i in ['stA','stA2','stA3','stB','stB2','stB3']:
    idx = [pairs[i][pairs[i].interacting_pair == j].index[0] for j in speclr if j in pairs[i].interacting_pair.values]
    if len(idx)>0:
        df = pairs[i].loc[idx,['interacting_pair']+pairs[i].columns[12:].to_list()]
        df = df.dropna(axis=1, how='all') #remove the nan columns
        df= df.set_index('interacting_pair')
        if i in ['stA','stB']: i=i+'1'
        cols = ['{}_{}'.format(i,c) for c in df.columns]
        df.columns = cols
        fdf = fdf.join(df)
fdf = fdf.reindex(fdf.count(axis=1).sort_values(ascending=False).index)
del idx,cols,df

pf.pl_heatmap_blue(fdf, 'top15_ints_st')


#-----------------------------------------------------------------------------------------------------
################## LGALS9 L-R pairs
# LGALS9 in all datasets
speclra = []
for i in ['stA','stA2','stA3','stB','stB2','stB3']:
    for j in pairs[i].interacting_pair:
        if 'LGALS9' in j:
            speclra.append(j)
del i,j
speclra = list(set(speclra))

df = pd.DataFrame(index=speclra,
                  columns=['stA','stA2','stA3','stB','stB2','stB3'])
for i in pairs.keys():
    for j in speclra:
        if len(pairs[i].loc[pairs[i].interacting_pair == j,pairs[i].columns[12:]].count(axis=1).values)>0:
            c = pairs[i].loc[pairs[i].interacting_pair == j,pairs[i].columns[12:]].count(axis=1).values[0]
            c= c/pairs[i].shape[0] #percentages of the number of occurance in the corresponding dataset
        else:
            c = 0
        df.loc[j,i] = c
for c in df.columns:
    df[c] = df[c].astype(float)
df.columns = ['stA1','stA2','stA3','stB1','stB2','stB3']
df = df[['stA1','stA2','stA3','stB1','stB2','stB3']]
(df != 0).sum(1)
del i,j,c
norm=plt.Normalize(df.min().min(),df.max().max())
cmap = mpl.colors.LinearSegmentedColormap.from_list("", ['#faf9f6','#64c6e3','#038d98'])

fig = plt.figure(figsize=(12,8))
plt.imshow(df, cmap=cmap, norm=norm)
plt.vlines(x=[r+0.5 for r in range(len(df.columns)-1)], ymin=-0.5, ymax=len(speclra)-0.5,
           color='black',lw=0.5)
fsize = 16
plt.xticks(range(len(df.columns)),df.columns.to_list(), fontsize=fsize, rotation=90)
plt.yticks(range(len(df.index)),df.index.to_list(), fontsize=fsize)
cb = plt.colorbar()
cb.ax.tick_params(labelsize=16)
fig.tight_layout()
fig.savefig('./figures/LGALS9_ints_all.png', dpi=400, bbox_inches='tight')
plt.show()







