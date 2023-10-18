# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 16:11:30 2023

@author: Gulben AVSAR
"""
import anndata as ad
import pandas as pd
import scanpy as sc
# from scripts import paper_label_SEQdata as lsc
from scripts import useful_funcs as uf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcdefaults()
sc.set_figure_params(scanpy=True, dpi_save=400)


def _itrc_mtx(df):
    #one-by-one, seperately
        c = list(dict.fromkeys([j.split('|')[0] for j in df.columns]))
        
        mtx = pd.DataFrame(index=c, columns=c)
        
        for k in df.columns:
            x = np.count_nonzero(~np.isnan(df[k]))
            
            mtx.loc[k.split('|')[0], k.split('|')[1]] = x
        
        for cs in mtx.columns:
            mtx[cs] = mtx[cs].astype(float)
        return(mtx)

def _pl_heatmap_indv(data, sName):
    # heatmap with the number of interactions
    # cells at y-axis are senders, cells at x-axis are receivers
    saveName = sName
    if len(saveName.split('_'))>1: sName=saveName.split('_')[1]
    df = data.copy()
    df = df[df.columns[::-1]]
    norm=plt.Normalize(df.min().min(),df.max().max())
    
    cmap = mpl.colors.LinearSegmentedColormap.from_list("", ['#fff6e2','#c8875e','#5e4c6c'])    
    # xticsk only at the top
    plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = True
    plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = False
    
    figSize=(10,10)
    if sName in ['scC','scA','scB']: figSize=(16,16)
    fig = plt.figure(figsize=figSize)
    plt.imshow(df, cmap=cmap, norm=norm)
    im_ratio = df.shape[1]/df.shape[0] # Calculate (width_of_image/height_of_image)
    # add the numbers on each cell of map
    if sName in ['scC','scA','scB']:
        for (i, j), z in np.ndenumerate(df):
            if sName in ['scA','scB']: fontSize=26
            else: fontSize=22
            plt.text(j, i, '{}'.format(int(z)), ha='center', va='center', fontsize=fontSize)
    else:
        for (i, j), z in np.ndenumerate(df):
            plt.text(j, i, '{}'.format(int(z)), ha='center', va='center', fontsize=20)
    # fsize = 20
    fsize = 26
    if sName not in ['scC','scA','scB']: fsize=24
    plt.xticks(range(len(df.columns)),df.columns.to_list(), fontsize=fsize, rotation=90)
    plt.yticks(range(len(df.columns)),df.index.to_list(), fontsize=fsize)
    cb = plt.colorbar(fraction=0.047*im_ratio)
    # labelsize=16
    cb.ax.tick_params(labelsize=22)
    plt.grid(False)
    fig.tight_layout()
    plt.show()
    fig.savefig('./figures/CCI_totints_{}.png'.format(saveName), dpi=400, bbox_inches='tight')
    return()


def _tot_itrc_mtx(df):
    #the total numb of interactions
        c = list(dict.fromkeys([j.split('|')[0] for j in df.columns]))
        
        mtx = pd.DataFrame(index=c, columns=c)
        tmtx = pd.DataFrame(index=c, columns=c)
        for k in df.columns:
            x = np.count_nonzero(~np.isnan(df[k]))
            
            mtx.loc[k.split('|')[0], k.split('|')[1]] = x
        
        for k in df.columns:
            if k.split('|')[0] != k.split('|')[1]:
                x = mtx.loc[k.split('|')[0], k.split('|')[1]] + mtx.loc[k.split('|')[1], k.split('|')[0]]
                tmtx.loc[k.split('|')[0], k.split('|')[1]] = x
                tmtx.loc[k.split('|')[1], k.split('|')[0]] = x
            else:
                tmtx.loc[k.split('|')[0], k.split('|')[1]] = mtx.loc[k.split('|')[0], k.split('|')[1]]
        
        for cs in tmtx.columns:
            tmtx[cs] = tmtx[cs].astype(float)
        
        return(tmtx)


def fnd_topN(pairsDict, dNames, topn):
    pairs = pairsDict.copy()
    top_n={} #top LRpairs that are most seen in whole dataset
    tot = pd.DataFrame(columns=['interacting_pair','counts'])
    for i in dNames:
        idx = [pairs[i].loc[pairs[i].interacting_pair == x].index[0]
               for x in pairs[i].interacting_pair if 'integrin' not in x.split('_')]
        df = pairs[i].loc[idx,:]
        counts = df[df.columns[12:].to_list()].dropna(axis=0,how='all').count(axis=1)
        idx_top = counts.sort_values(ascending=False)[:topn].index
        top_n[i] = pairs[i].loc[idx_top,['interacting_pair']+pairs[i].columns[12:].to_list()]
        p_counts = pd.DataFrame(df['interacting_pair'])
        p_counts['counts'] = counts
        tot = pd.concat([tot,p_counts])
    return(tot, top_n)


def pl_heatmap_blue(data, savetitle):
    fdf = data.copy()
    norm=plt.Normalize(fdf.min().min(),fdf.max().max())
    # https://www.color-hex.com/color-palette/1018449
    cmap = mpl.colors.LinearSegmentedColormap.from_list("", ['#faf9f6','#64c6e3','#038d98'])
    fig = plt.figure(figsize=(12,8))
    plt.imshow(fdf, cmap=cmap, norm=norm)
    # if sName not in ['GSE155698','seqA','seqB']:
    #     for (i, j), z in np.ndenumerate(fdf):
    #         plt.text(j, i, '{}'.format(int(z)), ha='center', va='center', fontsize=22)
    fsize = 10
    plt.xticks(range(len(fdf.columns)),fdf.columns.to_list(), fontsize=fsize, rotation=90)
    plt.yticks(range(len(fdf.index)),fdf.index.to_list(), fontsize=fsize)
    cb = plt.colorbar(cax = fig.add_axes([1.03, 0.22, 0.03, 0.44]))
    cb.ax.tick_params(labelsize=10)
    fig.tight_layout()
    fig.savefig('./figures/{}.png'.format(savetitle), dpi=400, bbox_inches='tight')
    plt.show()


















