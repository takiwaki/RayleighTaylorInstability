#!/usr/bin/python

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


Msun   = 1.989e33

input="s12.0_presn"

def Main():
    dfinp = ReadData(input)
    dfout = ReduceElement(dfinp)
    SaveData(dfout,"s12.0.txt")
    dfmod = ModifyData(dfout)
    SaveData(dfout,"s12.0mod.txt")
    CompareData(dfout,dfmod)
    
def ReadData(file):
    path = file
    columnnames = [ "grid", "mass", "radius", "velocity", "density", "temperature", "pressure","e_i","entropy","omega", "A_bar", "Y_e","stability","NETWORK","neutrons","H1","He3","He4","C12","N14","O16","Ne20","Mg24","Si28","S32","Ar36", "Ca40","Ti44","Cr48","Fe52","Fe54","Ni56","Fe56","Fe" ]
    df = pd.read_csv(path,sep='\\s+',skiprows=2,header=None,names=columnnames)
    df["Fe56"] =(pd.to_numeric(df["Fe56"],errors="coerce").fillna(0.0))
    #print(df)
    return df

def ReduceElement(df):
    H  = df["H1"]
    He = df["He3"] + df["He4"]
    CO = df["C12"] + df["N14"] + df["O16"] + df["Ne20"] + df["Mg24"] + df["Si28"] + df["S32"] + df["Ar36"] + df["Ca40"] 
    Fe = df["Ti44"] + df["Cr48"] + df["Fe52"] + df["Fe54"] + df["Ni56"] + df["Fe56"]

    dfout = pd.DataFrame({
        "radius": df["radius"],
        "mass": df["mass"],
        "density": df["density"],
        "pressure": df["pressure"],
        "H" : H,
        "He": He,
        "CO": CO,
        "Fe": Fe,
        })
    return dfout

def ModifyData(df):
    nmax = len(df)
    denmod = df["density"]
    for n in range(nmax):
        #print (df["radius"][n])
        if(df["radius"][n] > 1.0e10 ):
            denmod[n]= denmod[n]/10.0
            
    dfout = pd.DataFrame({
        "radius": df["radius"],
        "mass": df["mass"],
        "density": denmod,
        "pressure": df["pressure"],
        "H" : df["H"],
        "He": df["He"],
        "CO": df["CO"],
        "Fe": df["Fe"],
        })
    return dfout

def SaveData(df,output):
    nmax = len(df)
    with open(output,"w",encoding="utf-8") as f:
        f.write(f"# {nmax}\n")
        df.to_csv(f,sep=" ",index=False)

def CompareData(df1,df2):
    fnameforfig="sans-serif"
    fsizeforfig=14
    fsizeforlabel=16
    plt.rcParams['font.family'] = fnameforfig
    plt.rcParams['font.size'] = fsizeforfig
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams["xtick.minor.visible"] = True
    plt.rcParams["ytick.minor.visible"] = True
    plt.rcParams['xtick.top'] = True
    plt.rcParams['ytick.right'] = True

    fig = plt.figure()
    ax = fig.add_subplot(111)
    x = df1["radius"]
    y = df1["density"]*x**3
    ax.plot(x, y)
    x = df2["radius"]
    y = df2["density"]*x**3
    ax.plot(x, y)
    
    ax.set_xlabel("radius [cm]")
    ax.set_ylabel(r"$\rho r^3$ [g]")
    ax.set_xscale("log")
    ax.set_xlim(1.0e8,5.0e11)
    ax.set_yscale("log")
    fig.savefig("compare.png", dpi=150)

Main()
