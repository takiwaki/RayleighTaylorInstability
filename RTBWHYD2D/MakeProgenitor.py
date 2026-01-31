#!/usr/bin/python

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


Msun   = 1.989e33

input="s12.0_presn"

def Main():
    dfinp = ReadData(input)
    dfori = ReduceElement(dfinp)
    SaveData(dfori,"s12.0.txt")
    dfmod = ModifyData(dfori)
    SaveData(dfmod,"s12.0mod.txt")
    dfana = AnalyticModel()
    SaveData(dfana,"analytic.txt")
    CompareData(dfori,dfana)
    
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
    dfout = df.copy()
    mask = dfout["radius"] > 2.0e9
    dfout.loc[mask,"density"] /= 10.0
    dfout.loc[mask,"pressure"] /= 10.0
    return dfout

def AnalyticModel():
    r_max = 1.0e13 # stellar radius [cm]
    r_min = 0.0e0
    L = (r_max-r_min)
    nr=500
    drmin= L /nr /100
    
    r = 1.01
    eps = 1.0e-8
    max_iter = 30
    for _ in range(max_iter):
        # fn = (x1max-x1min) - dx1min*(r^n - 1)/(r - 1)
        fn = L - drmin * (r**nr - 1.0) / (r - 1.0)

        # dfndx1r = -nbl*dx1min*r^(nbl-1)/(r-1) + dx1min*(r^n - 1)/(r-1)^2
        df = (-nr * drmin * r**(nr - 1) / (r - 1.0)
              + drmin * (r**nr - 1.0) / (r - 1.0)**2)

        deltr = -fn / df
        err = abs(deltr / r)
        r = r + deltr

        if err < eps:
            break

    if err > eps:
        raise RuntimeError(f"coordinate error: {err} (did not converge)")

    dr = drmin * r**np.arange(nr)        # dx_k
    r_e = np.empty(nr + 1, dtype=np.float64)
    r_e[0] = r_min
    r_e[1:] = r_min + np.cumsum(dr)
    r_c = (r_e[:-1] + r_e[1:])/2

    
    # ===== component definition (0-based) =====
    ncomp = 4
    nFe, nCO, nHe, nH = 0, 1, 2, 3

    pi = np.arccos(-1.0)

    # ===== shell radii =====
    Rshell = np.zeros(ncomp)
    Rshell[nFe] = 2.0e8
    Rshell[nCO] = 1.0e9
    Rshell[nHe] = 1.0e10
    Rshell[nH ] = r_max

    # ===== density slopes =====
    slope = np.zeros(ncomp)
    slope[nFe] = 2.0
    slope[nCO] = 2.0
    slope[nHe] = 2.8
    slope[nH ] = 0.5

    # ===== densities at shell boundaries =====
    rho0 = np.zeros(ncomp)
    rho0[nFe] = 1.5e7
    rho0[nCO] = rho0[nFe] * (Rshell[nCO] / Rshell[nFe]) ** (-slope[nFe])
    rho0[nHe] = 1.0e-1 * rho0[nCO] * (Rshell[nHe] / Rshell[nCO]) ** (-slope[nCO])
    rho0[nH ] = 1.0e-1 * rho0[nHe] * (Rshell[nH ] / Rshell[nHe]) ** (-slope[nHe])

    # ===== arrays =====
    density  = np.empty(nr)
    pressure = np.empty(nr)
    velocity = np.zeros(nr)

    X = np.full((ncomp, nr), 1.0e-10)

    # ===== fill structure =====
    for i in range(nr):
        r = r_c[i]

        if r < Rshell[nFe]:  # Fe core
            density[i]  = rho0[nFe] * (r / Rshell[nFe]) ** (-slope[nFe])
            pressure[i] = 5.0e17 * density[i]
            X[nFe, i] = 1.0

        elif r < Rshell[nCO]:  # CO core
            density[i]  = rho0[nCO] * (r / Rshell[nCO]) ** (-slope[nCO])
            pressure[i] = 1e17 * density[i]
            X[nCO, i] = 1.0

        elif r < Rshell[nHe]:  # He core
            density[i]  = rho0[nHe] * (r / Rshell[nHe]) ** (-slope[nHe])
            pressure[i] = 1.0e16 * density[i]
            X[nHe, i] = 1.0

        else:  # H envelope
            density[i]  = rho0[nH] * (r / Rshell[nH]) ** (-slope[nH])
            pressure[i] = 0.5e16 * density[i]
            X[nH, i] = 1.0

    # ===== cell mass =====
    dV = 4.0 * pi * r_c**2 * dr
    mass = density * dV

    # ===== DataFrame output =====
    dfout = pd.DataFrame({
        "radius": r_c,
        "mass": mass,
        "density": density,
        "pressure": pressure,
        "H":  X[nH],
        "He": X[nHe],
        "CO": X[nCO],
        "Fe": X[nFe],
    })
    Msolar = 1.989e33
    # ===== summary (optional) =====
    Mshell = {
        "Fe": (mass * X[nFe]).sum() / Msolar,
        "CO": (mass * X[nCO]).sum() / Msolar,
        "He": (mass * X[nHe]).sum() / Msolar,
        "H" : (mass * X[nH ]).sum() / Msolar,
    }

    print("=== Progenitor summary [M_sun] ===")
    for k, v in Mshell.items():
        print(f"{k:>2} : {v:.4e}")

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
    fig.savefig("compareDen.png", dpi=150)
    plt.close()
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    x = df1["radius"]
    y = df1["pressure"]*x**3
    ax.plot(x, y)
    x = df2["radius"]
    y = df2["pressure"]*x**3
    ax.plot(x, y)
    
    ax.set_xlabel("radius [cm]")
    ax.set_ylabel(r"$p r^3$ [g]")
    ax.set_xscale("log")
    ax.set_xlim(1.0e8,5.0e11)
    ax.set_yscale("log")
    
    fig.savefig("comparePre.png", dpi=150)
    plt.close()
    
Main()
