#!/usr/bin/python

import sys
import glob
import re
import numpy as np
import matplotlib.pyplot as plt

dir_path = "./output/"
outputdatapath="./figures/"

def Main():
    global dir_path
    files=GetFileList()
    #files = [dir_path + "den00050.dat"]

    is_initial = True
    for file in files:
        print(file)
        time,x,y,rho = ReadData(file)
        
        match = re.search(r'\d+', file)
        fileindex = match.group()

        Plot2DData(fileindex,time,x,y,rho)

def mkdir(path):
    import os
    if not os.path.isdir(path):
        os.makedirs(path)

def GetFileList():
    global dir_path
    filenames= dir_path+"den*.dat"
    files = glob.glob(filenames)
    files = sorted(files)
    return files

def ReadData(file):
    input = file
    try:
        results = np.genfromtxt(input,skip_header=3,delimiter='   ') # Read numbers
#        print(results)
        inputf= open(input, 'r')
        header= inputf.readline() # Read the fisrt line 
        item= header.split()
        t = float(item[1])
        print("T="+str(t))
        header= inputf.readline() # Read the second line
        item= header.split()
#        print(item)
        nx = int(item[1])
        ny = int(item[2])
        print("nx,ny",nx,ny)
        inputf.close()

    except IOError:
        print(" cannot open " + input )
        sys.exit()

    x, y, rho =  np.split(results,3,1)
    x   =   x.reshape(ny,nx)
    y   =   y.reshape(ny,nx)
    rho = rho.reshape(ny,nx)
    return t, x, y, rho

from matplotlib import ticker, cm, colors
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

def Plot2DData(num,time,x,y,rho):
  #######################################
  # Space Time diagram
  #######################################
  global outputdatapath
  timetxt=r"$T=$"+str(time)

  outputfile=outputdatapath+ "den"+ num +".png"
  fig,ax = plt.subplots(figsize=(3.0,4.8))
  ax.set_title(r"$\rho$")
  ax.set_aspect(1.2)
  im=ax.pcolormesh( y,x,rho,cmap="viridis_r")
  ax.set_xlabel(r"$x$", fontsize=fsizeforlabel)
  ax.set_ylabel(r"$z$", fontsize=fsizeforlabel)
  ax.text(-1.0, 1.0,timetxt,transform=ax.transAxes)
  fig.colorbar(im)
  fig.tight_layout()
  print("output"+outputfile)
  fig.savefig(outputfile)
  plt.close()

#######################################
# Execute
#######################################
if __name__ == "__main__":
    Main()
