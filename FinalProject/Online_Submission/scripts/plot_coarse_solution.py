#!/usr/local/bin/python3

import pandas as pd
import matplotlib.pyplot as plt
import os
import math
import pgf

class TemperatureDist:
    R1 = 0.01
    R2 = 0.1
    T_D = 500
    Q = 100
    T_0 = T_D
    k = 1

    def __init__(self, path) -> None:
        self.df = pd.read_csv(path, sep=',', usecols=['Temperature','Points:0','Points:1'], decimal='.') #converters={"Temperature": lambda x: "NaN" if x == "nan" else x})
        self.df.dropna(subset=["Temperature"], inplace=True)
        for i in ['Temperature','Points:0','Points:1']:
            self.df[i].apply(lambda x: float(x))
        print(self.df)

    @classmethod 
    def AnalyticalSolution(self, r):
        r = abs(r)
        if r < (self.R1 - 1e-10):
            return self.T_0 - self.Q/(2.0 * math.pi * self.k) * (0.5 * ((r ** 2) / (self.R1 ** 2) - 1) + math.log(self.R1/self.R2))
        else:
            return self.T_0 - self.Q/(2.0 * math.pi * self.k) * math.log(r / self.R2)

    @staticmethod
    def PlotTmp(TmpArray, R2):
        fs = 15
        fig = plt.figure("2D Unsteady")
        ax = fig.gca()

        names = ['coarse', 'fine', 'finest']

        r = []
        for Tmp in TmpArray:
            if len(r) == 0:
                for i in range(0, len(Tmp.df['Points:0'])):
                    fac = 0
                    if Tmp.df.iloc[i]['Points:0'] < 0 and Tmp.df.iloc[i]['Points:1'] < 0:
                        fac = -1
                    else:
                        fac = 1
                    r.append(fac * math.sqrt(Tmp.df.iloc[i]['Points:0']**2 + Tmp.df.iloc[i]['Points:1']**2 ))
                    
            #for i in Tmp.df["Temperature"]:
            ax.plot(r, Tmp.df['Temperature'], label=names[TmpArray.index(Tmp)])

        TmpA = []
        rA = []
        n = 10000
        dr = 2 * R2 / n
        for i in range(int(-n/2),int(n/2)):
            rA.append(i * dr)
            TmpA.append(TemperatureDist.AnalyticalSolution(i * dr))
        ax.plot(rA, TmpA, label = "analytic")

        ax.set_xlabel("r in [m]", fontsize=fs)
        ax.set_ylabel("T in [K]", fontsize=fs, rotation=0, labelpad = 20)
        ax.set_xlim(xmin = -0.1, xmax = 0.1)
        leg = ax.legend(fontsize=fs,framealpha=1.0)
        leg.get_frame().set_edgecolor('k')
        ax.grid()
        
        plt.savefig(os.getcwd() + "/TemperatureDist.eps")
        pgf.savePgf("serial/TemperatureDist.pgf", factor=0.8)

directory = os.fsencode(os.getcwd())
Tmp = []    
for file in os.listdir(directory):
     filename = os.fsdecode(file)
     if filename.endswith(".csv"): 
        Tmp.append(TemperatureDist(filename))
     else:
         continue

TemperatureDist.PlotTmp(Tmp, TemperatureDist.R2)
