#!/usr/local/bin/python3

from calendar import firstweekday
import os
import re
import sys
import matplotlib.pyplot as plt
import statistics as st
import pgf

# Extract timings
timings = {}
averageTimes = {}
speedup = {}
eff = {}
All_Sp = {}
All_Eff = {}

meshes = ["coarse", "medium", "fine"]

fs = 15

# Plot Scaling
def PlotScaling():
    x = []
    y = []
    for i in sorted(timings.items()):
        x.append(i[0])
        y.append(st.mean([float(k) for k in i[1]]))
    
    print(x)
    print(y)
    plt.xlabel("number of processors ")
    plt.ylabel("time [s]")
    plt.plot(x,y)
    plt.savefig(os.getcwd() + "/Scaling.eps")

def PlotSpeedUp(name, all=False, Dir=None, label=None):
    fig = plt.figure(name)
    if not all:
        fig.clear()
    ax = fig.gca()
    fig.subplots_adjust(bottom=0.2, left=0.2)

    p = []
    Sp = []
    for i in sorted(speedup.items()):
        p.append(i[0])
        Sp.append(i[1])

    ax.set_xlabel("number of threads", fontsize=fs)
    ax.set_ylabel(r"$S_p = \frac{T_1}{T_p}$", rotation=0, labelpad=30, fontsize=fs)
    ax.set_xlim(xmin=1, xmax = max(sorted(speedup.keys())))
    ax.set_ylim(ymin=0, ymax = max(sorted(speedup.values())))
    ax.grid()
    ax.xaxis.set_major_locator(plt.MultipleLocator(1))
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    if not all:
        ax.plot(p,Sp)
        plt.savefig(os.getcwd() + "/Speedup_" + name + ".eps")
    else:
        ax.plot(p,Sp, label=label)
        ax.legend(loc="upper left", framealpha=1.0, fontsize=fs)
        plt.savefig(Dir + "/SpeedUp.eps")

def PlotEfficiency(name, all=False, Dir=None, label=None):
    fig = plt.figure(name)
    if not all:
        fig.clear()
    ax = fig.gca()
    fig.subplots_adjust(bottom=0.2, left=0.2)

    p = []
    E = []
    for i in sorted(eff.items()):
        p.append(i[0])
        E.append(i[1])

    ax.set_xlabel("number of threads", fontsize=fs)
    ax.set_ylabel(r"$E_p = \frac{S_p}{p}$", rotation=0, labelpad=30, fontsize=fs)
    ax.set_xlim(xmin=1, xmax = max(sorted(eff.keys())))
    ax.set_ylim(ymin=0, ymax = max(sorted(eff.values())))
    ax.grid()
    ax.xaxis.set_major_locator(plt.MultipleLocator(1)) 
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    if not all:
        ax.plot(p,E) 
        plt.savefig(os.getcwd() + "/Efficiency_" + name + ".eps")
    else:
        ax.plot(p,E, label=label) 
        ax.legend(loc="upper right", framealpha=1.0, fontsize=fs)
        plt.savefig(Dir + "/Efficiency.eps")

def PlotAllInOne(Dir, l):
    PlotEfficiency("all1", True, Dir, l)
    PlotSpeedUp("all2",True, Dir, l)
   
def main():
    directory = os.getcwd()
    
    for Dir in meshes:
        print("plotting " + Dir +" ... ")
        os.chdir(directory + "/" + Dir)

        with open("timings.txt", 'r') as f:
            cpu = 0
            for line in f:
                match1 = re.match(r'(\d+)\t(\d+),(\d+)\n', line)
                match2 = re.match(r'\t(\d+),(\d+)\n', line)

                if match1:
                    cpu = int(match1.group(1))
                    if match1.group(1) not in timings.keys():
                        timings[cpu] = []
                    timings[cpu].append(float(match1.group(2) + "." + match1.group(3)))
                elif match2:
                    timings[cpu].append(float(match2.group(1) + "." + match2.group(2)))

        print(timings)

        #calc average times
        for i in sorted(timings.items()):
            averageTimes[i[0]] = st.mean(i[1])

        print(averageTimes)

        #calc speedup and efficiency
        T1 = averageTimes[1]
        for i in averageTimes.items():
            speedup[i[0]] = T1 / i[1]
            eff[i[0]] = speedup[i[0]] / i[0]
        
        print(speedup)
        print(eff)
        PlotSpeedUp(Dir)
        PlotEfficiency(Dir)
        PlotAllInOne(directory, Dir)

if __name__ == "__main__":
    main()