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

loc = {'S' : '', 'E' : '', 'R' : ''}

def setLoc(Ptype):

    if Ptype == "OpenMP":
        loc['S'] = "upper left"
        loc['E'] = "upper right"
        loc['R'] = "upper right"
        
    elif Ptype == "MPI":
        loc['S'] = "upper right"
        loc['E'] = "upper right"
        loc['R'] = "upper right"

def PlotOverProcessors(dict_, name, xlabel, ylabel, filename, loc, all=False, Dir=None, label= None, PlotOpt=False):
    fig = plt.figure(name)
    if not all:
        fig.clear()
    ax = fig.gca()
    fig.subplots_adjust(bottom=0.2, left=0.2)

    x = []
    y = []
    for i in sorted(dict_.items()):
        x.append(i[0])
        y.append(i[1])

    ax.set_xlabel(xlabel, fontsize=fs)
    ax.set_ylabel(ylabel, rotation=0, labelpad=30, fontsize=fs)
    
    ax.set_xlim(xmin=1, xmax = max(sorted(dict_.keys())))
    
    Max = max(sorted(dict_.values()))
    if PlotOpt and filename == "Efficiency":
        Max += 0.1 * Max
    ax.set_ylim(ymin=0, ymax = Max)
    
    ax.grid()
    ax.xaxis.set_major_locator(plt.MultipleLocator(1))
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)

    if not all:
        if PlotOpt and filename == "Efficiency":
            ax.plot(x,y, label=label)
            ax.plot(x, [1 for i in x], label = 'optimum')
            ax.legend(loc='center right', framealpha=1.0, fontsize=fs)
        else:
            ax.plot(x,y)
        plt.savefig(os.getcwd() + "/" + filename + "_" + name + ".eps")
    else:
        ax.plot(x,y, label=label)
        if PlotOpt and filename == "Efficiency" and label == "fine":
            ax.plot(x, [1 for i in x], label = 'optimum')
            loc = 'center right'
        ax.legend(loc = loc, framealpha=1.0, fontsize=fs)
        plt.savefig(Dir + "/" + filename + ".eps")

def PlotAllInOne(Dir, l, PlotOpt=False):
    PlotOverProcessors(speedup, "all1", "number of threads", r"$S_p = \frac{T_1}{T_p}$", "SpeedUp", loc['S'], True, Dir, l)
    PlotOverProcessors(eff, "all2", "number of threads", r"$E_p = \frac{S_p}{p}$", "Efficiency", loc['E'],True, Dir, l, PlotOpt=PlotOpt)
    PlotOverProcessors(averageTimes, "all3", "number of threads", r"$T [s]$", "Runtime", loc['R'],True, Dir, l)
    
def main():
    directory = os.getcwd()

    setLoc(sys.argv[1])

    PlotOpt = False
    if sys.argv[1] == "MPI":
        PlotOpt = True
    
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
        PlotOverProcessors(speedup, Dir, "number of threads", r"$S_p = \frac{T_1}{T_p}$", "SpeedUp", loc['S'])
        PlotOverProcessors(eff, Dir, "number of threads", r"$E_p = \frac{S_p}{p}$", "Efficiency", loc['E'], label = Dir, PlotOpt=PlotOpt)
        PlotOverProcessors(averageTimes, Dir, "number of threads", r"$T [s]$", "Runtime", loc['R'])
        
        PlotAllInOne(directory, Dir, PlotOpt=PlotOpt)

if __name__ == "__main__":
    main()