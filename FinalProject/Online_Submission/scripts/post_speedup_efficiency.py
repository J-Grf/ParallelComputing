#!/usr/local/bin/python3

# TODO
# - add pgf format
# - plot scattered points and get value at point

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

pointS = 12

loc = {'S' : '', 'E' : '', 'R' : ''}

f = 0.8

def setLoc(Ptype):

    if Ptype == "OpenMP":
        loc['S'] = "upper left"
        loc['E'] = "upper right"
        loc['R'] = "upper right" 
        
    elif Ptype == "MPI":
        loc['S'] = "upper right"
        loc['E'] = "upper right"
        loc['R'] = "upper right"

def PlotOverProcessors(dict_, name, xlabel, ylabel, filename, loc, pointOffSetX, pointOffSetY, all=False, Dir=None, label= None, PlotOpt=False, PType=None):
    fig = plt.figure(name)
    if not all:
        fig.clear()
    ax = fig.gca()
    fig.subplots_adjust(bottom=0.2, left=0.2)

    #print(pointOffSetY)

    x = []
    y = []
    for i in sorted(dict_.items()):
        x.append(i[0])
        y.append(i[1])

    ax.set_xlabel(xlabel, fontsize=fs)
    ax.set_ylabel(ylabel, rotation=0, labelpad=30, fontsize=fs)
    
    x_O = 1
    if x[-1] > 12:
        x_O = 2
    ax.set_xlim(xmin=0, xmax = max(sorted(dict_.keys())) + x_O)
    
    Max = max(sorted(dict_.values()))
    #if PlotOpt and filename == "Efficiency":
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
            ax.scatter(x,y)
            for index in range(len(x)):
                ax.text(x[index] + pointOffSetX[filename], y[index] + pointOffSetY[filename][label], round(y[index],2), size = pointS)
            
            ax.legend(loc='center right', framealpha=1.0, fontsize=fs)
        else:
            ax.plot(x,y)
            ax.scatter(x,y)
            for index in range(len(x)):
                ax.text(x[index] + pointOffSetX[filename], y[index] + pointOffSetY[filename][label], round(y[index],2), size = pointS)
        plt.savefig(os.getcwd() + "/" + filename + "_" + name + ".eps")
        pgf.savePgf(PType + "/" + filename +"_" + name + ".pgf", factor=f)
    else:
        ax.plot(x,y, label=label)
        ax.scatter(x,y)
        for index in range(len(x)):
                ax.text(x[index] + pointOffSetX[filename], y[index] + pointOffSetY[filename]['fine'], round(y[index],2), size = pointS)
        if PlotOpt and filename == "Efficiency" and label == "fine":
            ax.plot(x, [1 for i in x], label = 'optimum')
            loc = 'center right'
        ax.legend(loc = loc, framealpha=1.0, fontsize=fs)
        plt.savefig(Dir + "/" + filename + ".eps")
        pgf.savePgf(PType + "/" + filename + ".pgf", factor=f)

def PlotAllInOne(Dir, l, xLabel, pointOX, pointOY, PlotOpt=False, PType=None):
    PlotOverProcessors(speedup, "all1", xLabel, r"$S_p = \frac{T_1}{T_p}$", "SpeedUp", loc['S'], pointOX, pointOY, all = True, Dir = Dir, label = l, PType=PType)
    PlotOverProcessors(eff, "all2", xLabel, r"$E_p = \frac{S_p}{p}$", "Efficiency", loc['E'], pointOX, pointOY, all = True, Dir = Dir, label = l, PlotOpt=PlotOpt, PType=PType)
    PlotOverProcessors(averageTimes, "all3", xLabel, r"$T [s]$", "Runtime", loc['R'], pointOX, pointOY, all = True, Dir = Dir, label = l, PType=PType)
    
def main():
    directory = os.getcwd()

    setLoc(sys.argv[1])

    PlotOpt = False
    
    pointOffSetX = 0
    pointOffSetY = {'SpeedUp' : {}, 'Efficiency' : {}, 'Runtime' : {}}

    PType = ""

    xLabel = ""
    if sys.argv[1] == "MPI":
        PlotOpt = True
        PType = "MPI"
        pointOffSetX = {'SpeedUp' : 0.25, 'Efficiency' : 0.0, 'Runtime' : -0.5}
        pointOffSetY['SpeedUp'] = {'coarse' : 0.025, 'medium' : 0.04, 'fine' : 0.05}
        pointOffSetY['Efficiency'] = {'coarse' : 0.025, 'medium' : 0.025, 'fine' : 0.025} 
        pointOffSetY['Runtime'] = {'coarse' : -0.7, 'medium' : -10, 'fine' : -50}
        xLabel = "number of cores"
    elif sys.argv[1] == "OpenMP":    
        PType = "OpenMP"
        pointOffSetX = {'SpeedUp' : 0.0, 'Efficiency' : 0.0, 'Runtime' : 0.0}
        pointOffSetY['SpeedUp'] = {'coarse' : 0.025, 'medium' : 0.1, 'fine' : 0.15}
        pointOffSetY['Efficiency'] = {'coarse' : 0.025, 'medium' : 0.025, 'fine' : 0.025} 
        pointOffSetY['Runtime'] = {'coarse' : 0.04, 'medium' : 1.4, 'fine' : 15}
        xLabel = "number of threads"

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
        PlotOverProcessors(speedup, Dir, xLabel, r"$S_p = \frac{T_1}{T_p}$", "SpeedUp", loc['S'], pointOffSetX, pointOffSetY, label = Dir, PType=PType)
        PlotOverProcessors(eff, Dir, xLabel, r"$E_p = \frac{S_p}{p}$", "Efficiency", loc['E'], pointOffSetX, pointOffSetY, label = Dir, PlotOpt=PlotOpt, PType=PType)
        PlotOverProcessors(averageTimes, Dir, xLabel, r"$T [s]$", "Runtime", loc['R'], pointOffSetX, pointOffSetY, label = Dir, PType=PType)
        
        PlotAllInOne(directory, Dir, xLabel, pointOffSetX, pointOffSetY, PlotOpt=PlotOpt, PType=PType)

if __name__ == "__main__":
    main()