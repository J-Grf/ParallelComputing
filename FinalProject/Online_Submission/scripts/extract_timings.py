#!/usr/bin/python3

import os
import re
import sys

reOutput = ""
if sys.argv[1] == "mpi":
    reOutput = r"output_(?:coarse|medium|fine)_(\d+)CPU.*"
elif sys.argv[1] == "openmp":
    reOutput = r"out_(\d+)CPU_.*"

directory = os.fsencode(os.getcwd())
timings = {}
for file in os.listdir(directory):
     filename = os.fsdecode(file)
     match = re.match(reOutput, filename)
     if match: 
        with open(file,'r') as f:
            for line in f:
                if "Elapsed time" in line:
                    match2 = re.match(r"Elapsed time is (?:(\d+).(\d+))", line)
                    print(line)
                    if match.group(1) not in timings:
                        timings[match.group(1)] = []
                    timings[match.group(1)].append(match2.group(1)+","+match2.group(2))
     else:
         continue

with open('timings.txt', 'w') as f:
    for CPU in sorted(timings.keys()):
        f.write(CPU)
        for time in timings[CPU]:
            f.write("\t" + str(time) + "\n")
        f.write("\n")
