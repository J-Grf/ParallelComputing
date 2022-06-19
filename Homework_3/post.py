#!/usr/local/bin/python3

import os
import re
import sys
import matplotlib.pyplot as plt

# Extract timings
timings = {}

os.chdir(os.getcwd() + "/" + sys.argv[1])
directory = os.fsencode(os.getcwd())
for file in os.listdir(directory):
     filename = os.fsdecode(file)
     match = re.match("out_(\d+)CPU\.\d+\.txt",filename)
     if match: 
        cpu = int(match.group(1))
        #print("search in file: out_" + str(cpu) +"CPU...")
        with open(filename,'r') as f:
            data = f.readlines()
            for line in data:
                match2 = re.match("Elapsed time is (\d+\.\d+)", line)
                if match2:
                    #print(match2.group(1))
                    timings[cpu] = match2.group(1) 
     else:
         continue

print(timings)

# Plot Scaling
x = []
y = []
for i in sorted(timings.items()):
    x.append(int(i[0]))
    y.append(float(i[1]))
print(x)
print(y)
plt.plot(x,y)
plt.savefig(os.getcwd() + "/Scaling.eps")