#!/usr/bin/python

import matplotlib.pyplot as plt
import csv
import math
import numpy as np
import scipy as sp
from   numpy import array
from   scipy import stats
import os
import urllib2
import sys

######################################################
# Settings #
NUMBER_OF_TESTCASES = 30
NUMBER_OF_TTIS = 20
USE_XKCD = False
BAR_WIDTH = 0.8
SPACE_BETWEEN_GROUPS = BAR_WIDTH * 0.5
COLORSP1 = ["#a3c1fe","#9dbb61","#a5a5a5","#FD0006"]
COLORS = COLORSP1
######################################################

# read csv
res = csv.reader(open('lte.csv'), delimiter=';')
res.next()

if USE_XKCD:
    plt.xkcd();

results = {"sequential":[],"horizon":[],"parLbts":[],"parIBS":[],"comLbts":[],"comIBS":[],"par12LPLbts":[],"par12LPIBS":[],}

def intod(value):
    try:
        return int(value)
    except:
        return 0

def snapr(data):
    return [intod(ele) for ele in data.split(',')]

rowindex = 1
for col in res:
    runtype = ""
    if rowindex <= NUMBER_OF_TESTCASES:
        runtype = "horizon"
    elif rowindex <= 2 * NUMBER_OF_TESTCASES:
        runtype = "sequential"
    elif rowindex <= 3 * NUMBER_OF_TESTCASES:
        runtype = "parLbts"
    elif rowindex <= 4 * NUMBER_OF_TESTCASES:
        runtype = "parIBS"
    elif rowindex <= 5 * NUMBER_OF_TESTCASES:
        runtype = "comLbts"
    elif rowindex <= 6 * NUMBER_OF_TESTCASES:
        runtype = "comIBS"
    elif rowindex <= 7 * NUMBER_OF_TESTCASES:
        runtype = "par12LPLbts"
    elif rowindex <= 8 * NUMBER_OF_TESTCASES:
        runtype = "par12LPIBS"

    results[runtype].append({"runtime":col[1]})

    rowindex += 1

def mapEntries(li, func):
    return [func(elem) for elem in li]

def buildBarPlot(title, yaxis, data, mapFunction = lambda x: x, offsetx = 0, xlim=None, ylim=None,filename=None, drawTicks=True, drawLegend=True, numberOffset=0.5):
    opacity = 0.9
    fig, ax = plt.subplots()

    rects = []
    ind = max([len(elem[1]) for elem in data])
    width = 1.0/len(data)

    inds =  [0 for i in range(len(data)+1)]
    indss = [0 for i in range(len(data))]

    for i in range(len(data)):
        inds[i+1] = len([0 for ele in data[i][1] if len(ele[1]) > 0])
        if i > 0:
            inds[i+1] += inds[i]

    for i in range(ind):
        means = [array(mapEntries(elem[1][i][1], mapFunction)).mean() for elem in data if len(elem[1][i][1]) > 0]
        median = [calcMedian(array(mapEntries(elem[1][i][1], mapFunction))) for elem in data if len(elem[1][i][1]) > 0]
        #stder = [array(mapEntries(elem[1][i][1], mapFunction)).std()  for elem in data if len(elem[1][i][1]) > 0]
        conf = [calcConfInt(array(mapEntries(elem[1][i][1], mapFunction)))  for elem in data if len(elem[1][i][1]) > 0]

        idxs = [elem for elem in range(len(data)) if len(data[elem][1][i][1]) > 0]

        for j in idxs:
            indss[j] += 1

        rect = ax.bar([offsetx + inds[idxs[ele]]*BAR_WIDTH + idxs[ele]*SPACE_BETWEEN_GROUPS + indss[idxs[ele]]*BAR_WIDTH for ele in range(len(means))], means, BAR_WIDTH, alpha=opacity, color=COLORS[i], yerr=conf)
        
        for r in rect:
            ax.text(r.get_x()+r.get_width()/2., numberOffset+r.get_height(), '%.2f'%r.get_height(), ha='center', va='bottom')
        rects.append(rect)

    ax.set_ylabel(yaxis)
    ax.set_title(title)

    if drawTicks:
        ax.set_xticks([offsetx + inds[ele]*BAR_WIDTH + ele*SPACE_BETWEEN_GROUPS + BAR_WIDTH + indss[ele]*BAR_WIDTH/2.0 for ele in range(len(data))])
        labels = ax.set_xticklabels([elem[0] for elem in data])

    for i in ax.get_xticklabels():
        if(i.get_text() == 'Distributed Horizon'):
            i.set_text('DISTRIBUTED HORIZON')
            i.set_variant('small-caps')
            i.set_family('serif')

    lgd = None
    if(drawLegend):
        lgd = ax.legend([plt.Rectangle((0, 0), 1, 1, fc=COLORS[i]) for i in range(len(data[0][1]))], [elem[0] for elem in data[0][1]], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    if(xlim != None):
        plt.xlim(xlim[0], xlim[1])
    if(ylim != None):
        plt.ylim(ylim[0], ylim[1])

    if filename != None:
        plt.savefig(filename, format='pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    else:
        plt.show()

def buildBoxPlot(title, yaxis, data, mapFunction = lambda x: x, offsetx = 0, xlim=None, ylim=None, BAR_WIDTH=BAR_WIDTH, filename=None):
    opacity = 0.9
    fig, ax = plt.subplots()

    rects = []
    ind = max([len(elem[1]) for elem in data])
    width = 1.0/len(data)

    inds =  [0 for i in range(len(data)+1)]
    indss = [0 for i in range(len(data))]

    for i in range(len(data)):
        inds[i+1] = len([0 for ele in data[i][1] if len(ele[1]) > 0])
        if i > 0:
            inds[i+1] += inds[i]

    for i in range(ind):
        nonzeroelements = [mapEntries(elem[1][i][1], mapFunction) for elem in data if len(elem[1][i][1]) > 0]

        idxs = [elem for elem in range(len(data)) if len(data[elem][1][i][1]) > 0]

        for j in idxs:
            indss[j] += 1

        rect = ax.boxplot(nonzeroelements, 
            positions=[offsetx + inds[idxs[ele]]*BAR_WIDTH + idxs[ele]*SPACE_BETWEEN_GROUPS + indss[idxs[ele]]*BAR_WIDTH for ele in range(len(nonzeroelements))], 
            widths=[BAR_WIDTH for ele in range(len(nonzeroelements))])
        
        rects.append(rect)

    ax.set_ylabel(yaxis)
    ax.set_title(title)

    ax.set_xticks([offsetx + inds[ele]*BAR_WIDTH + ele*SPACE_BETWEEN_GROUPS + (indss[ele]+1)*BAR_WIDTH/2.0 for ele in range(len(data))])
    ax.set_xticklabels([elem[0] for elem in data])

    if(xlim != None):
        plt.xlim(xlim[0], xlim[1])
    if(ylim != None):
        plt.ylim(ylim[0], ylim[1])

    if filename != None:
        plt.savefig(filename)
    else:
        plt.show()

def linePlot(title, yaxis, xaxis, data):
    for i in range(len(data)):
        plt.plot(range(len(data[i][1])), data[i][1], color=COLORS[i])

    plt.legend([plt.Rectangle((0, 0), 1, 1, fc=COLORS[i]) for i in range(len(data))], [elem[0] for elem in data])

    plt.title(title)
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)

    plt.show()

def calcMedian(values):
    arr = array(values)

    return np.ma.median(arr)

def calcConfInt(values, interval=0.95):
    arr = array(values)

    return stats.norm.interval(interval, loc=arr.mean(), scale=arr.std()/math.sqrt(len(arr)))[1]-arr.mean()

def calcValues(values, interval=0.95):
    arr = array(values)

    return np.ma.median(arr), (arr.min(),arr.max()), stats.norm.interval(interval, loc=arr.mean(), scale=arr.std()/math.sqrt(len(arr)))

def combineLists(l1, l2):
    return [(l1[i],l2[i]) for i in range(len(l1))]

def maketbl(xlabels,ylabels,data):
    print('\\begin{tabularx}{\\textwidth}{|X|'),
    for xl in xlabels:
        print('c|'),
    print '}'
    print "\\hline \\rowcolor{slightgray}"
    for xl in xlabels:
        print("&" + xl),
    print "\\\\"

    for i in range(len(ylabels)):
        print "\\hline\\cellcolor{slightgray}"+ylabels[i]
        for dat in data[i]:
            print("&"+str(dat)),
        print "\\\\"

    print "\\hline\\end{tabularx}"

buildBarPlot("Speedup of different approaches on four 12 core machines", "Speedup",
    [
    ("Parsim",  (("LBTS", combineLists(results["sequential"], results["parLbts"])), ("DBS", combineLists(results["sequential"], results["parIBS"])))),
    ("Distributed Horizon", (("LBTS", combineLists(results["sequential"], results["comLbts"])), ("DBS", combineLists(results["sequential"], results["comIBS"]))))

    ], mapFunction=lambda ele:(float(ele[0]["runtime"])/float(NUMBER_OF_TTIS))/(float(ele[1]["runtime"])/float(NUMBER_OF_TTIS)), offsetx=0.9, filename="lte-eval-4machines.pdf", numberOffset=0.63, ylim=[0,48])

buildBarPlot("Speedup of different approaches on one 12 core machine", "Speedup",
    [
    ("Parsim", ((("LBTS", combineLists(results["sequential"], results["par12LPLbts"])), ("DBS", combineLists(results["sequential"], results["par12LPIBS"])),("Horizon", []),))),
    ("Horizon", (("LBTS", []), ("DBS", []), ("Horizon", combineLists(results["sequential"], results["horizon"])))),

    ], mapFunction=lambda ele:(float(ele[0]["runtime"])/float(NUMBER_OF_TTIS))/(float(ele[1]["runtime"])/float(NUMBER_OF_TTIS)), offsetx=0.3, filename="lte-eval-1machine.pdf", numberOffset=0.26)


exit(1)
buildBarPlot("Runtimes with Different Approaches", "Events/s",
    [
    ("Sequential",          (("C++", cppdata["sequential"]), ("YASiL Optimized(C++ Fading)", yasiloptifading["sequential"]), 
                                                             ("YASiL Optimized(WO C++ Fading)", yasilopti["sequential"]),
                                                             ("YASiL", yasilpure["sequential"]))),

    ("Horizon",             (("C++", cppdata["horizon"]),    ("YASiL Optimized(C++ Fading)", yasiloptifading["horizon"]), 
                                                             ("YASiL Optimized(WO C++ Fading)", yasilopti["horizon"]),
                                                             ("YASiL", []))),

    ("YASiL Unconditional", (("C++", cppdata["yasil1"]),     ("YASiL Optimized(C++ Fading)", yasiloptifading["yasil1"]), 
                                                             ("YASiL Optimized(WO C++ Fading)", yasilopti["yasil1"]),
                                                             ("YASiL", yasilpure["yasil1"]))),

    ("YASiL Conditional",   (("C++", cppdata["yasil2"]),     ("YASiL Optimized(C++ Fading)", yasiloptifading["yasil2"]), 
                                                             ("YASiL Optimized(WO C++ Fading)", yasilopti["yasil2"]),
                                                             ("YASiL", yasilpure["yasil2"]))),

    ("YASiL Unconditional IO", (("C++", cppdata["yasil1io"]),     ("YASiL Optimized(C++ Fading)", yasiloptifading["yasil1io"]), 
                                                             ("YASiL Optimized(WO C++ Fading)", yasilopti["yasil1io"]),
                                                             ("YASiL", yasilpure["yasil1io"]))),
    ], mapFunction=lambda ele:1000*ele["events"]/float(ele["runtime"]))

def secdiv(x, y):
    return float(x) / y if y != 0 else 0


overalll = reduce(lambda x,y: map(lambda w,z:w+z,x,y),[map(lambda x,y:x+y, ele["2ndYes"], ele["2ndUnd"]) for ele in yasilpure["yasil2"]])
overall = reduce(lambda x,y:x+y, overalll)

linePlot("Conditional Parallelization", "Probability", "Number of Workers", [
    ("Successful",map(lambda x,y: secdiv(y,x), overalll, reduce(lambda x,y: map(lambda w,z:w+z,x,y),[ele["2ndYes"] for ele in yasilpure["yasil2"]]))),
    ("Asked",[float(e)/overall for e in reduce(lambda x,y: map(lambda w,z:w+z,x,y),[map(lambda x,y:x+y, ele["2ndYes"], ele["2ndUnd"]) for ele in yasilpure["yasil2"]])])
    ])

overalll = reduce(lambda x,y: map(lambda w,z:w+z,x,y),[map(lambda x,y:x+y, ele["2ndYes"], ele["2ndUnd"]) for ele in yasilpure["yasil2io"]])
overall = reduce(lambda x,y:x+y, overalll)

linePlot("Conditional Parallelization IO RNG", "Probability", "Number of Workers", [
    ("Successful",map(lambda x,y: secdiv(y,x), overalll, reduce(lambda x,y: map(lambda w,z:w+z,x,y),[ele["2ndYes"] for ele in yasilpure["yasil2io"]]))),
    ("Asked",[float(e)/overall for e in reduce(lambda x,y: map(lambda w,z:w+z,x,y),[map(lambda x,y:x+y, ele["2ndYes"], ele["2ndUnd"]) for ele in yasilpure["yasil2io"]])])
    ])
