import matplotlib.pyplot as plt
import pickle
import numpy

def readIn(ins):
    fin = open(ins+'.pkl','rb')
    maxes = pickle.load(fin)
    maxs = []
    maxavg = []
    nummax = 0
    for i in range(len(maxes[0])):
        maxs.append(numpy.max(maxes[0][i]))
        maxavg.extend(maxes[0][i])
        nummax += len(maxes[0][i])
    fin = open(ins+'_v2.pkl','rb')
    maxes = pickle.load(fin)
    for i in range(len(maxes[0])):
        maxs.append(numpy.max(maxes[0][i]))
        maxavg.extend(maxes[0][i])
        nummax += len(maxes[0][i])
    theMean = numpy.mean(maxavg)
    print theMean
    print numpy.sum((maxavg-theMean)**2)/float(nummax)
    print numpy.mean((maxavg-theMean)**2)
    print numpy.sum((maxavg-theMean)**3)/float(nummax)
    print numpy.mean((maxavg-theMean)**3)
    print numpy.sum((maxavg-theMean)**4)/float(nummax)
    print numpy.mean((maxavg-theMean)**4)
    print numpy.sum((maxavg-theMean)**5)/float(nummax)
    print numpy.mean((maxavg-theMean)**5)
    print nummax