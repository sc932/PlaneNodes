import random, numpy
import matplotlib.pyplot as plt
import time
import scipy.ndimage
import math
import pickle

h = 0.5

def Run1(waves,times,reals):
    #f = open('./Run1_' + str(waves) + 'waves_' + str(times) + 'times_' + str(reals) + '.dat','w')
    output = open('./Run1_' + str(waves) + 'waves_' + str(times) + 'times_' + str(reals) + '_v2.pkl', 'wb')
    bigbig = []
    totaltime = 0
    #f.write('Run1_' + str(waves) + 'waves_' + str(times) + 'times_' + str(reals) + '\n')
    #f.write(str(h) + 'h, 1024g\n')
    for i in range(times):
        numMax, avgMax, biggest, avgH, timeTaken, numWaves, gridSize, hf, real = RandomNodes2D(waves,512,reals)
        bigbig.append(biggest)
        #f.write(str(numMax) + ',' + str(avgMax) + ',' + str(biggest) + ',' + str(avgH) + ',' + str(timeTaken) + '\n')
        print i, timeTaken
        totaltime += timeTaken
    #print bigbig
    pickle.dump([bigbig],output)
    output.close()
    #f.write(str(totaltime) + 's total')
    #f.close()
                
def RandomNodes2D(numWaves, gridSize, real):
    theMat = numpy.zeros((gridSize,gridSize))
    t1 = time.time()
    t0 = t1
    if real:
        if numWaves == 1000:
            theMat = numpy.fromfunction(thousandWaves, (gridSize,gridSize))
        elif numWaves == 2000:
            theMat = numpy.fromfunction(twoThousandWaves, (gridSize,gridSize))
        elif numWaves == 4000:
            theMat = numpy.fromfunction(fourThousandWaves, (gridSize,gridSize))
        elif numWaves == 12:
            theMat = numpy.fromfunction(dozenWaves, (gridSize,gridSize))
        else:
            print "invalid number of waves, try 1000, 2000 or 4000"
        #print time.time()-t1
    else:
        if numWaves == 1000:
            theMat = numpy.fromfunction(thousandWavesi, (gridSize,gridSize))
        elif numWaves == 2000:
            theMat = numpy.fromfunction(twoThousandWavesi, (gridSize,gridSize))
        elif numWaves == 4000:
            theMat = numpy.fromfunction(fourThousandWavesi, (gridSize,gridSize))
        else:
            print "invalid number of waves, try 1000, 2000 or 4000"
        #print time.time()-t1
    t1 = time.time()
    numMax = 0
    maxH = 0.0
    avgH = 0.0
    biggest = []
    for x in range(1,gridSize-1):
        for y in range(1,gridSize-1):
            avgH += theMat[x][y]/(gridSize*gridSize)
            if theMat[x][y] > max(theMat[x-1][y],theMat[x-1][y-1],theMat[x-1][y+1],theMat[x][y-1],theMat[x+1][y],theMat[x+1][y-1],theMat[x+1][y+1],theMat[x][y+1]):
                #print x,y
                numMax += 1
                maxH += theMat[x][y]
                biggest.append(theMat[x][y])
    #print numMax
    #print maxH/numMax
    #print biggest
    #print avgH
    #print time.time()-t1
    #plt.matshow(theMat)
    return numMax, maxH/numMax, biggest, avgH, time.time()-t0, numWaves, gridSize, h, real
    
def thousandWaves(x,y):
    ans = 0.0
    for i in range(1000):
        direction = random.uniform(0,2*numpy.pi)
        phase = random.uniform(0,2*numpy.pi)
        ans += numpy.sin(numpy.sin(direction)*x*h+numpy.cos(direction)*y*h+phase)
    return ans**2
    
def thousandWaves2(x,y):
    numWaves = 2
    direction = []
    phase = []
    for i in range(numWaves):
        direction.append(random.uniform(0,2*numpy.pi))
        phase.append(random.uniform(0,2*numpy.pi))
    direction = numpy.array(direction)
    phase = numpy.array(phase)
    phase = numpy.reshape(phase,(numWaves,1,1))
    phase = phase*numpy.ones(numpy.shape(x))
    direction = numpy.reshape(direction,(1,1,numWaves))
    sDir = numpy.sin(direction)
    cDir = numpy.cos(direction)
    x = numpy.reshape(x,(numpy.shape(x)[0],numpy.shape(x)[1],1))
    y = numpy.reshape(y,(numpy.shape(y)[0],numpy.shape(y)[1],1))
    ans = sum(numpy.sin(sDir*x*h+cDir*y*h+phase),axis=0)
    return ans*ans
    
def twoThousandWaves(x,y):
    ans = 0.0
    direction = numpy.random.uniform(0,2*numpy.pi,2000)
    phase = numpy.random.uniform(0,2*numpy.pi,2000)
    for i in range(2000):
        #direction = random.uniform(0,2*numpy.pi)
        #phase = random.uniform(0,2*numpy.pi)
        ans += numpy.sin(numpy.sin(direction[i])*x*h+numpy.cos(direction[i])*y*h+phase[i])
    return ans**2
    
def twoThousandWaves2(x,y):
    direction = numpy.random.uniform(0,2*numpy.pi,2000)
    phaseA = numpy.random.uniform(0,2*numpy.pi,2000)
    phaseA = numpy.reshape(phaseA,(2000,1,1))
    direction = numpy.reshape(direction,(1,1,2000))
    directionS = numpy.sin(numpy.reshape(direction,(2000,1,1)))
    directionC = numpy.cos(numpy.reshape(direction,(2000,1,1)))
    xA = numpy.reshape(x,(1,int(math.sqrt(numpy.size(x))),int(math.sqrt(numpy.size(x)))))
    yA = numpy.reshape(y,(1,int(math.sqrt(numpy.size(y))),int(math.sqrt(numpy.size(y)))))
    ans1 = numpy.sin(directionS*xA*h+directionC*yA*h+phaseA*numpy.reshape(numpy.ones((int(math.sqrt(numpy.size(x))),int(math.sqrt(numpy.size(x))))),(1,int(math.sqrt(numpy.size(x))),int(math.sqrt(numpy.size(x))))))
    ansSum = numpy.sum(ans1,axis=0)
    return ansSum**2
    return ans*ans
    
def fourThousandWaves(x,y):
    ans = 0.0
    for i in range(4000):
        direction = random.uniform(0,2*numpy.pi)
        phase = random.uniform(0,2*numpy.pi)
        ans += numpy.sin(numpy.sin(direction)*x*h+numpy.cos(direction)*y*h+phase)
    return ans**2
    
def thousandWavesi(x,y):
    ansR = 0.0
    ansI = 0.0
    for i in range(1000):
        directionR = random.uniform(0,2*numpy.pi)
        phaseR = random.uniform(0,2*numpy.pi)
        ansR += numpy.sin(numpy.sin(directionR)*x*h+numpy.cos(directionR)*y*h+phaseR)
        directionI = random.uniform(0,2*numpy.pi)
        phaseI = random.uniform(0,2*numpy.pi)
        ansI += numpy.sin(numpy.sin(directionI)*x*h+numpy.cos(directionI)*y*h+phaseI)
    return ansR**2 + ansI**2
    
def twoThousandWavesi(x,y):
    ansR = 0.0
    ansI = 0.0
    for i in range(2000):
        directionR = random.uniform(0,2*numpy.pi)
        phaseR = random.uniform(0,2*numpy.pi)
        ansR += numpy.sin(numpy.sin(directionR)*x*h+numpy.cos(directionR)*y*h+phaseR)
        directionI = random.uniform(0,2*numpy.pi)
        phaseI = random.uniform(0,2*numpy.pi)
        ansI += numpy.sin(numpy.sin(directionI)*x*h+numpy.cos(directionI)*y*h+phaseI)
    return ansR**2 + ansI**2
    
def fourThousandWavesi(x,y):
    ansR = 0.0
    ansI = 0.0
    for i in range(4000):
        directionR = random.uniform(0,2*numpy.pi)
        phaseR = random.uniform(0,2*numpy.pi)
        ansR += numpy.sin(numpy.sin(directionR)*x*h+numpy.cos(directionR)*y*h+phaseR)
        directionI = random.uniform(0,2*numpy.pi)
        phaseI = random.uniform(0,2*numpy.pi)
        ansI += numpy.sin(numpy.sin(directionI)*x*h+numpy.cos(directionI)*y*h+phaseI)
    return ansR**2 + ansI**2
    
def singleWave(x,y):
    print x
    print y
    phase = random.uniform(0,2*numpy.pi)
    direction = random.uniform(0,2*numpy.pi)
    return numpy.sin(numpy.sin(direction)*x*h+numpy.cos(direction)*y*h+phase)

def dozenWaves(x,y):
    phase = []
    direction = []
    for i in range(1000):
        phase.append(random.uniform(0,2*numpy.pi))
        direction.append(random.uniform(0,2*numpy.pi))
    phaseA = numpy.reshape(numpy.array(phase),(len(phase),1,1))
    directionS = numpy.sin(numpy.reshape(numpy.array(direction),(len(direction),1,1)))
    directionC = numpy.cos(numpy.reshape(numpy.array(direction),(len(direction),1,1)))
    xA = numpy.reshape(x,(1,int(math.sqrt(numpy.size(x))),int(math.sqrt(numpy.size(x)))))
    yA = numpy.reshape(y,(1,int(math.sqrt(numpy.size(y))),int(math.sqrt(numpy.size(y)))))
    #print numpy.shape(phaseA)
    #print numpy.shape(directionS)
    #print numpy.shape(directionC)
    #print numpy.shape(xA)
    #print numpy.shape(yA)
    #print numpy.shape(directionS*xA*h)
    #print numpy.shape(directionC*yA*h)
    #print numpy.shape(phaseA*numpy.reshape(numpy.ones((int(math.sqrt(numpy.size(x))),int(math.sqrt(numpy.size(x))))),(1,int(math.sqrt(numpy.size(x))),int(math.sqrt(numpy.size(x))))))
    ans1 = numpy.sin(directionS*xA*h+directionC*yA*h+phaseA*numpy.reshape(numpy.ones((int(math.sqrt(numpy.size(x))),int(math.sqrt(numpy.size(x))))),(1,int(math.sqrt(numpy.size(x))),int(math.sqrt(numpy.size(x))))))
    print numpy.shape(ans1)
    ansSum = numpy.sum(ans1,axis=0)
    return ansSum**2