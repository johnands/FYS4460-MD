import numpy as np
import scipy.fftpack
import matplotlib.pyplot as plt

def extract(filename):

    infile = open(filename, 'r')

    numberOfBins = int(infile.readline())
    infile.readline()							# skip comment line
			
    # initialize list
    bins = []
    radialDistribution = []

    # extract data from file
    i = 0
    for line in infile:
        words = line.split() 
        if i < numberOfBins:
		    bins.append(float(words[0]))
        radialDistribution.append(float(words[1]))
        i += 1

    infile.close()

    return np.array(radialDistribution), np.array(bins), numberOfBins


def smoothing(x, y, N):

   w = scipy.fftpack.rfft(y)
   f = scipy.fftpack.rfftfreq(N, x[1]-x[0])
   spectrum = w**2

   cutoff_idx = spectrum < (spectrum.max()/5)
   w2 = w.copy()
   w2[cutoff_idx] = 0

   y2 = scipy.fftpack.irfft(w2)

   return y2



def treatData():
    radialDistribution, bins, numberOfBins = extract('radialDistribution.txt')

    numberOfTimeSteps = len(radialDistribution) / len(bins)


    # find time average of radial distribution
    cumulativeRadialDistribution = np.zeros(numberOfBins)
    for t in xrange(numberOfTimeSteps):
        for i in xrange(numberOfBins):
            cumulativeRadialDistribution[i] += radialDistribution[numberOfBins*t+i]

    # average over time
    radialDistTimeAverage = cumulativeRadialDistribution / numberOfTimeSteps
    
    # smooth function
    radialDistTimeAverageSmooth = smoothing(bins, radialDistTimeAverage, numberOfBins)

    plt.plot(bins, radialDistTimeAverageSmooth, 'b-', bins, np.zeros(numberOfBins)+1, 'r-')
    plt.show()
    

# ----- main -----

treatData()

