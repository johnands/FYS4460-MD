import numpy as np
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



def treatData():
    radialDistribution, bins, numberOfBins = extract('radialDistribution.txt')
    print numberOfBins

    numberOfTimeSteps = len(radialDistribution) / len(bins)
    binSize = bins[1] - bins[0]


    # find time average of radial distribution
    cumulativeRadialDistribution = np.zeros(numberOfBins)
    for t in xrange(numberOfTimeSteps):
        for i in xrange(numberOfBins):
            cumulativeRadialDistribution[i] += radialDistribution[numberOfBins*t+i]

    # average over time
    radialDistTimeAverage = cumulativeRadialDistribution / numberOfTimeSteps

    plt.plot(bins, radialDistTimeAverage, 'b-', bins, np.zeros(numberOfBins)+1, 'r-')
    plt.show()
    

# ----- main -----

treatData()

