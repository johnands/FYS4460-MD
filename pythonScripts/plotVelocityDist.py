import numpy as np
import matplotlib.pyplot as plt

def extract(filename):

    infile = open(filename, 'r')

    n = int(infile.readline()) 		# number of atoms
    infile.readline()				# skip comment line
			
    # initialize lists
    vx = []; vy = []; vz = []

    # extract velocities from file
    i = 0
    for line in infile:
        if i >= n:
            print i
            break
        words = line.split() 
        vx.append(float(words[4]))
        vy.append(float(words[5]))
        vz.append(float(words[6]))
        i += 1


    infile.close()

    return vx, vy, vz, n



def treat():
    vx, vy, vz, n = extract('movie.xyz')

    # convert to numpy arrays
    vx = np.array(vx); vy = np.array(vy); vz = np.array(vz)

    # magnitude of velocites
    v = np.sqrt(vx**2 + vy**2 + vz**2)
    
    # compute mean
    mean_vx = sum(vx)/n
    mean_vy = sum(vy)/n
    mean_vz = sum(vz)/n

    print "Mean vx = ", mean_vx
    print "Mean vy = ", mean_vy 
    print "Mean vz = ", mean_vz

    # compute standard deviation
    sigma_x = np.sqrt(sum((vx-mean_vx)**2)/n)
    sigma_y = np.sqrt(sum((vy-mean_vy)**2)/n)
    sigma_z = np.sqrt(sum((vz-mean_vz)**2)/n)
    print "St.dev. vx = ", sigma_x
    print "St.dev. vy = ", sigma_y
    print "St.dev. vz = ", sigma_z

    
    plt.hist(vx, bins=100)
    plt.show()

    plt.hist(v, bins=100)
    plt.show()
    

# main

treat()
