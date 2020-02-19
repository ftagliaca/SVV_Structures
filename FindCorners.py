from numpy import *


# Input : 2 axis the contain the data, the coordinates of the point you're tryin to locate

def FindCorners(gridx,gridy,x,y):

    # find the two values on each axis that are closest to the coordinate corresponding to the axis
    i = argpartition(abs(gridx-x),2)[:2] 
    j = argpartition(abs(gridy-y),2)[:2]

    return i,j

# Output : 2 arrays containing the 2 pairs of indeces
    

