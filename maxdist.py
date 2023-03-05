# Maximum distance algorithm in 2-D - Alain R. Marchand - February 2023
import os
import math
import numpy as np
import time

# Parameters
Size_params= 100, 20_000, 1.05                                    # min, max, geometric progression
Dimension=3
Shape='square'
Repeat=1                                                                     # create a new set each time 
Precision=1+1.e-12
Step=1                                                                         # only 1 or a prime number, e.g. 911
Integers=False                                                              # Integers favor duplicate points 

# Debugging options
Verify=False
Test=False                                                                  # if True, perform unit tests
Test_set = None
# [(27, 1), (27, 10), (-7, 6), (3, 5)]                                # list of test points
if Test_set: Size_list, Repeat = [len(Test_set)], 1
Print_points=False

#====================================== create_set
def create_set(size):
    """
    Draw 2D points from a Gaussian distribution
    Parameters: size, integer
    Returns: points, a list of tuples of 2 float
    """
# uniform square
    shape=(size, Dimension)
    array=np.zeros(shape, dtype=float)
    if Shape=='square':
        for d in range(Dimension):
            array[:,d]=list(np.random.random(size)*10) 

    points=[]
    for i in range(size):
        points+=[tuple(array[i,:])]

    if Print_points: print_points(points)

    return points

#====================================== print_points
def print_points(points):
    """
    Print points and stop
    """
    for x,y in points:
        print(x, "\t", y)
    os._exit(1)

#====================================== d_squared
def d_squared(A, B):
    """
    return squared distance between points A and B
    increment global count distances
    """
    global distances
    distances+=1                                                                       # counts as 1 operation
    d=0
    for i in range(Dimension):
        d+=(A[i]-B[i])**2
    return d
    
#====================================== inc_position
def inc_position(n, position, step, size):
    """
    add n steps to position, modulo size
    """                 
    return (position+n*step) % size  

#====================================== inside
def inside(A, O, r_squared):
    """
    check whether point A is in circle (O, r_squared)
    Returns True if inside, False if outside
    """                 
    return d_squared(A, O)<=r_squared*Precision                  # allow for rounding errors                        

#====================================== eliminate
def eliminate(points, position, size):
    """
    eliminate a point by placing it at the end
    Parameters: points, a list of points, tuples of 2 float
                  position, an integer between 0 and size-1
                  size: max size to explore
    Returns: points
                  size: new size
    """
    points[position], points[size] = points[size], points[position]
    return points, size                        

#====================================== point_outside
def point_outside(sample, size, position, O, r_squared):
    """
    look for a point outside circle
    cyclic search in sample
    stop if found
    Parameters: sample, a list of 3 points, tuples of 2 float
                  size: max size to explore
                  position, an integer between 0 and size-1
                  O, a point, center of the circle
                  r_squared, a float, square of circle radius
    Returns: E, point if found, else None
                  position: new position, an integer
    """
    start=position
    while True:
        E=sample[position]
        position=inc_position(1, position, Step, size)
        if not inside(E, O, r_squared):
            return E, position                                   # a point outside
        if position==start:
            return None, position                            # no point outside
        
#====================================== circle
def circle(A, B):
    """
    calculate circle based on 2 points
    Parameters: A, B, two points to use as a diameter
    Returns: O, a point, tuple of 2 float
                  r_squared, a float
    """
    # calculate circle
    O=[0]*Dimension
    for i in range(Dimension):
        O[i] = (A[i]+B[i])/2                                              # middle of AB
    O=tuple(O)
    r_squared=d_squared(A, O)
    return O, r_squared

#====================================== unit_tests
def unit_tests():
    """
    test each module
    """
    global Dimension, steps, distances, position
    Dimension=2
    steps=0                                                              # number of subsets examined
    distances=0                                                        # number of distances measured
    position=0                                                          # index to points
    sample=[(3,2), (3,7), (2,1), (8,2), (9,9)]
    
    assert(d_squared((0,0), (3,4))==25)                   # d_squared
    assert(distances==1)                                          # distances
    assert(inc_position(2, position, 7, 11) == 3)       # inc_position
    
    assert(inside((0,0), (3,4), 25))                              # inside
    assert(inside((3,4), (0,0), 26))                             
    assert(not inside((0,3), (4,0), 24))

    assert(point_outside(sample, 5, 0, (5,5), 25)==((9,9), 0))  # point_outside & inc_position
    assert(point_outside(sample, 5, 3, (5,5), 32)==(None, 3))

    assert(circle((1,2), (1,4))==((1,3), 1) )                 # circle

    print("\nTests OK")    
    os._exit(1)
    
#====================================== main program
"""
Draw sample of 2-D points and determine the most distant points
evaluate performance for each sample
"""
if Test: unit_tests()

if Step>1:
    print("Warning: with a Step>1 you can avoid randomization, but it is your responsibility")
    print("to verify that Step is a prime number in order to avoid avoid subsampling")
    print("Examples of prime numbers: 53, 97, 179, 271, 367, 557, 719, 911, 1777")
    print("Your Step:", Step, "\n")

# compute increasing size list
if not Test_set:
    size, max_size, coef=Size_params
    Size_list=[]
    while size<=max_size:
        Size_list+=[size]
        size=int(size*coef)
    print(Shape, "\tdim", Dimension, "\nSet_size\tDist./N\tSteps\tTime(ms)\tt/pt(µs)")
else:
    Size_list=[len(Test_set)]

for index, Set_size in enumerate(Size_list):   
    for rep in range(Repeat):

        # performance stats
        steps=0                                                              # number of subsets examined
        distances=0                                                        # number of distances measured
        
        # initialize set of points
        points=create_set(Set_size)
        if Test: points=Test

        # compute minimum circle
        position=2

        diameter=points[:2]
        size=Set_size
        steps=1
        
        # start time
        t0=time.perf_counter()

##################################################
        
        while True:
            O,r_squared=circle(*diameter)    
            
            C, position =point_outside(points, size, position, O, r_squared)
            if not C: break                                                 # AB is solution
            D, position =point_outside(points, size, position, C, 4*r_squared)
            if not D:
                position=(position-1)%size
                points, size = eliminate(points, position, size-1)  # eliminate C
                if position==size: position=0
                continue
            diameter=[C, D]
            steps+=1

##################################################

        # end time
        t1=time.perf_counter()
        rt=(t1-t0)*1000
        
        # print stats
        print(f"{Set_size}\t{distances/Set_size:3.3}\t{steps}\t{rt:9.2f}\t{rt*1000/Set_size:3.3}")

        # just a verification
        if Verify:
            r_squared=d_squared(*diameter)
            for A in points:
                for B in points:
                    d=d_squared(A, B)
                    if d>r_squared:                     
                        print("\n*** Not solved ! ***\n")
                        os._exit(1)

print("\n")

