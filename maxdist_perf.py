# Maximum distance algorithm in 2-D - Alain R. Marchand - February 2023
import os
import math
import numpy as np
import time

# Parameters
Size_list=[100_000]
Shape='square'
Repeat=1                                                                     # create a new set each time 
Precision=1+1.e-12
Step=1                                                                         # only 1 or a prime number, e.g. 911
Integers=False                                                              # Integers favor duplicate points 

# Debugging options
Test=False                                                                  # if True, perform unit tests
Test_set = None
# [(27, 1), (27, 10), (-7, 6), (3, 5)]                                # list of test points
if Test_set: Size_list, Repeat = [len(Test_set)], 1
Print_points=False

# Combinations
C_3_2=[(0,1,2), (2,0,1), (1,2,0)]                                                         # first 2 out of 3
C_3_2_plus=[(0,1,2,3), (2,0,1,3), (1,2,0,3)]                                         # first 2 out of 3 (ignore last)
C_4_2=[(0,1,2,3), (0,2,1,3), (0,3,1,2), (1,2,0,3), (1,3,0,2), (2,3,0,1)]     # first 2 out of 4
C_4_3=[(0,1,2,3), (3,0,1,2), (2,3,0,1), (1,2,3,0)]                                  # first 3 out of 4

#====================================== create_set
def create_set(size):
    """
    Draw 2D points from a Gaussian distribution
    Parameters: size, integer
    Returns: points, a list of tuples of 2 float
    """
# uniform square
    if Shape=='square':
        x=list(np.random.random(size)*10) 
        y=list(np.random.random(size)*10)         

# uniform disk
    elif Shape=='circle':
        x=list(np.random.random(2*size)*10) 
        y=list(np.random.random(2*size)*10)                  # start with a square
        xy=[(a-5,b-5) for a,b in zip(x,y) if (a-5)**2+(b-5)**2<=25]
        x, y = list(zip(*xy[:size]))                                        
    
    if Integers: x, y = list(map(int, x)), list(map(int, y))

    points=list(zip(x,y))
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
    return (A[0]-B[0])**2 + (A[1]-B[1])**2
    
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
    Returns points
    """
    points[position], points[size] = points[size], points[position]
    return points[:-1], size                        

#====================================== point_outside
def point_outside(sample, position, O, r_squared):
    """
    look for a point outside circle
    cyclic search in sample
    stop if found
    Parameters: sample, a list of 3 points, tuples of 2 float
                  position, an integer between 0 and len(sample)-1
                  O, a point, center of the circle
                  r_squared, a float, square of circle radius
    Returns: E, point if found, else None
                  new position, an integer
    """
    start=position
    while True:
        E=sample[position]
        position=inc_position(1, position, Step, len(sample))
        if not inside(E, O, r_squared):
            return E, position                                   # a point outside
        if position==start:
            return None, position                            # no point outside
        
#====================================== circle
def circle(A, B):
    """
    calculate circles based on 2 points
    Parameters: A, B, two points to use as a diameter
                  remain: a set of remaining points to test for inclusion in the circle
    Returns: O, a point, tuple of 2 float
                  r_squared, a float
                  [A, B], a list of 2 points forming a diameter of the circle
                  True if all points are in the circle, False otherwise
    """
    # calculate circle
    O = (A[0]+B[0])/2, (A[1]+B[1])/2                          # middle of AB
    r_squared=d_squared(A, O)
    return O, r_squared

#====================================== unit_tests
def unit_tests():
    """
    test each module
    """
    global steps, distances, position
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

    assert(obtuse(sample[:3])==True)                        # obtuse
    assert(obtuse([(0,3), (8,3), (4,-5)])==False)                        

    assert(point_outside(sample, 0, (5,5), 25)==((9,9), 0))  # point_outside & inc_position
    assert(point_outside(sample, 3, (5,5), 32)==(None, 3))

    assert(extremes(sample)==[(2,1), (9,9), (2,1), (9,9)])

    assert(circle2((1,2), (1,4), [(1,3)])==((1,3), 1, True))  # circle2
    assert(circle2((1,2), (1,3), [(1,4)])==(None, None, False)) 

    assert(circle3(sample[:3])==((-0.5, 4.5), 18.5))     # circle3
    assert(circle3([(1,2), (1,3), (1,4)])==((1,3), 1))
    assert(circle3([(0,3), (8,3), (4,-5)])==((4,0), 25))
    assert(circle3([(3,0), (3,8), (-5,4)])==((0,4), 25))

    assert(minimal2(sample[1:])==((5.5, 5.0), 28.25, [(2, 1), (9, 9)], True))   # minimal2 (4 pts)
    assert(minimal2(sample[:3])==((2.5, 4.0), 9.25, [(3, 7), (2, 1)], True))     # minimal2 (3 pts)
    assert(minimal2(sample[1:4])==(None, None, [None, None], False))     # minimal2 (3 pts)

    assert(minimal3([(0,3), (8,3), (4,-5)])==((4.0, 0.0), 25.0, [(0, 3), (8, 3), (4, -5)])) # minimal3 (3 pts)
    assert(minimal3([(0,3), (8,3), (4,-5), (6,2)])==((4.0, 0.0), 25.0, [(0, 3), (8, 3), (4, -5)])) # minimal3 (4 pts)

    subset, O, r_squared = minicircle([(0,3), (8,3), (4,-5), (1,1), (6,2)], [(0,3), (8,3), (4,-5), (6,2)], 0)     # minicircle (4 pts) 
    assert((set(subset), O,r_squared)==(set([(0,3), (8,3), (4,-5)]), (4.0, 0.0), 25.0))                              # order is unreliable        

    print("\nTests OK")    
    os._exit(1)
    
#====================================== main program
"""
Draw sample of 2-D points and determine their minimal enclosing circle
evaluate performance for each sample
"""
if Test: unit_tests()

if Step>1:
    print("Warning: with a Step>1 you can avoid randomization, but it is your responsibility")
    print("to verify that Step is a prime number in order to avoid avoid subsampling")
    print("Examples of prime numbers: 53, 97, 179, 271, 367, 557, 719, 911, 1777")
    print("Your Step:", Step, "\n")

# compute increasing size list
size=100
Size_list=[]
while size<=1_000_000:
    Size_list+=[size]
    size=int(size*1.05)
print(Shape, "\nSet_size\tDist./N\tSteps\tTime(ms)")

for index, Set_size in enumerate(Size_list):   
    for rep in range(Repeat):

        # performance stats
        steps=0                                                              # number of subsets examined
        distances=0                                                        # number of distances measured
        
        # initialize set of points
        points=create_set(Set_size)
        if Test: points=Test

        # compute minimum circle
        position=0

##        points=points[:8]                                        # debug
##        Set_size=len(points)
##        print(points, '\n')
        
        diameter=points[:2]
        size=Set_size
        steps=1
        
        # start time
        t0=time.perf_counter()

##################################################
        
        while True:
            O,r_squared=circle(*diameter)    # avoid repeating this to optimize
##            print(position, diameter, r_squared)
            
            C, position =point_outside(points, position, O, r_squared)
            if not C: break                                                 # AB is solution
##            print("found", C, "at position", (position-1)%size)
            D, position =point_outside(points, position, C, 4*r_squared)
            if not D:
                position=(position-1)%size
                points, size = eliminate(points, position, size-1)  # eliminate C
                if position==size: position=0
##                print("eliminating", C, "at position", position)
##                print(points, '\n')
                continue
##            print("adding", D)
            diameter=[C, D]
            steps+=1

##################################################

        # end time
        t1=time.perf_counter()

        # print stats
        print(Set_size, f"\t{distances/Set_size:3.3}\t{steps}\t{(t1-t0)*1000:9.2f}")

        # just a verification
        if point_outside(points, 0, O, r_squared)[0]:                     
            print("\n*** Not solved ! ***\n")
            os._exit(1)

print("\n")


